##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Caliper placement shared by the '-timed' transformations.

One table, CAPTURED_REGIONS, names the seven loops the Kokkos transformation
captures. Every timed build puts a caliper on those seven and on nothing else,
so a row in one build's timer.txt covers the same executions as the row with
the same name in another's. That is the whole purpose of this module: the
comparison is only a comparison if both sides bracket the same work.

WHY THIS EXISTS RATHER THAN psyclone_tools.profile_loops

psyclone_tools already places PSyData calipers, and this module deliberately
does not call it. Two reasons, and only the second is about Kokkos.

1. profile_loops names a region by position:

       region_name = invoke_name + ":" + k_name + "_k" + str(count)

   where 'count' indexes into the invoke's coded kernels and increments only
   for loops that were profiled. A region name carrying an ordinal is not
   stable across two builds that generate different numbers of loops -- which
   is precisely what a Kokkos build is. Capturing one loop as a Kokkos region
   removes it from coded_kernels(), so every later kernel in that invoke would
   be renamed, and a Fortran row and a Kokkos row describing the same kernel
   would not join.

2. profile_loops iterates loops and skips anything without coded_kernels(). A
   captured region is a Call, not a Loop, so it would be skipped in silence --
   the build would succeed and report a time for every region except the ones
   being measured.

psyclone_tools.py is lfric_core infrastructure shared by every application, so
it is read and not edited: teaching it about Kokkos regions would put one
experiment's knowledge into everybody's build. The constraint it enforces is
right and is reproduced here rather than relaxed -- a caliper landing inside an
OpenMP region is refused below with the same wording, so that a wrong ordering
fails the build instead of mismeasuring it.

WHY ONLY THE CAPTURED CALL SITES ARE TIMED, AND NOT EVERY LOOP

Two measured facts, neither of them a matter of taste.

*The timer table is not big enough.* gungho_model runs 298 distinct coded
kernels, counted from the generated PSy layer of a built 'minimum':

    grep -rhoiE '^ *CALL +[a-z0-9_]+_code *\\(' --include='*_psy.f90' . \\
        | sed -E 's/^ *[Cc][Aa][Ll][Ll] +//; s/ *\\($//' | sort -u | wc -l

timer_mod holds 300 -- 'integer(i_def), parameter :: num_subs = 300' -- and
LFRic's own timing_mod claims one before any of these are reached. A caliper on
every coded-kernel loop runs the table out and stops the model at

    call log_event( "Run out of timers, increase num_subs", LOG_LEVEL_ERROR )

which is fatal, not a warning. Raising num_subs would be an lfric_core change
made for a measurement harness, and it would not fix the second problem:
timer() finds a region by scanning the table linearly and comparing 128-byte
names,

    do k = 1, num_tim_in_use
       if( lowname == routine_name(k) ) exit
    end do

so a full table costs some hundreds of string comparisons on every caliper
entry and every caliper exit. Against a region whose measured minimum is 36
microseconds that is not negligible, and it falls on the baseline and the
Kokkos build alike -- inflating both sides of the comparison this stage exists
to make.

*A kernel is not a call site.* Four of the seven captured kernels are invoked
from more than one algorithm:

    fv_difference_z      ffsl_vert_alg_mod, flux_precomputations_mod
    inject_wt_to_sh_w3   end_of_transport_step_alg_mod, ffsl_control_alg_mod,
                         mol_consistent_alg_mod,
                         transport_rho_times_field_alg_mod
    moist_dyn_gas        init_gungho_lbcs_alg_mod, moist_dyn_factors_alg_mod
    sample_w3_to_wtheta  map_fd_to_prognostics_alg_mod,
                         physics_mappings_alg_mod, si_operators_alg_mod

The Kokkos build captures one call site of each. A caliper placed by kernel name
alone would, in that build, accumulate the captured Kokkos executions and the
uncaptured Fortran ones into a single row -- so the row would not measure the
Kokkos region, and the difference against the baseline would be diluted by an
unknown amount. Selecting by (module, invoke, kernel) is what makes the row
mean one thing.

HOW A REGION IS NAMED, AND WHY IT JOINS ACROSS BUILDS

    ('region', <kernel name with a trailing '_code' removed>)

The LFRic timer concatenates the pair with a colon, so apply_helmholtz_operator
appears in timer.txt as 'region:apply_helmholtz_operator' in every build. The
colon also keeps these rows distinct from LFRic's own, which are bare
subroutine names. The kernel name alone is enough to be unique because
CAPTURED_REGIONS names each kernel once; that is checked below rather than
assumed.

The name comes from the *Fortran* kernel name in both builds, which is what
makes it join. The two builds do not agree on the name of the thing being run:
the Fortran PSy layer calls tri_solve_code and the Kokkos PSy layer calls
tri_solve_r_single_kokkos, because LFRicKokkosTrans has to select one
implementation of a kind-polymorphic interface where Fortran can defer it.
Deriving a common key by stripping precision infixes off the Kokkos name would
be guesswork about a naming convention. It is not needed: the table already
states the Fortran kernel, so the two sides agree by construction.

WHERE THIS FILE SITS

At the root of optimisation/, above the transformation directories, because all
three timed transformations import it and the region names have to agree
between them. Three copies of this table could drift, and a drift would show up
as a table that silently fails to join rather than as a build error.

The scripts reach it by walking their own parents for the directory holding
this file, rather than by a fixed number of '.parent's, so that a copy placed
inside a transformation tree is found first. That matters if a timed tree is
ever bound out-of-tree through LFRIC_OPTIMISATION_PATH, which binds the
transformation directory alone and would leave this one outside the container.
'''
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.nodes import (
    Call,
    Container,
    OMPDoDirective,
    OMPParallelDirective,
    OMPParallelDoDirective,
    ProfileNode,
)
from psyclone.psyir.transformations import ProfileTrans
from psyclone.transformations import TransformationError

#: The module half of every region identifier this module places. The LFRic
#: timer joins it to the region half with a colon, so every row this harness
#: produces is greppable out of a timer.txt that also holds LFRic's own.
REGION_MODULE = 'region'

#: The call sites the Kokkos transformation captures, as
#: (PSy-layer module, invoke, coded kernel). Each row is the TARGET_INVOKE and
#: TARGET_KERNEL of one script under kokkos/psykal/, with the PSy module the
#: algorithm generates; the Kokkos scripts check themselves against it, so
#: capturing an eighth region without adding it here fails that build rather
#: than producing a Kokkos time with no baseline beside it.
#:
#: Every row must name a distinct kernel, since the kernel alone is the region
#: name. Asserted by _region_names() below.
CAPTURED_REGIONS = (
    ('pressure_operator_alg_mod_psy',
     'invoke_apply_h',
     'apply_helmholtz_operator_code'),
    ('ffsl_vert_alg_mod_psy',
     'invoke_0_fv_difference_z_kernel_type',
     'fv_difference_z_code'),
    ('sci_geometric_constants_mod_psy',
     'invoke_12_get_dz_w3_kernel_type',
     'get_dz_w3_code'),
    ('ffsl_control_alg_mod_psy',
     'invoke_1_inject_wt_to_sh_w3_kernel_type',
     'inject_wt_to_sh_w3_code'),
    ('moist_dyn_factors_alg_mod_psy',
     'invoke_2',
     'moist_dyn_gas_code'),
    ('si_operators_alg_mod_psy',
     'invoke_6_sample_w3_to_wtheta_kernel_type',
     'sample_w3_to_wtheta_code'),
    ('pressure_precon_alg_mod_psy',
     'invoke_0_tri_solve_kernel_type',
     'tri_solve_code'),
)

#: Loop types that carry a coded kernel but are the inner half of a colouring,
#: so timing them would time a fraction of the work and would not correspond to
#: anything in an uncoloured build. Copied from psyclone_tools.profile_loops
#: because the two must agree about what a coloured loop looks like.
INNER_COLOUR_LOOPS = ('cells_in_colour', 'tiles_in_colour', 'cells_in_tile')


def region_name(kernel):
    '''
    Turns a coded kernel's name into the region half of its identifier.

    :param str kernel: the name of the coded kernel, e.g. 'tri_solve_code'.

    :returns: the region name, e.g. 'tri_solve'.
    :rtype: str

    '''
    name = kernel.lower()
    suffix = '_code'
    if name.endswith(suffix):
        name = name[:-len(suffix)]
    return name


def _region_names():
    '''
    The region names CAPTURED_REGIONS produces, checked to be distinct.

    :returns: one region name per row of the table.
    :rtype: list[str]

    :raises RuntimeError: if two rows would produce the same region name, in
        which case their times would be accumulated into one row of timer.txt
        and neither would mean anything.

    '''
    names = [region_name(kernel) for _, _, kernel in CAPTURED_REGIONS]
    if len(set(names)) != len(names):
        raise RuntimeError(
            f"CAPTURED_REGIONS names a kernel twice: {sorted(names)}. The "
            "kernel is the region name, so two rows would share one timer.")
    return names


def _options(kernel):
    '''
    Builds the ProfileTrans options that name a region after its kernel.

    :param str kernel: the name of the coded kernel.

    :returns: options selecting an explicit region identifier.
    :rtype: dict

    '''
    return {'region_name': (REGION_MODULE, region_name(kernel))}


def _refuse_inside_openmp(node):
    '''
    Refuses to place a caliper inside an OpenMP region.

    PSyData's calipers are not thread-safe and timer_mod's accumulators are
    module variables toggled by name, so a caliper entered by every thread
    would corrupt every row rather than only its own.

    :param node: the node a caliper is about to wrap.
    :type node: :py:class:`psyclone.psyir.nodes.Node`

    :raises TransformationError: if the node is inside an OpenMP region, which
        means this was called after openmp_parallelise_loops.

    '''
    if (node.ancestor(OMPParallelDirective)
            or node.ancestor(OMPParallelDoDirective)
            or node.ancestor(OMPDoDirective)):
        raise TransformationError(
            "Must apply the timed_region calipers BEFORE "
            "openmp_parallelise_loops function in optimisation script.")


def _module_name(node):
    '''
    The name of the PSy-layer module a node belongs to.

    :param node: any node below the module's Container.
    :type node: :py:class:`psyclone.psyir.nodes.Node`

    :returns: the lower-case module name, or None if there is no Container.
    :rtype: str or None

    '''
    container = node.ancestor(Container)
    return container.name.lower() if container else None


def _target_loop(schedule, kernel):
    '''
    Finds the single loop in a schedule that runs a named coded kernel.

    The outer loop of a colouring is returned and its inner half is not, so
    that a coloured build and an uncoloured build bracket the same work under
    the same name.

    :param schedule: the invoke schedule to search.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`
    :param str kernel: the name of the coded kernel, e.g. 'tri_solve_code'.

    :returns: the loop to time.
    :rtype: :py:class:`psyclone.psyir.nodes.Loop`

    :raises RuntimeError: if the schedule does not hold exactly one such loop,
        rather than timing whichever loop happens to be first.

    '''
    wanted = kernel.lower()
    loops = [loop for loop in schedule.loops()
             if loop.loop_type not in INNER_COLOUR_LOOPS
             and [called.name.lower() for called in loop.coded_kernels()]
             == [wanted]]
    if len(loops) != 1:
        raise RuntimeError(
            f"expected one '{kernel}' loop in '{schedule.name}', "
            f"found {len(loops)}")
    return loops[0]


def time_captured_loops(psyir):
    '''
    Places a caliper round each captured call site this PSy layer holds.

    Called from the global script of a Fortran timed transformation, which
    PSyclone runs once per algorithm. Most algorithms hold none of the seven
    and nothing happens; where one is held, the loop is bracketed under the
    same region name the Kokkos build gives its captured region.

    A row of CAPTURED_REGIONS naming this module must be found. A table entry
    that no longer matches the model is otherwise invisible -- the region is
    simply absent from timer.txt, which looks the same as a region that was
    never entered -- so it fails the build instead.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :raises RuntimeError: if a row naming this module has no such invoke, or
        that invoke does not hold exactly one loop running the named kernel.
    :raises TransformationError: if a caliper would land inside an OpenMP
        region, which means this was called after openmp_parallelise_loops.

    '''
    _region_names()
    profile_trans = ProfileTrans()

    schedules = {schedule.name.lower(): schedule
                 for schedule in psyir.walk(InvokeSchedule)}
    here = {_module_name(schedule) for schedule in schedules.values()}

    for module, invoke, kernel in CAPTURED_REGIONS:
        if module.lower() not in here:
            continue
        schedule = schedules.get(invoke.lower())
        if schedule is None:
            raise RuntimeError(
                f"CAPTURED_REGIONS names invoke '{invoke}' in '{module}', "
                f"which holds {sorted(schedules)}")
        loop = _target_loop(schedule, kernel)
        _refuse_inside_openmp(loop)
        profile_trans.apply(loop, options=_options(kernel))


def time_kokkos_call(psyir, invoke, kernel):
    '''
    Places a caliper round the call a captured Kokkos region left behind.

    LFRicKokkosTrans replaces a loop with a call to the generated region, so
    there is no loop left for time_captured_loops to find. The call is located
    by the kernel it came from rather than by being the only call in the
    schedule -- capturing one loop in gungho_model leaves fourteen calls where
    there were none, so "the only call" is not a way to find anything.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param str invoke: the name of the schedule holding the captured region.
    :param str kernel: the name of the coded kernel that was captured, e.g.
        'tri_solve_code'. The same string the Kokkos script passes to
        kokkos_region.capture, so the caliper cannot end up named for a
        different kernel than the one that was captured.

    :raises RuntimeError: if this call site is not in CAPTURED_REGIONS, if the
        PSy layer does not hold exactly one such schedule, or if that schedule
        does not hold exactly one call to the captured region.

    '''
    _region_names()
    base = region_name(kernel)

    schedules = [schedule for schedule in psyir.walk(InvokeSchedule)
                 if schedule.name.lower() == invoke.lower()]
    if len(schedules) != 1:
        raise RuntimeError(
            f"expected one '{invoke}' schedule, found {len(schedules)}")

    # The Fortran builds place their calipers from the table alone, so a call
    # site missing from it would be timed here and nowhere else -- a Kokkos
    # row with no baseline beside it, which is worse than no row at all.
    site = (_module_name(schedules[0]), invoke.lower(), kernel.lower())
    known = {(module.lower(), name.lower(), called.lower())
             for module, name, called in CAPTURED_REGIONS}
    if site not in known:
        raise RuntimeError(
            f"{site} is not in CAPTURED_REGIONS, so the Fortran builds would "
            "place no caliper on it and this region would have no baseline. "
            "Add it to timed_region.CAPTURED_REGIONS.")

    # The generated routine is named for the kernel and may carry a precision
    # infix between the two, tri_solve_code becoming tri_solve_r_single_kokkos.
    # Matching both ends rather than the whole name accommodates that without
    # predicting which infix appears.
    calls = [call for call in schedules[0].walk(Call)
             if call.routine
             and call.routine.name.lower().startswith(base)
             and call.routine.name.lower().endswith('_kokkos')]
    if len(calls) != 1:
        found = sorted(call.routine.name for call in schedules[0].walk(Call)
                       if call.routine)
        raise RuntimeError(
            f"expected one call to a '{base}' Kokkos region in '{invoke}', "
            f"found {len(calls)}. The schedule holds: {found}")

    _refuse_inside_openmp(calls[0])
    ProfileTrans().apply(calls[0], options=_options(kernel))


def report(psyir):
    '''
    Prints the regions this module placed, so a build log records them.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    for node in psyir.walk(ProfileNode):
        print(f"timed_region: placed {node.module_name}:{node.region_name} "
              f"in {_module_name(node)}")
