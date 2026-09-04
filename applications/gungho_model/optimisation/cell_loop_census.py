##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
One shared caliper on every cell loop the captured eight do not already cover.

timed_region.py answers 'how long does this one region take'. This module
answers the question underneath it: how much of a timestep is spent in cell
loops at all, so that a speed-up on one region can be read against the work
that is available to speed up. The two are designed to be added rather than
compared, and the arithmetic they make is

    cell_loops = region:cell_loops  +  sum of the eight region:<kernel> rows

with region:cell_loops placed here and the eight placed by timed_region.

WHY ONE NAME AND NOT 645

'A caliper on every cell loop' cannot be done one-caliper-one-name. From the
generated PSy layer of a built 'minimum':

    grep -rhoiE '^ *DO +cell *=' --include='*_psy.f90' . | wc -l      -> 645
    grep -rhoiE '^ *CALL +[a-z0-9_]+_code *\\(' --include='*_psy.f90' . \\
        | sed -E 's/^ *[Cc][Aa][Ll][Ll] +//; s/ *\\($//' | wc -l      -> 668
    ... the same pipeline with 'sort -u' before 'wc -l'                -> 298
    find . -name '*_psy.f90' | wc -l                                   -> 129

timer_mod declares 'integer(i_def), parameter :: num_subs = 300' and running
out is fatal rather than a warning:

    call log_event( "Run out of timers, increase num_subs", LOG_LEVEL_ERROR )

So a name per call site (645) does not fit, and neither does a name per kernel
(298, plus the timers LFRic claims before any of these are reached, plus
timed_region's eight). timed_region.py's docstring records the 298 and the
ceiling already; what it never had to work out, because it only ever places
eight calipers, is what to do instead.

One shared name fits in one slot. timer(cname) finds a row by name, so 637 call
sites all naming 'region:cell_loops' accumulate into a single row, num_calls
counts the entries, and the table stays at about fifteen rows. That also
disposes of the second objection timed_region.py raises against a full table --
timer() scans the table linearly comparing 128-byte names, so a 300-row table
costs hundreds of string comparisons on every caliper entry and exit -- because
the table never grows.

The PSyData wrapper permits this and does not have to be persuaded to.
psyclone/lib/profiling/lfric_timer/profile_psy_data_mod.F90 builds
module_name//":"//region_name per instance and calls timer(name); no registry
anywhere assumes a name is claimed once, and
PSyDataTrans.get_unique_region_name returns a user-supplied name verbatim,
saying in its own docstring that it is then up to the user to guarantee
uniqueness. 637 instances with one name is 637 matched toggles of one row.

WHY A NESTED CALIPER IS REFUSED RATHER THAN TOLERATED

Sharing the name is what makes nesting dangerous, and the danger is silence.
timer(cname) is a toggle on start_stop(k): an odd call starts the row, an even
call stops it and accumulates. output_timer refuses a row left open --
'Timer for routine X not closed.' at LOG_LEVEL_ERROR -- so an *odd* number of
stray toggles is fatal and would be found.

A nested pair is even, and is not:

    outer caliper enters at t0    starts
    inner caliper enters at t1    STOPS, accumulating t1-t0
    inner caliper exits  at t2    STARTS
    outer caliper exits  at t3    stops, accumulating t3-t2

The row closes cleanly and holds the outer duration *minus* the inner one. No
check anywhere in LFRic catches it: the build succeeds, the model runs, and the
number is wrong by an amount nothing reports. _refuse_nesting below is
therefore a build-time failure and not a warning, in the way
_refuse_inside_openmp is.

WHY THE CAPTURED EIGHT ARE EXCLUDED RATHER THAN NESTED

The obvious arrangement -- a shared caliper on all 645 with timed_region's
eight nested inside -- is exactly the arrangement the table above forbids for
a shared name. It would also make 'captured / cell_loops' a ratio of a
measured number to a difference, so an error in either term would move the
ratio twice.

This module partitions instead. The shared caliper goes on every outermost cell
loop EXCEPT the eight named in timed_region.CAPTURED_REGIONS, which
time_captured_loops keeps placing exactly as it does today. Both terms of the
sum are then measured directly, no PSyData region is ever inside another, and
timed_region.py is not edited -- so a census build and a plain '-timed' build
place the same eight calipers on the same eight loops, and their rows join.

The exclusion is derived from CAPTURED_REGIONS through
timed_region._target_loop rather than from what is already wrapped, so it does
not depend on the order the two are called in. A candidate that turns out to
be inside a ProfileNode anyway is a fact this module does not understand, and
it is refused rather than skipped.

WHERE THIS FILE SITS

Beside timed_region.py at the root of optimisation/, for the reason given at
the end of that file: the census transformations import both, and the region
names have to agree between them. The scripts reach it by walking their own
parents for the directory holding it, rather than by a fixed number of
'.parent's.
'''
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.nodes import Loop, ProfileNode
from psyclone.psyir.transformations import ProfileTrans

from timed_region import (
    CAPTURED_REGIONS,
    INNER_COLOUR_LOOPS,
    REGION_MODULE,
    _module_name,
    _refuse_inside_openmp,
    _target_loop,
)

#: The region half of the one identifier this module places. Every caliper it
#: places names it, so timer.txt holds one 'region:cell_loops' row however many
#: call sites were wrapped, and num_calls on that row is the number of entries.
CENSUS_REGION = 'cell_loops'


def _captured_loops(schedule):
    '''
    The loops in a schedule that timed_region times, and this module must not.

    Derived from CAPTURED_REGIONS rather than from what is already wrapped, so
    that the census places the same calipers whether it runs before or after
    time_captured_loops.

    :param schedule: the invoke schedule to search.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`

    :returns: the loops timed_region would bracket in this schedule.
    :rtype: list[:py:class:`psyclone.psyir.nodes.Loop`]

    :raises RuntimeError: if a row of CAPTURED_REGIONS names this module and
        invoke but the schedule does not hold exactly one loop running its
        kernel, which would mean the partition has a term it cannot measure.

    '''
    module = _module_name(schedule)
    name = schedule.name.lower()
    return [_target_loop(schedule, kernel)
            for row_module, row_invoke, kernel in CAPTURED_REGIONS
            if row_module.lower() == module and row_invoke.lower() == name]


def _census_loops(schedule):
    '''
    The loops in a schedule this module places a shared caliper on.

    A loop qualifies when it runs at least one coded kernel, is not the inner
    half of a colouring, and is not one of the loops timed_region brackets. The
    inner half of a colouring is excluded for the reason psyclone_tools
    excludes it: timing it would time a fraction of the work and would not
    correspond to anything in an uncoloured build.

    :param schedule: the invoke schedule to search.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`

    :returns: the loops to wrap, outermost first in schedule order.
    :rtype: list[:py:class:`psyclone.psyir.nodes.Loop`]

    :raises RuntimeError: if a row of CAPTURED_REGIONS naming this schedule
        cannot be resolved to exactly one loop.

    '''
    # Identity rather than equality throughout. PSyIR nodes compare
    # structurally, so two loops running the same kernel over the same bounds
    # are equal without being the same loop -- and 'loop not in captured' would
    # then exclude both of them.
    captured = [id(loop) for loop in _captured_loops(schedule)]
    return [loop for loop in schedule.loops()
            if loop.loop_type not in INNER_COLOUR_LOOPS
            and loop.coded_kernels()
            and id(loop) not in captured]


def _refuse_nesting(candidates, schedule):
    '''
    Refuses a census in which one shared caliper would sit inside another.

    Both calipers would name CENSUS_REGION, so the inner pair would stop and
    restart the outer one's row: it would close cleanly and hold the outer
    duration minus the inner one, and nothing in LFRic would report it. See the
    toggle table in this module's docstring.

    :param candidates: the loops about to be wrapped.
    :type candidates: list[:py:class:`psyclone.psyir.nodes.Loop`]
    :param schedule: the invoke schedule they belong to, named in the message.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`

    :raises RuntimeError: if one candidate is an ancestor of another.

    '''
    wanted = {id(loop) for loop in candidates}
    for loop in candidates:
        ancestor = loop.ancestor(Loop)
        while ancestor is not None:
            if id(ancestor) in wanted:
                raise RuntimeError(
                    f"in '{schedule.name}', a '{loop.loop_type}' cell loop is "
                    f"inside a '{ancestor.loop_type}' one and both would be "
                    f"timed as '{REGION_MODULE}:{CENSUS_REGION}'. The LFRic "
                    "timer toggles a row by name, so the inner pair would "
                    "subtract itself from the outer row and the build would "
                    "not fail. Decide which of the two to time and exclude "
                    "the other.")
            ancestor = ancestor.ancestor(Loop)


def _refuse_wrapped(loop, schedule):
    '''
    Refuses a candidate that is already inside a PSyData region.

    The eight captured loops are excluded by name before this is reached, so a
    candidate inside a ProfileNode is one this module does not know about. It
    is refused rather than skipped: skipping it would drop it out of the census
    silently, and the census total is the only evidence that every cell loop
    was reached.

    :param loop: the loop a caliper is about to wrap.
    :type loop: :py:class:`psyclone.psyir.nodes.Loop`
    :param schedule: the invoke schedule it belongs to, named in the message.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`

    :raises RuntimeError: if the loop has a ProfileNode ancestor.

    '''
    node = loop.ancestor(ProfileNode)
    if node is not None:
        raise RuntimeError(
            f"in '{schedule.name}', a cell loop the census would time is "
            f"already inside '{node.module_name}:{node.region_name}'. Nothing "
            "in timed_region.CAPTURED_REGIONS puts it there, so either that "
            "table is out of date or another transformation placed it. The "
            "census refuses to nest rather than to drop the loop.")


def census_cell_loops(psyir):
    '''
    Places the shared census caliper on every uncaptured cell loop.

    Called from the global script of a census transformation, which PSyclone
    runs once per algorithm. Every caliper names the same region, so the whole
    model contributes one row of timer.txt whose accumulated time is the cell
    loop total and whose num_calls is the number of entries.

    Must be called before openmp_parallelise_loops and, where both are used,
    after time_captured_loops -- not because the exclusion depends on the
    order, which it does not, but because _target_loop looks for a loop that
    census calipers would have moved inside a ProfileNode.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :raises RuntimeError: if a row of CAPTURED_REGIONS naming a schedule here
        cannot be resolved, if one candidate is nested inside another, or if a
        candidate is already inside a PSyData region.
    :raises TransformationError: if a caliper would land inside an OpenMP
        region, which means this was called after openmp_parallelise_loops.

    '''
    profile_trans = ProfileTrans()
    options = {'region_name': (REGION_MODULE, CENSUS_REGION)}

    for schedule in psyir.walk(InvokeSchedule):
        candidates = _census_loops(schedule)
        _refuse_nesting(candidates, schedule)
        for loop in candidates:
            _refuse_inside_openmp(loop)
            _refuse_wrapped(loop, schedule)
        # Applied in a second pass so that no caliper is placed at all if any
        # candidate in the schedule is refused. A half-instrumented schedule
        # would still build and would still produce a row.
        for loop in candidates:
            profile_trans.apply(loop, options=options)


def census_report(psyir):
    '''
    Prints what the census placed and what it left to timed_region.

    One line per PSy-layer module. The build log is the only record of the
    count: timer.txt reports one row however many call sites fed it, so the
    number of calipers placed cannot be recovered from a run.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    placed = {}
    skipped = {}
    for schedule in psyir.walk(InvokeSchedule):
        module = _module_name(schedule)
        placed.setdefault(module, 0)
        skipped.setdefault(module, 0)
        skipped[module] += len(_captured_loops(schedule))

    for node in psyir.walk(ProfileNode):
        if (node.module_name == REGION_MODULE
                and node.region_name == CENSUS_REGION):
            module = _module_name(node)
            placed[module] = placed.get(module, 0) + 1

    for module in sorted(placed):
        print(f"cell_loop_census: placed {placed[module]} "
              f"{REGION_MODULE}:{CENSUS_REGION} calipers in {module}, "
              f"leaving {skipped.get(module, 0)} to timed_region")
