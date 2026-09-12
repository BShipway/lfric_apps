##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Capture every cell loop LFRicKokkosTrans accepts as a Kokkos region.

The 'kokkos' transformation captures eight loops, each named by a local script
that sits under psykal/ at its algorithm's position. That is the right shape
for eight and the wrong shape for six hundred: a global script reaches every
algorithm, so the 'kokkos-all' trees have one global.py that calls capture()
below and no local scripts at all.

WHAT IS CAPTURED

Every coded-kernel cell loop whose LFRicKokkosTrans.validate passes, in every
invoke of the algorithm PSyclone is processing. A loop validate refuses is left
as Fortran and is not recorded here: the coverage survey in psy-ir-aidev
records why, per loop, and this module does not duplicate that census. What
is recorded is the other three outcomes -- captured, skipped on purpose, and
the one that should not happen: validate passed and apply then failed.

WHERE THE GENERATED C++ GOES, AND WHY NOT BESIDE THE PSY LAYER

A local script writes <algorithm>_kokkos.cpp beside <algorithm>_psy.f90, which
works because each algorithm captures at most one region. A whole-model capture
breaks that twice over. One algorithm may capture several loops, and -- the
harder case -- two algorithms may capture the same kernel, and the region is an
extern "C" symbol named for the kernel implementation: dg_matrix_vector_code
captured from both diagnostic_alg_mod and semi_implicit_solver_alg_mod would
be two definitions of dg_matrix_vector_kokkos, which does not link.

So the translation units go to WORKING_DIR/kokkos_regions/<symbol>.cpp, one
file per region symbol, however many call sites reach it. compile.mk collects
every *.cpp under WORKING_DIR with find, so the directory needs no make rule.
The first call site to a symbol writes the file; every later one reads it back
and requires the generated text to be identical, because two call sites that
generate different code for one symbol would be a real finding -- a region
whose shape depends on its caller -- and must fail the build rather than let
the last writer win. PSyclone runs once per algorithm and make runs those in
parallel, so the write is a rename of a temporary and the comparison tolerates
two processes producing the same text at once.

THE MANIFEST

Each algorithm writes WORKING_DIR/kokkos_regions/<psy module>.tsv holding one
row per call site it captured:

    <psy module>\t<invoke>\t<kernel>\t<region symbol>\t<written|shared>

and one comment row per site it did not capture for a reason worth recording:

    # skipped\t<psy module>\t<invoke>\t<kernel>\t<reason from SKIP>
    # unmodelled\t<psy module>\t<invoke>\t<kernel>\t<the error apply raised>

psy-ir-aidev's bin/capture-manifest concatenates the fragments, and the tests
that used to hand-list captured regions read that instead. The fragments live
in kokkos_regions/ so that the build's PSyclone-regeneration purge can remove
the whole directory at once: a stale fragment would otherwise report a region
that no longer exists.

SKIP

A site that must not be captured, with the reason written beside it. It exists
for one purpose: when a whole-model build's checksums move, the divergence is
bisected by skipping the newly captured sites in groups and regenerating the
PSy layer, which is minutes rather than the hours a model rebuild takes. An
entry with no reason is refused at import, so nothing can be skipped quietly,
and the phase-5 exit criterion required this table to be empty.

Phase 6 put two entries in it, and they are not a bisection left behind. The
two pointwise_convert_xyz2llr sites in lfric_xios_setup_mod convert the nodal
coordinates of a function space to longitude and latitude, once, while the
I/O context is set up, and hand them to XIOS as a domain's lonvalue and
latvalue. XIOS then matches those against the domain's cell bounds, which
lfric_xios_setup_mod computes on the host through the same Fortran xyz2llr,
by hashing the coordinates. A kernel run on a GPU computes atan with the
device's libm, which differs from glibc's in the last bit for some arguments,
and the hashes of one bit apart do not match: the first whole-model runs on
an H100 (psy-ir-aidev, validation/device/README.md, 2026-09-11 and -12)
segfaulted inside xios::CMesh::createMeshEpsilon on the first write, and ran
ten timesteps to a checksum within acceptance with these two sites skipped.
The consumer requires the producer's bits to equal the host's, which no
device transcendental can promise, so the two sites stay in Fortran. They
run once and take no measurable time.

COLOURED BUILDS

capture(psyir, coloured=True) colours a loop before capturing it, whenever
its kernel writes a field two cells share. That is the alternative answer to
a shared write: LFRicKokkosTrans generates a Kokkos::atomic_* for such a
write by default, and takes a plain statement when the loop it is given is
already coloured, because the cells of one colour meet at no dof.

WHICH LOOPS ARE OFFERED THE COLOURING

Every loop whose kernel LFRicKokkosTrans._shared_arguments answers for, which
is the transformation's own definition of a write two cells may make to one
dof and is asked of it rather than restated here. Two things put an argument
in that answer: an access of gh_inc or gh_readinc, under which the cells add
to a dof, and a written field on a space this library does not call
discontinuous, under which the cells each store one. The second half arrived
with task D4, which stopped refusing a continuous gh_write and gave it
Kokkos::atomic_store; until task D4.1 this module asked only about the first,
so those stores reached the coloured build as atomics -- 18 regions of them,
at 27 call sites -- and the device comparison the coloured build exists for
would have compared atomics with atomics at every one.

WHAT THE COLOURED BUILD CAPTURES AND THE DEFAULT ONE DOES NOT

An atomic update answers a read-modify-write of one element and nothing else.
A kernel that updates a shared field with a whole-array expression, or with a
statement that is not one of the shapes ATOMIC_UPDATES names, is refused by
LFRicKokkosTrans.validate on the atomic arm with 'Colour the loop instead' --
twenty-four call sites of the model, in nine kernels. Colouring is the answer
the refusal names: it puts the cells that meet at a dof in different launches,
so no shape is required of the update at all and validate asks nothing about
it. So this build asks the coloured question first for a loop whose kernel
has a shared write: it colours a copy of the invoke schedule, validates the
copy's coloured loop, and colours the real loop only when that answered yes.
A loop the coloured arm refuses, and one LFRicColourTrans will not colour,
falls back to the uncoloured validate and the atomic arm, so nothing the
default build captures is lost.

Ordering the two arms this way is why the two builds no longer capture the
same sites: 'kokkos-all-coloured' captures those twenty-four and 'kokkos-all'
does not. The continuous stores are a second and different difference: both
builds capture them, and only the coloured build captures them without an
atomic.

WHY 'kokkos-all' IS NOT CHANGED TO MATCH

'kokkos-all' is the tree whose checksums are bit-identical to the Fortran
reference, and it has to stay so: it is the control every other capture
profile is measured against, and a control that moved when the capture grew
would prove nothing about the growth. Colouring changes the order the
contributions to a shared dof are summed, which is a floating-point
difference and therefore a checksum difference; that is expected of a
coloured build and would be a regression in the default one. So the default
build keeps the atomic arm, keeps the sites it captured before, and refuses
the twenty-four, and the coloured build is where the wider capture lives.

Colouring a store is the one case that moves no arithmetic. A gh_write to a
continuous space is legal LFRic because its author promises that every cell
reaching a shared dof stores the same value there, so the order the stores
happen in cannot change what the dof ends up holding. Neither this module nor
the transformation relies on that promise -- the default build still emits
Kokkos::atomic_store for every one of them -- but it is why the coloured
build's checksums are expected to move only where an accumulation was
recoloured.

THE INVARIANT: NO COLOURED FORTRAN IS LEFT BEHIND

A loop this build does not capture must be left as the Fortran loop it was.
Colouring only after a dry run on a copy is what keeps that true for a loop
validate refuses: the real loop is never coloured unless the identical
colouring has already been shown to validate. The one case that survives is
a loop validate accepts and apply then fails on, which is recorded as
'# unmodelled' and is the outcome this module already treats as a defect to
be reported; _check_no_colouring_left asserts at the end of the capture that
every coloured loop still standing is one of those and raises if it is not.

TIMED BUILDS

capture(psyir, timed=True) places a caliper on each captured call whose site is
in timed_region.CAPTURED_REGIONS, under the same region name the 'minimum-timed'
and 'kokkos-timed' builds use, so the three join. It does not caliper every
captured region: timer_mod holds 300 timers and LFRic's own claim some, and
timed_region.py's docstring records that a caliper on every kernel runs the
table out and stops the model. What a whole-model timed build measures is the
timestep, from LFRic's own timers, with the nine known regions beside it for
comparison against the other timed builds.
'''
import os
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.lfric_builtins import LFRicBuiltIn
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.backend.kokkos_staging import HEADER_NAME, header_text
from psyclone.psyir.nodes import Call, Container
from psyclone.psyir.transformations import TransformationError
from psyclone.transformations import LFRicColourTrans

#: The loop type a colouring leaves the cells of one colour in. Read from the
#: transformation rather than restated, so that the two cannot drift.
# pylint: disable-next=protected-access
COLOURED_LOOP_TYPE = LFRicKokkosTrans._COLOURED_LOOP_TYPE

#: Where the generated translation units and the manifest fragments go, under
#: WORKING_DIR. Named once here and in psy-ir-aidev's bin/capture-manifest and
#: bin/build-gungho-model, which read and purge it.
REGIONS_DIR = 'kokkos_regions'

#: Sites that are not captured, as (psy module, invoke, kernel) -> reason. All
#: three names lower-case. Empty is the intended state; see the docstring.
SKIP = {
    ('lfric_xios_setup_mod_psy', 'invoke_2_pointwise_convert_xyz2llr_kernel_type',
     'pointwise_convert_xyz2llr_code'):
        'XIOS hashes these longitudes and latitudes against bounds the host '
        'computes with glibc atan; a device atan one bit apart does not '
        'match (phase 6, H4)',
    ('lfric_xios_setup_mod_psy', 'invoke_4_pointwise_convert_xyz2llr_kernel_type',
     'pointwise_convert_xyz2llr_code'):
        'XIOS hashes these longitudes and latitudes against bounds the host '
        'computes with glibc atan; a device atan one bit apart does not '
        'match (phase 6, H4)',
}


def _check_skip():
    '''
    Refuses a SKIP entry that carries no reason.

    :raises ValueError: if any entry's reason is empty or not a string, or if
        any key is not a three-tuple of lower-case names.

    '''
    for site, reason in SKIP.items():
        if (not isinstance(site, tuple) or len(site) != 3
                or any(not isinstance(part, str) or part != part.lower()
                       for part in site)):
            raise ValueError(
                f"SKIP key {site!r} must be a (module, invoke, kernel) "
                "tuple of lower-case names.")
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError(
                f"SKIP entry {site} has no reason. A site is skipped for a "
                "reason that is written down, or it is not skipped.")


_check_skip()


def working_dir():
    '''
    The PSyclone output directory the build passes as WORKING_DIR.

    :returns: the working directory.
    :rtype: :py:class:`pathlib.Path`

    :raises RuntimeError: if WORKING_DIR is unset.

    '''
    value = os.environ.get('WORKING_DIR')
    if not value:
        raise RuntimeError(
            "WORKING_DIR must name the PSyclone output directory: the "
            "generated Kokkos source has nowhere else to go.")
    return Path(value)


def _module_name(schedule, psyir):
    '''
    The lower-case name of the PSy-layer module holding a schedule.

    :param schedule: the invoke schedule.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`
    :param psyir: the file container, whose name is the fallback.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :returns: the module name.
    :rtype: str

    '''
    container = schedule.ancestor(Container)
    return (container.name if container else psyir.name).lower()


def _write_region(regions, symbol, source, site):
    '''
    Writes one region's translation unit, or checks it against the existing.

    :param regions: the kokkos_regions directory.
    :type regions: :py:class:`pathlib.Path`
    :param str symbol: the region's extern "C" name.
    :param str source: the generated C++.
    :param tuple[str, str, str] site: the (module, invoke, kernel) capturing it,
        for the message if the text disagrees.

    :returns: 'written' if this call site created the file, 'shared' if an
        earlier one had and the text agrees.
    :rtype: str

    :raises RuntimeError: if the file exists with different text, which means
        two call sites generate different code for one symbol.

    '''
    path = regions / f'{symbol}.cpp'
    if path.exists():
        existing = path.read_text(encoding='utf-8')
        if existing != source:
            raise RuntimeError(
                f"{'/'.join(site)} generated {symbol} differently from an "
                f"earlier call site; see {path}. One extern \"C\" symbol "
                "cannot have two bodies, and a region whose code depends on "
                "its caller is a finding, not something to overwrite.")
        return 'shared'
    temporary = regions / f'.{symbol}.{os.getpid()}.tmp'
    temporary.write_text(source, encoding='utf-8')
    os.replace(temporary, path)
    return 'written'


def _write_staging_header(regions):
    '''
    Writes the staging header the generated regions include.

    Every region PSyclone generates now obtains its Views through
    lfric_kokkos::stage(), declared in lfric_kokkos_staging.hpp, so the header
    has to sit beside them: compile.mk collects the .cpp files under
    WORKING_DIR and compiles each in its own directory, which is where the
    quoted include looks first.

    The text comes from PSyclone rather than from a copy kept here, for the
    same reason COLOURED_LOOP_TYPE is read from the transformation: a header
    that drifted from the regions generated against it would be a mismatch no
    build reports. Written the same way a region is -- a temporary and a
    rename -- because make runs PSyclone once per algorithm, in parallel, and
    every one of them writes this.

    :param regions: the kokkos_regions directory.
    :type regions: :py:class:`pathlib.Path`

    '''
    path = regions / HEADER_NAME
    source = header_text()
    if path.exists() and path.read_text(encoding='utf-8') == source:
        return
    temporary = regions / f'.{HEADER_NAME}.{os.getpid()}.tmp'
    temporary.write_text(source, encoding='utf-8')
    os.replace(temporary, path)


def _new_region_call(schedule, before):
    '''
    Finds the launch call apply() left where the loop was.

    :param schedule: the invoke schedule that held the loop.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`
    :param set[int] before: ids of the Calls in the schedule before apply().

    :returns: the region's extern "C" symbol, lower-case.
    :rtype: str

    :raises RuntimeError: if apply() left anything other than exactly one new
        call to a *_kokkos routine.

    '''
    calls = [call for call in schedule.walk(Call)
             if id(call) not in before and call.routine
             and call.routine.name.lower().endswith('_kokkos')]
    if len(calls) != 1:
        raise RuntimeError(
            f"expected apply() to leave one new *_kokkos call in "
            f"'{schedule.name}', found {len(calls)}")
    return calls[0].routine.name.lower()


def _write_manifest(regions, module, rows):
    '''
    Writes this algorithm's manifest fragment.

    :param regions: the kokkos_regions directory.
    :type regions: :py:class:`pathlib.Path`
    :param str module: the PSy module, which names the fragment.
    :param list[tuple] rows: the rows, each a tuple of strings; a first
        element starting with '#' is a comment row.

    '''
    path = regions / f'{module}.tsv'
    text = ''.join('\t'.join(row) + '\n' for row in rows)
    temporary = regions / f'.{module}.{os.getpid()}.tsv.tmp'
    temporary.write_text(text, encoding='utf-8')
    os.replace(temporary, path)


def _shares_a_write(kernel):
    '''
    Whether two cells of one launch would write the same dof.

    The question is asked of the transformation rather than restated here,
    the way COLOURED_LOOP_TYPE is read from it and the way psy-ir-aidev's
    survey runs its rules: the coloured arm has to be offered to exactly the
    kernels LFRicKokkosTrans calls shared, or the build colours a kernel that
    needs no colouring and -- the failure this call site was written for --
    leaves an atomic on one that does.

    Two different things make a write shared and they are read from different
    places. An access of gh_inc or gh_readinc says that cells add to one dof;
    the function space of a written field says that cells each store one.
    This asked only the first until task D4.1, so the coloured build gave
    every continuous gh_write a Kokkos::atomic_store -- which is the answer
    the coloured build exists to be compared against.

    :param kernel: the coded kernel the loop calls.
    :type kernel: :py:class:`psyclone.domain.lfric.LFRicKern`

    :returns: whether any argument is written by more than one cell of a
        launch, whether by accumulating into it or by storing to it.
    :rtype: bool

    '''
    # LFRicKokkosTrans._shared_arguments is the transformation's own answer
    # to this question and there is no public statement of it. The survey
    # reads _validate_shared_updates the same way and for the same reason.
    # pylint: disable-next=protected-access
    return bool(LFRicKokkosTrans._shared_arguments(kernel))


def _colour(loop):
    '''
    Colours a loop in place and returns the loop over the cells of one colour.

    The colouring replaces the loop with a loop over colours holding a loop
    over that colour's cells, at the position the loop had. Finding the inner
    loop by that position rather than by the identity of the kernel is what
    lets the same helper colour a copy of a schedule, where no node is the
    same object as the one it stands for.

    :param loop: the loop to colour.
    :type loop: :py:class:`psyclone.domain.lfric.LFRicLoop`

    :returns: the cells-in-colour loop.
    :rtype: :py:class:`psyclone.domain.lfric.LFRicLoop`

    :raises TransformationError: if LFRicColourTrans refuses the loop, which
        the caller answers by leaving the loop uncoloured.
    :raises RuntimeError: if the colouring did not leave exactly one loop over
        the cells of a colour where the loop was.

    '''
    parent, position = loop.parent, loop.position
    LFRicColourTrans().apply(loop)
    inner = [each for each in parent.children[position].walk(LFRicLoop)
             if each.loop_type == COLOURED_LOOP_TYPE]
    if len(inner) != 1:
        raise RuntimeError(
            f"colouring left {len(inner)} loops over the cells of a colour "
            f"where one loop was, at position {position}")
    return inner[0]


def _refuses_coloured(loop, schedule):
    '''
    Whether the transformation would refuse this loop once it was coloured.

    Asked before the loop itself is coloured, and answered by colouring a copy
    of the PSy layer and validating the copy's coloured loop. That is
    what keeps the module docstring's invariant: the real loop is coloured only
    when the identical colouring has already been shown to validate, so a loop
    the transformation refuses is left as the Fortran loop it was rather than
    as a coloured one no region replaces.

    The copy is of the whole PSy layer rather than of the invoke schedule
    alone, and the loop is found again by its position, for the reason
    LFRicKokkosTrans._rooted_copy gives about a kernel schedule: a schedule
    copied on its own is detached from the FileContainer above it, and
    LFRicKern.get_callees looks for a local implementation of the kernel in
    that ancestor Container before it reads the kernel file. Validating a loop
    inside a detached copy raises out of PSyclone rather than answering.

    Colouring the copy calls colourmap_init on the invoke, which is not itself
    copied. That reads the invoke's real schedule for kernels that are
    coloured, and the loop this is asked about is not one, so a dry run adds
    no colourmap symbols to a PSy layer that would then not use them.

    :param loop: the uncoloured loop.
    :type loop: :py:class:`psyclone.domain.lfric.LFRicLoop`
    :param schedule: the invoke schedule holding it.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`

    :returns: the refusal, or None if the coloured loop would be captured.
    :rtype: Optional[str]

    '''
    root = schedule.root
    index = root.walk(LFRicLoop).index(loop)
    try:
        inner = _colour(root.copy().walk(LFRicLoop)[index])
    except TransformationError as err:
        return f'LFRicColourTrans: {err}'
    try:
        LFRicKokkosTrans().validate(inner)
    except TransformationError as err:
        return str(err)
    return None


def _check_no_colouring_left(psyir, colours, captured, unmodelled):
    '''
    Assert the invariant: no loop this build left as Fortran was coloured.

    A loop the transformation refuses is never coloured, because the colouring
    is tried on a copy first, and a captured one has had its cells-in-colour
    loop replaced by the launch call. So the only coloured loop that can still
    be standing is one validate accepted and apply then failed on, which is
    recorded as '# unmodelled' and reported as the defect it is. Anything else
    is a colouring on a loop this build left as Fortran, which would make the
    coloured tree's PSy layer differ from the default tree's at a site neither
    of them captured.

    Asked twice, of two different things. The bookkeeping asks which sites were
    coloured and did not end as a capture or an unmodelled site; the walk asks
    the generated tree whether a cells-in-colour loop is standing anywhere it
    should not be. The first is immune to a coloured loop having been lowered
    out of LFRicLoop by a sibling's capture and the second is immune to the
    bookkeeping being wrong, and each is cheap.

    :param psyir: the PSyIR of the PSy layer, after the capture.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param set[tuple[str, str, str]] colours: the sites that were coloured.
    :param set[tuple[str, str, str]] captured: the sites that were captured.
    :param set[tuple[str, str, str]] unmodelled: the sites apply() failed on.

    :raises RuntimeError: if a colouring is standing at any other site.

    '''
    left = {'/'.join(site) for site in colours - captured - unmodelled}
    for schedule in psyir.walk(InvokeSchedule):
        module = _module_name(schedule, psyir)
        for loop in schedule.walk(LFRicLoop):
            if loop.loop_type != COLOURED_LOOP_TYPE:
                continue
            for kernel in loop.kernels():
                site = (module, schedule.name.lower(), kernel.name.lower())
                if site not in unmodelled:
                    left.add('/'.join(site))
    if left:
        raise RuntimeError(
            "the coloured capture left a colouring behind at a site it did "
            f"not capture: {', '.join(sorted(left))}. A loop this build does "
            "not capture must be the Fortran loop it was.")


def capture(psyir, timed=False, coloured=False):
    '''
    Captures every loop the transformation accepts, in every invoke.

    Called from the global script after the global transformations, which
    is the order the 'kokkos' profile's local scripts also use: the redundant
    computation transformation moves loop bounds, and a loop must be captured
    with the bounds it will run with.

    :param psyir: the PSyIR of the PSy layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param bool timed: place a caliper on each captured call whose site is in
        timed_region.CAPTURED_REGIONS.
    :param bool coloured: colour a loop whose kernel writes a field two cells
        share before capturing it, so that the region takes the coloured
        answer to a shared write rather than the atomic one. The colouring is
        tried first and the atomic arm is the fallback, which is what makes a
        coloured build capture the sites an atomic has no shape for; see the
        module docstring.

    :raises RuntimeError: if two call sites generate one symbol differently,
        if apply() leaves other than one new call, if a coloured build leaves
        a colouring on a loop it did not capture, or if a timed build finds
        a CAPTURED_REGIONS site in this module that was not captured -- the
        timed builds must bracket the same work under the same names.

    '''
    if timed:
        # Imported here rather than at the top so that a plain build never
        # depends on the timing helper being importable.
        import timed_region  # pylint: disable=import-outside-toplevel

    regions = working_dir() / REGIONS_DIR
    regions.mkdir(parents=True, exist_ok=True)
    # Before any region is written, so that a build interrupted between the
    # two leaves a region with no header rather than a header with no region:
    # the first does not compile and is seen, the second compiles and is not.
    _write_staging_header(regions)

    rows_by_module = {}
    captured_sites = set()
    coloured_sites = set()
    unmodelled_sites = set()
    for schedule in psyir.walk(InvokeSchedule):
        module = _module_name(schedule, psyir)
        rows = rows_by_module.setdefault(module, [])
        # A snapshot, because apply() replaces loops while this iterates.
        for loop in list(schedule.walk(LFRicLoop)):
            kernels = loop.kernels()
            if not kernels or any(isinstance(kernel, LFRicBuiltIn)
                                  for kernel in kernels):
                continue
            site = (module, schedule.name.lower(), kernels[0].name.lower())
            if site in SKIP:
                rows.append(('# skipped',) + site + (SKIP[site],))
                print(f"Kokkos: skipped {'/'.join(site)}: {SKIP[site]}")
                continue
            transformation = LFRicKokkosTrans()
            # The coloured arm is asked first for a loop with a shared write,
            # because it is the arm that takes the updates no atomic answers.
            # It is asked of a copy, and the loop itself is coloured only
            # where the copy validated: a colouring is never left on a loop
            # this build then leaves as Fortran.
            take_colour = coloured and _shares_a_write(kernels[0])
            if take_colour:
                refusal = _refuses_coloured(loop, schedule)
                if refusal is None:
                    loop = _colour(loop)
                    coloured_sites.add(site)
                else:
                    # Either arm may take it; the atomic one is asked next,
                    # so nothing the default build captures is lost here.
                    take_colour = False
                    print(f"Kokkos: uncoloured {'/'.join(site)}: {refusal}")
            if not take_colour:
                try:
                    transformation.validate(loop)
                except TransformationError:
                    continue
            before = {id(call) for call in schedule.walk(Call)}
            try:
                source = transformation.apply(loop)
            except Exception as err:            # pylint: disable=broad-except
                reason = f'{err.__class__.__name__}: {err}'
                rows.append(('# unmodelled',) + site + (reason,))
                unmodelled_sites.add(site)
                print(f"Kokkos: UNMODELLED {'/'.join(site)}: {reason}")
                continue
            symbol = _new_region_call(schedule, before)
            status = _write_region(regions, symbol, source, site)
            rows.append(site + (symbol, status))
            captured_sites.add(site)
            print(f"Kokkos: captured {symbol} at {'/'.join(site)} ({status})")

    if coloured:
        _check_no_colouring_left(psyir, coloured_sites, captured_sites,
                                 unmodelled_sites)

    if timed:
        for module, invoke, kernel in timed_region.CAPTURED_REGIONS:
            site = (module.lower(), invoke.lower(), kernel.lower())
            if site[0] not in rows_by_module:
                continue
            if site not in captured_sites:
                raise RuntimeError(
                    f"timed_region.CAPTURED_REGIONS names {'/'.join(site)}, "
                    "which this build did not capture, so the timed builds "
                    "would not bracket the same work. Capture it or take it "
                    "out of SKIP.")
            timed_region.time_kokkos_call(psyir, invoke, kernel)

    for module, rows in rows_by_module.items():
        _write_manifest(regions, module, rows)
