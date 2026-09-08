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
and the phase-5 exit criterion requires this table to be empty.

COLOURED BUILDS

capture(psyir, coloured=True) colours a loop before capturing it, whenever
its kernel writes a field two cells share. That is the alternative answer to
a shared write: LFRicKokkosTrans generates a Kokkos::atomic_add for such an
update by default, and takes a plain read-modify-write when the loop it is
given is already coloured, because the cells of one colour meet at no dof.
The two builds capture the same sites and differ only in what the region
does with the update, which is what makes their checksums comparable.

Only a loop the transformation has already accepted is coloured, so a
colouring never changes a PSy layer this build then leaves as Fortran, and a
colouring LFRicColourTrans refuses leaves the loop to the atomic arm rather
than dropping the capture.

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

from psyclone.core import AccessType
from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.lfric_builtins import LFRicBuiltIn
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.nodes import Call, Container
from psyclone.psyir.transformations import TransformationError
from psyclone.transformations import LFRicColourTrans

#: The accesses that make two cells of one launch update one dof. Read from
#: the kernel's metadata, which is what LFRicKokkosTrans reads too.
SHARED_ACCESSES = (AccessType.INC, AccessType.READINC)

#: Where the generated translation units and the manifest fragments go, under
#: WORKING_DIR. Named once here and in psy-ir-aidev's bin/capture-manifest and
#: bin/build-gungho-model, which read and purge it.
REGIONS_DIR = 'kokkos_regions'

#: Sites that are not captured, as (psy module, invoke, kernel) -> reason. All
#: three names lower-case. Empty is the intended state; see the docstring.
SKIP = {
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
    Whether two cells of one launch would update the same dof.

    :param kernel: the coded kernel the loop calls.
    :type kernel: :py:class:`psyclone.domain.lfric.LFRicKern`

    :returns: whether any argument is incremented rather than written.
    :rtype: bool

    '''
    return any(argument.access in SHARED_ACCESSES
               for argument in kernel.arguments.args)


def _colour(loop, schedule, kernel):
    '''
    Colours a loop and returns the inner loop, over the cells of one colour.

    A colouring LFRicColourTrans refuses is not an error here. It leaves the
    loop as it was, and LFRicKokkosTrans then generates the atomic update it
    generates for every uncoloured loop, so the site is still captured.

    :param loop: the loop to colour.
    :type loop: :py:class:`psyclone.domain.lfric.LFRicLoop`
    :param schedule: the invoke schedule holding it.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`
    :param kernel: the loop's kernel, which the colouring moves into the
        inner loop and which is how that loop is found again.
    :type kernel: :py:class:`psyclone.domain.lfric.LFRicKern`

    :returns: the cells-in-colour loop, or the original loop if the
        colouring was refused.
    :rtype: :py:class:`psyclone.domain.lfric.LFRicLoop`

    '''
    try:
        LFRicColourTrans().apply(loop)
    except TransformationError as err:
        print(f"Kokkos: not coloured, {kernel.name.lower()}: {err}")
        return loop
    return [inner for inner in schedule.walk(LFRicLoop)
            if inner.loop_type == 'cells_in_colour'
            and any(each is kernel for each in inner.kernels())][0]


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
        answer to a shared write rather than the atomic one.

    :raises RuntimeError: if two call sites generate one symbol differently,
        if apply() leaves other than one new call, or if a timed build finds
        a CAPTURED_REGIONS site in this module that was not captured -- the
        timed builds must bracket the same work under the same names.

    '''
    if timed:
        # Imported here rather than at the top so that a plain build never
        # depends on the timing helper being importable.
        import timed_region  # pylint: disable=import-outside-toplevel

    regions = working_dir() / REGIONS_DIR
    regions.mkdir(parents=True, exist_ok=True)

    rows_by_module = {}
    captured_sites = set()
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
            try:
                transformation.validate(loop)
            except TransformationError:
                continue
            # After validate() and not before it: a colouring is only ever
            # applied to a loop that is about to be captured, so a build
            # never leaves a coloured Fortran loop behind.
            if coloured and _shares_a_write(kernels[0]):
                loop = _colour(loop, schedule, kernels[0])
            before = {id(call) for call in schedule.walk(Call)}
            try:
                source = transformation.apply(loop)
            except Exception as err:            # pylint: disable=broad-except
                reason = f'{err.__class__.__name__}: {err}'
                rows.append(('# unmodelled',) + site + (reason,))
                print(f"Kokkos: UNMODELLED {'/'.join(site)}: {reason}")
                continue
            symbol = _new_region_call(schedule, before)
            status = _write_region(regions, symbol, source, site)
            rows.append(site + (symbol, status))
            captured_sites.add(site)
            print(f"Kokkos: captured {symbol} at {'/'.join(site)} ({status})")

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
