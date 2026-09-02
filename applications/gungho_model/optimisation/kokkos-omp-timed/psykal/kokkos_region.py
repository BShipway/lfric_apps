##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Shared body of the local PSyclone scripts that capture a Kokkos region, in a
build that also colours and OpenMP-threads everything the capture leaves, with
a timing caliper round each captured region.

The 'kokkos-omp-timed' copy of 'kokkos-omp'. It differs from it in this file
and in the profiling-library marker and nowhere else -- global.py and the seven
per-algorithm scripts are byte-identical to their counterparts, so

    diff -r optimisation/kokkos-omp/psykal optimisation/kokkos-omp-timed/psykal

reports this file alone. That is a smaller difference than the one between
'kokkos' and 'kokkos-timed', which is two files, because the calipers here go
where the capture is and global.py has none either way.

WHY THIS TREE EXISTS SEPARATELY FROM 'kokkos-omp'

'kokkos-omp' answers the coexistence question and must stay untimed to answer
it: its checksums are the evidence, and they are compared against 'meto-ex1a'
built untimed. This tree is the measurement, and it stands to 'kokkos-omp'
exactly as 'production-timed' stands to 'meto-ex1a' -- the same transformation
with calipers added. The pair exists so that neither role has to be argued
from the other.

THE ORDER IS FOUR CONSTRAINTS AT ONCE, AND ONLY ONE ORDER SATISFIES THEM

    redundant computation -> capture -> colour -> caliper -> OpenMP

Each arrow is forced, and each is enforced by something that raises:

1. Capture before colouring. LFRicKokkosTrans._validate_loop refuses a loop
   that is already coloured -- "LFRicKokkosTrans supports only an uncoloured
   cell-column loop".
2. Colouring before the caliper. psyclone_tools.colour_loops walks a Routine's
   immediate children and raises TransformationError -- "Must apply
   colour_loops BEFORE profile_loops function in optimisation script" -- on
   finding a ProfileNode among them, because a caliper already in place would
   hide the loop it was meant to colour. The caliper this file places is on
   the captured Call, which is an immediate child of the schedule, so it would
   trip that check exactly as a caliper on a loop would.
3. The caliper before OpenMP. timed_region._refuse_inside_openmp declines a
   caliper with an OpenMP ancestor, because PSyData is not thread-safe and
   timer_mod's table is module state toggled by name.
4. Colouring before OpenMP, which is psyclone_tools' own ordering and the one
   'meto-ex1a' uses.

Constraints 1 and 2 pin the capture and the colouring either side of each
other, which is what makes this order the only one available rather than a
preference. Every one of the three refusals names the function to move, so a
wrong order fails the build rather than quietly mismeasuring it.

WHAT THE CALIPER BRACKETS, AND WHY IT JOINS 'production-timed'

The Call, not the region function -- the same node 'kokkos-timed' brackets, so
the crossing is inside the figure on both. 'production-timed' brackets the
outer loop of the colouring for the same kernel under the same region name, so
a row here and a row there cover the same work in the same units. See
timed_region for the naming and for why psyclone_tools.profile_loops is not
used.

The seven captured algorithms reach this file through a local script beside
each of them; global.py is what meets every other algorithm, and how
psyclone_psykal.mk pairs each script with its algorithm.

'''
import os
import sys
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops)

# The caliper helper is shared with the two Fortran timed transformations and
# so lives above this tree, where a single CAPTURED_REGIONS table can be the
# one both sides agree on. The first parent holding it is used, so a copy
# placed inside a transformation tree would be found ahead of the shared one.
_HELPER = [p for p in Path(__file__).resolve().parents
           if (p / 'timed_region.py').is_file()]
if not _HELPER:
    raise RuntimeError(f"{__file__} has no timed_region.py above it")
sys.path.insert(0, str(_HELPER[0]))

from timed_region import time_kokkos_call, report    # noqa: E402


def psykal_root(script):
    '''
    Finds the 'psykal' directory a local transformation script sits under.

    :param str script: the __file__ of the calling script.

    :returns: the root of the transformation tree.
    :rtype: :py:class:`pathlib.Path`

    :raises RuntimeError: if the script has been moved out of that tree, in
        which case the algorithm it transforms cannot be identified.

    '''
    for parent in Path(script).resolve().parents:
        if parent.name == 'psykal':
            return parent
    raise RuntimeError(
        f"{script} is not below a 'psykal' directory, so the position of the "
        "algorithm it transforms cannot be determined.")


def sidecar_path(script):
    '''
    Works out where the generated translation unit belongs.

    A local script sits under the 'psykal' directory at the same relative path
    as the algorithm it transforms, which is how psyclone_psykal.mk pairs the
    two. That relative path is therefore also the algorithm's position in
    WORKING_DIR, so it is taken from the script rather than restated.

    :param str script: the __file__ of the calling script.

    :returns: the path to write the generated C++ to.
    :rtype: :py:class:`pathlib.Path`

    :raises RuntimeError: if WORKING_DIR is unset.

    '''
    working_dir = os.environ.get('WORKING_DIR')
    if not working_dir:
        raise RuntimeError(
            "WORKING_DIR must name the PSyclone output directory: the "
            "generated Kokkos source has nowhere else to go.")

    resolved = Path(script).resolve()
    relative = resolved.relative_to(psykal_root(script))
    return Path(working_dir) / relative.with_name(f'{resolved.stem}_kokkos.cpp')


def target_loop(psyir, invoke, kernel):
    '''
    Finds the single loop the prototype captures.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param str invoke: the name of the schedule holding the loop.
    :param str kernel: the lower-case name of the coded kernel in it.

    :returns: the loop to capture.
    :rtype: :py:class:`psyclone.domain.lfric.LFRicLoop`

    :raises RuntimeError: if the algorithm no longer holds exactly one such
        loop, rather than transforming whichever loop happens to be first.

    '''
    schedules = [schedule for schedule in psyir.walk(InvokeSchedule)
                 if schedule.name == invoke]
    if len(schedules) != 1:
        raise RuntimeError(
            f"expected one '{invoke}' schedule, found {len(schedules)}")

    loops = [loop for loop in schedules[0].walk(LFRicLoop)
             if [called.name.lower() for called in loop.kernels()] == [kernel]]
    if len(loops) != 1:
        raise RuntimeError(
            f"expected one '{kernel}' loop in '{invoke}', found {len(loops)}")
    return loops[0]


def capture(psyir, script, invoke, kernel):
    '''
    Captures the target loop, colours and calipers, then threads the rest.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param str script: the __file__ of the calling script.
    :param str invoke: the name of the schedule holding the loop.
    :param str kernel: the lower-case name of the coded kernel in it.

    '''
    redundant_computation_setval(psyir)

    source = LFRicKokkosTrans().apply(target_loop(psyir, invoke, kernel))

    # Four constraints, one order; see the module docstring. Every other loop
    # in the algorithm is transformed exactly as global.py transforms an
    # uncaptured algorithm's, so the only difference between this build and
    # 'production-timed' is the captured regions.
    colour_loops(psyir)
    time_kokkos_call(psyir, invoke, kernel)
    report(psyir)
    openmp_parallelise_loops(psyir)

    sidecar = sidecar_path(script)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(source)
    print(f"Kokkos: generated {sidecar}")
