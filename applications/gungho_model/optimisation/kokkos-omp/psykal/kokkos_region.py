##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Shared body of the local PSyclone scripts that capture a Kokkos region, in a
build that also colours and OpenMP-threads everything the capture leaves.

The 'kokkos-omp' copy of 'kokkos'. It differs from it in this file and in
global.py and nowhere else -- the seven per-algorithm scripts are
byte-identical to their counterparts, so

    diff -r optimisation/kokkos/psykal optimisation/kokkos-omp/psykal

is the whole of the difference between the two transformations, as it is for
'kokkos-timed'.

'kokkos-omp-timed' is this tree with a caliper round each captured region. It
shares global.py verbatim and differs from this file only in the caliper and
the order the caliper forces; see its copy for the four constraints.

THE ORDER IS THE WHOLE EXPERIMENT

    redundant computation -> capture -> colour -> OpenMP

The capture must come first. LFRicKokkosTrans._validate_loop refuses a loop
that is already coloured -- "LFRicKokkosTrans supports only an uncoloured
cell-column loop" -- so colouring first would make the capture impossible
rather than merely different.

What was not known before this transformation was built is whether the other
direction holds: whether colour_loops and openmp_parallelise_loops walk a
schedule that now holds a Call where a Loop used to be, without complaint.

They do, and the reason is worth stating because it is a property of
psyclone_tools rather than a coincidence. colour_loops iterates a Routine's
immediate children and acts only on nodes that are `isinstance(child, Loop)`;
openmp_parallelise_loops iterates `subroutine.loops()`, which yields Loop nodes
only. A Call is neither, so both walk past the captured region and transform
the loops around it. Neither raises, and neither has to be told to skip it.

colour_loops does have one immediate child it refuses -- a ProfileNode, on
which it raises "Must apply colour_loops BEFORE profile_loops function in
optimisation script". Nothing in this tree places one, so the refusal is never
reached here; it is what fixes the order in 'kokkos-omp-timed', where the
caliper wraps the captured Call and so becomes an immediate child itself.

That the mechanism is 'the Call is not a Loop' also says how it could break: a
future psyclone_tools that walked `subroutine.children` looking for kernels
rather than loops, or that asserted every child was a Loop, would refuse this
build. The refusal would be a real one, not a bug in this file.

WHAT COEXISTENCE DOES AND DOES NOT BUY

The captured region is threaded by Kokkos, on its own thread pool, and the
loops around it by OpenMP. They are not nested and they do not run at the same
time, so this is not a question of oversubscription within one call. It does
mean a build carries two runtimes' worth of threads, and the checksums are the
evidence that the answer is still the model's answer.

'''
import os
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops)


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
    Captures the target loop, then colours and threads what is left of the
    algorithm.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param str script: the __file__ of the calling script.
    :param str invoke: the name of the schedule holding the loop.
    :param str kernel: the lower-case name of the coded kernel in it.

    '''
    redundant_computation_setval(psyir)

    source = LFRicKokkosTrans().apply(target_loop(psyir, invoke, kernel))

    # After the capture, and in this order, for the reasons in the module
    # docstring. Every other loop in the algorithm is transformed exactly as
    # global.py transforms an uncaptured algorithm's, so the only difference
    # between this build and 'production-timed' is the captured regions.
    colour_loops(psyir)
    openmp_parallelise_loops(psyir)

    sidecar = sidecar_path(script)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(source)
    print(f"Kokkos: generated {sidecar}")
