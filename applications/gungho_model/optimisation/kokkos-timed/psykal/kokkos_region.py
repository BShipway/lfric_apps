##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Shared body of the local PSyclone scripts that capture a Kokkos region, with a
timing caliper round each captured region.

The 'kokkos-timed' copy of 'kokkos'. It differs from it in this file and in
global.py and nowhere else -- the seven per-algorithm scripts are byte-identical
to their counterparts, so

    diff -r optimisation/kokkos/psykal optimisation/kokkos-timed/psykal

is the whole of the difference between the two transformations, and a script
added to one tree can be copied to the other unchanged.

Each of those scripts names one invoke and one coded kernel within it. The
rest is the same wherever a region is captured: repeat the global
transformations, replace the loop, wrap the call that replaced it in a PSyData
caliper, and write the generated C++ into WORKING_DIR beside the PSy layer,
where compile.mk finds it.

The global transformations are repeated because psyclone_psykal.mk uses either
the local script or the global one for a given algorithm, never both. Omitting
them would leave a captured algorithm's built-ins generated differently from
the rest of the application, which would confound any comparison against a
build without the Kokkos regions.

The caliper is applied here rather than through psyclone_tools.profile_loops
because there is no loop left to profile: the region is a Call. See
timed_region for that and for how its region name is made to join with the
Fortran builds'.

This module sits at the root of the 'psykal' directory rather than beside the
scripts that use it, because their position under that root is not free: it is
how psyclone_psykal.mk pairs each script with its algorithm.

'''
import os
import sys
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule

from psyclone_tools import redundant_computation_setval

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
    Applies the global transformations, captures the target loop, and brackets
    the region it leaves behind with a timing caliper.

    The caliper goes on after the capture, not before: the loop it would have
    wrapped no longer exists once LFRicKokkosTrans has run, and a caliper
    placed first would be bracketing a node the transformation is about to
    replace.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    :param str script: the __file__ of the calling script.
    :param str invoke: the name of the schedule holding the loop.
    :param str kernel: the lower-case name of the coded kernel in it.

    '''
    redundant_computation_setval(psyir)

    source = LFRicKokkosTrans().apply(target_loop(psyir, invoke, kernel))

    time_kokkos_call(psyir, invoke, kernel)
    report(psyir)

    sidecar = sidecar_path(script)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(source)
    print(f"Kokkos: generated {sidecar}")
