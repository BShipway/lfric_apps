##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The moist_dyn_gas kernel of the traditional-moisture invoke is replaced by a
single call with C linkage, and the C++ that implements it is written into
WORKING_DIR beside the generated PSy layer. The build compiles and links every
generated .cpp it finds there.

The global transformations are repeated here because psyclone_psykal.mk uses
either the local script or the global one for a given algorithm, never both.
Omitting them would leave this algorithm's built-ins generated differently
from the rest of the application, which would confound any comparison against
a build without the Kokkos region.

'''
import os
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.transformations import LFRicKokkosTrans
from psyclone.psyGen import InvokeSchedule

from psyclone_tools import redundant_computation_setval

# The invoke for moisture_formulation_traditional, and the one coded kernel
# within it that the prototype captures.
TARGET_INVOKE = 'invoke_2'
TARGET_KERNEL = 'moist_dyn_gas_code'


def sidecar_path():
    '''
    Works out where the generated translation unit belongs.

    This script sits under the 'psykal' directory at the same relative path as
    the algorithm it transforms, which is how psyclone_psykal.mk pairs the two.
    That relative path is therefore also the algorithm's position in
    WORKING_DIR, so it is taken from here rather than restated.

    :returns: the path to write the generated C++ to.
    :rtype: :py:class:`pathlib.Path`

    :raises RuntimeError: if WORKING_DIR is unset or this script has been
        moved out of the 'psykal' directory.

    '''
    working_dir = os.environ.get('WORKING_DIR')
    if not working_dir:
        raise RuntimeError(
            "WORKING_DIR must name the PSyclone output directory: the "
            "generated Kokkos source has nowhere else to go.")

    script = Path(__file__).resolve()
    for parent in script.parents:
        if parent.name == 'psykal':
            relative = script.relative_to(parent)
            break
    else:
        raise RuntimeError(
            f"{script} is not below a 'psykal' directory, so the position of "
            "the algorithm it transforms cannot be determined.")

    return Path(working_dir) / relative.with_name(f'{script.stem}_kokkos.cpp')


def target_loop(psyir):
    '''
    Finds the single loop the prototype supports.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :returns: the loop to capture.
    :rtype: :py:class:`psyclone.domain.lfric.LFRicLoop`

    :raises RuntimeError: if the algorithm no longer holds exactly one such
        loop, rather than transforming whichever loop happens to be first.

    '''
    schedules = [schedule for schedule in psyir.walk(InvokeSchedule)
                 if schedule.name == TARGET_INVOKE]
    if len(schedules) != 1:
        raise RuntimeError(
            f"expected one '{TARGET_INVOKE}' schedule, found "
            f"{len(schedules)}")

    loops = [loop for loop in schedules[0].walk(LFRicLoop)
             if [kernel.name.lower() for kernel in loop.kernels()]
             == [TARGET_KERNEL]]
    if len(loops) != 1:
        raise RuntimeError(
            f"expected one '{TARGET_KERNEL}' loop in '{TARGET_INVOKE}', "
            f"found {len(loops)}")
    return loops[0]


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)

    source = LFRicKokkosTrans().apply(target_loop(psyir))

    sidecar = sidecar_path()
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(source)
    print(f"Kokkos: generated {sidecar}")
