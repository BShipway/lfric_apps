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

The C16_MG configuration does not reach this invoke: its moisture formulation
is selected three levels above the kernel. It is captured all the same, as the
contract the Kokkos backend was built against.

'''
import sys
from pathlib import Path

# The position of this file under 'psykal' is how psyclone_psykal.mk pairs it
# with its algorithm, so it cannot sit beside the helper it shares with the
# other Kokkos scripts. PSyclone puts only this directory on sys.path, so the
# 'psykal' root that holds the helper is added here.
_PSYKAL = [p for p in Path(__file__).resolve().parents if p.name == 'psykal']
if not _PSYKAL:
    raise RuntimeError(f"{__file__} is not below a 'psykal' directory")
sys.path.insert(0, str(_PSYKAL[0]))

from kokkos_region import capture           # noqa: E402  needs the path above

# The invoke for moisture_formulation_traditional, and the one coded kernel
# within it that the prototype captures.
TARGET_INVOKE = 'invoke_2'
TARGET_KERNEL = 'moist_dyn_gas_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
