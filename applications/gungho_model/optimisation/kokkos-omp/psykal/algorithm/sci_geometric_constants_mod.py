##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The get_dz_w3 kernel computes the vertical spacing of the W3 cells from the
heights at the W2 points either side of them. It was chosen because C16_MG
measurably enters it, and because the existing backend can already express it:
the loop is over owned cells, its written field is on W3, and its body is
arithmetic on scalars and dofmap-indexed fields with no array sections and no
calls out.

The invoke runs once per mesh, and the field it writes is held in the
geometric-constants inventory and read by vertical FFSL transport throughout
the run. A wrong answer from the generated region therefore shows up in the
model's own checksums rather than only in a counter.

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

# The sole invoke of get_dz_w3, and the one coded kernel within it.
TARGET_INVOKE = 'invoke_12_get_dz_w3_kernel_type'
TARGET_KERNEL = 'get_dz_w3_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
