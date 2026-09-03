##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The fv_difference_z kernel takes the vertical finite-volume difference of a
mass flux: one W3 column is the difference of the W2v fluxes on the faces
either side of it. It writes a whole column as a single array section, which
LFRicKokkosTrans lowers to an explicit loop before generating the region.

It was chosen because it is the busiest loop C16_MG runs that no other
capability blocks. This invoke is entered 18 times per timestep, against
get_dz_w3's 2 times per run, so it is the first captured region to be inside
the timestep loop at all rather than in the model's setup.

The kernel is called from two algorithms. This is the hotter of the two: the
FFSL vertical update, which computes an increment from the flux at every
split step. The other, in flux_precomputations_alg_mod, is left alone, so
that a comparison between builds attributes a difference to one region.

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

# The first of this algorithm's invokes of fv_difference_z, and the one coded
# kernel within it.
TARGET_INVOKE = 'invoke_0_fv_difference_z_kernel_type'
TARGET_KERNEL = 'fv_difference_z_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
