##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The inject_wt_to_sh_w3 kernel moves a field from the Wtheta space of the
ordinary mesh onto the W3 space of the vertically shifted mesh the
finite-volume transport scheme works on. Each column is a copy of the interior
levels with a 3:1 weighted average at the bottom and top, so the region carries
a short prologue, a level loop and an epilogue over two function spaces and two
dofmaps.

It was chosen because it is the busiest loop C16_MG runs that phase 3 stage 3
unblocks and nothing else blocks: 60 entries over the run, 6 per timestep. It
is not the loop the stage was planned around -- see the divergence recorded
under Task 3.7 of docs/plans/2026-08-29-phase-3-coverage.md in psy-ir-aidev.
tri_solve, at 1872 entries, turned out to hold two automatic local arrays,
which is a second pattern the transformation does not model and this stage does
not add.

It is the first captured region whose kernel is kind-polymorphic.
sci_inject_wt_to_sh_w3_kernel_mod writes inject_wt_to_sh_w3_code as a generic
interface over inject_wt_to_sh_w3_code_r_single and
inject_wt_to_sh_w3_code_r_double, and this algorithm's actual arguments are
r_tran_field_type. LFRic's precision map makes r_tran 8 bytes in this
configuration, so LFRicKokkosTrans resolves the call to the r_double
implementation and the generated region computes in double.

That resolution is the second procedure of the interface, not the first. A
transformation that took whichever schedule PSyclone happened to present first
would capture the r_single body here, compute the same columns in float and
still present a consistent boundary; only the whole-model checksum shows the
difference. The region's bind(C) interface additionally asserts at compile time
that r_tran really is 8 bytes, so a build configured with R_TRAN=32 fails to
compile rather than running two precisions against each other.

This algorithm holds one invoke of the kernel, so unlike ffsl_vert_alg_mod
there is no second call site left uncaptured for a comparison to attribute a
difference to.

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

# This algorithm's only invoke of inject_wt_to_sh_w3, and the one coded kernel
# within it. The name is the one the generated PSy layer uses, read from it
# rather than predicted.
TARGET_INVOKE = 'invoke_1_inject_wt_to_sh_w3_kernel_type'
TARGET_KERNEL = 'inject_wt_to_sh_w3_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
