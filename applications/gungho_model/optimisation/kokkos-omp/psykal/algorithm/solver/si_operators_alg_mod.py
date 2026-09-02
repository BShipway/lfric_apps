##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The sample_w3_to_wtheta kernel samples a W3 field onto Wtheta. The interior
levels are a height-weighted linear interpolation between the two W3 levels
either side; the topmost Wtheta level is extrapolated in the logarithm of the
field rather than in the field, so that a pressure-like quantity stays
positive above the last W3 level.

It is the region phase 3 stage 4 adds, and it is the first that needs the C
writer to translate a mathematical function at all. The extrapolation is

    log_top_value = weight_upper / weight_denom * log(field_w3(...)) &
                  + weight_lower / weight_denom * log(field_w3(...))
    top_value     = exp(log_top_value)

so the region does not exist until log and exp can be written, which before
this stage they could not. KokkosWriter emits Kokkos::log and Kokkos::exp,
qualified so that the device overload is chosen rather than the host one.

The second half of the stage is what makes the kernel's casts safe. Every
weight is formed as

    weight_denom = real(height_w3(...) - height_w3(...), r_single)

a two-argument REAL whose kind argument this stage stopped discarding. It is
load-bearing here: height_wt and height_w3 are r_def fields and field_wt and
field_w3 are r_single ones, so this region carries two widths across the same
ABI and the cast is a genuine narrowing of double to float. A writer that
took the cast target from its own type map would have written double and
changed where the rounding happens.

Like inject_wt_to_sh_w3 and tri_solve this kernel is kind-polymorphic.
sci_sample_w3_to_wtheta_kernel_mod writes sample_w3_to_wtheta_code as a
generic interface over sample_w3_to_wtheta_code_r_single and
sample_w3_to_wtheta_code_r_double, and this algorithm's field arguments are
r_solver_field_type, which is real32 in this configuration. The transformation
therefore selects the interface's first member and generates
sample_w3_to_wtheta_r_single_kokkos.

The invoke is guarded by element_order_h == 0 .and. element_order_v == 0, so
it is entered only in the lowest-order configuration. C16_MG is one, and runs
it 40 times: 4 per timestep over 10 timesteps.

This algorithm holds one invoke of the kernel, so there is no second call site
left uncaptured for a comparison to attribute a difference to. The invoke is
at line 546 of the .x90, inside si_operators_alg_init.

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

# This algorithm's only invoke of sample_w3_to_wtheta, and the one coded
# kernel within it. The name is the one the generated PSy layer uses, read
# from it rather than predicted.
TARGET_INVOKE = 'invoke_6_sample_w3_to_wtheta_kernel_type'
TARGET_KERNEL = 'sample_w3_to_wtheta_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
