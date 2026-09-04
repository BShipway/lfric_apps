##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The apply_mixed_u_operator kernel forms the horizontal-velocity row of the
semi-implicit left-hand side,

    lhs_uv = norm_u * (Mu*u - grad*p),

for each column: one loop over the horizontal dofs applies the gradient to
the pressure increment, and two nested loops accumulate the generalised mass
matrix against the horizontal and vertical wind increments. It runs 117 times
over a 10-timestep C16_MG run.

It is the first captured region to take an LMA operator, and it takes TWO --
Mu on W2 x W2 and grad on W2 x W3. Each reaches the kernel as an integer
ncell_3d and a rank-3 array over (ncell_3d, ndf_from, ndf_to), so the region
needs no argument machinery of its own for either: both extents are formals
already, and the two operators become Views of different shape,

    Kokkos::View<const float***, ...> mu_cd(mu_cd_data, ncell1, ndf_w2, ndf_w2)
    Kokkos::View<const float***, ...> grad(grad_data, ncell2, ndf_w2, ndf_w3)

from actuals the PSy layer already had -- mm_vel_proxy%ncell_3d with
mm_vel_local_stencil, and div_star_proxy%ncell_3d with div_star_local_stencil.

Taking an operator also gives the kernel a leading 'cell' argument, and the
generated interface does not carry it. The PSy layer fills that argument with
the cell-column loop's own counter, which is the region's launch index, so
'cell' is declared in the launch body from the index the launch already has:

    const int cell_1 = team.league_rank();
    const int cell = cell_1 + 1;

The interface is therefore one argument shorter than the kernel's signature,
while the arithmetic the kernel writes over it -- ij = (cell-1)*nlayers + 1,
by which it finds its column's slice of both operators -- is generated
unchanged. This is also why the region's own index is named cell_1 rather than
cell: the collision is with the cell position it declares, where the
apply_helmholtz_operator region's was with a local the kernel declared for
itself. Same renaming, a different thing to avoid.

It is generated in single precision. R_SOLVER_PRECISION is not one of the
widths bin/build-gungho-model sets, so it takes lfric.mk's default of 32 and
r_solver reaches C++ as float, as it does in the two solver regions before it.

The kernel holds no local arrays, so nothing goes to scratch, and the region
takes the hierarchical launch on the ordinary ground: its level loops carry no
recurrence. TeamPolicy(ncells, Kokkos::AUTO), one team per cell, each of the
three level loops spread as a TeamVectorRange with a team_barrier after it.
The barriers are not decoration here -- the second and third loops read the
lhs_uv the first wrote.

The lhs is on W2broken, written GH_WRITE. That is what makes the loop
capturable at all: the same computation on continuous W2 would be a shared
write over cells, which is why the assemble step that follows it in the invoke
is not captured.

This algorithm holds one invoke of apply_mixed_u_operator_kernel_type, the
named apply_split_mixed_operator at line 270 of the .x90, so there is no
second call site left uncaptured for a comparison to attribute a difference
to. Within that invoke it is a different story, and this is the first captured
region for which it is: TWO coded kernels stay in Fortran beside it,

    assemble_w2h_from_w2hb_code   halo-depth, shared-write, continuous-write
    apply_mixed_wp_operator_code  array-bound, local-array

along with a setval_c built-in. Neither is blocked by taking an operator --
apply_mixed_wp_operator_code takes five of them -- so neither is waiting on
this stage. The invoke therefore runs one Kokkos region and two Fortran cell
loops in sequence, which is the arrangement a real incremental port would
have, and it is worth knowing that the generated PSy layer expresses it
without complaint.

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

# This algorithm's only invoke of apply_mixed_u_operator, and the first of the
# three coded kernels within it. The invoke is named in the .x90, so the
# generated schedule takes that name rather than a numbered one.
TARGET_INVOKE = 'invoke_apply_split_mixed_operator'
TARGET_KERNEL = 'apply_mixed_u_operator_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
