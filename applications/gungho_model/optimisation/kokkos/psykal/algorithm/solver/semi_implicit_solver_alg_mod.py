##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The dg_matrix_vector kernel applies an LMA operator to a field on a
discontinuous space, lhs = matrix * x, column by column: it zeroes the
column's slice of lhs over the output space's dofs, then for each dof of the
input space and each dof of the output accumulates the operator's column
block against the input. Here the operator is the divergence, the input is
f_star -- the density-weighted wind increment the invoke's first built-in
forms -- and the output is div_u, which the invoke's last built-in folds into
the density increment. It runs 40 times over a 10-timestep C16_MG run, four
per timestep.

It is the pattern stage 10 cleared -- an LMA operator reaching the kernel as
an integer ncell_3d and a rank-3 array over (ncell_3d, ndf_to, ndf_from) --
in its simplest form: one operator, one input field, one output field, and no
local arrays, so nothing goes to scratch. The kernel writes its level
arithmetic as array sections, lhs(i1:i1+nl) and matrix(ik:ik+nl, df1, df2),
which LFRicKokkosTrans lowers to explicit level loops before generating; the
region therefore takes the hierarchical launch, TeamPolicy(ncells,
Kokkos::AUTO) with the lowered level loops spread as TeamVectorRange.

It is generated in single precision. The kernel is kind-polymorphic over
r_single and r_double, and the algorithm's actuals are r_solver_field_type;
R_SOLVER_PRECISION takes lfric.mk's default of 32, so the r_single member is
selected, as it is for tri_solve and apply_mixed_u_operator.

This is not the model's only call site of the kernel: diagnostic_alg_mod
invokes it twice over the run, from invoke_0_dg_matrix_vector_kernel_type,
and that site stays Fortran. The captured site is the invoke at lines 711-714
of the .x90, the middle of three steps that compute the rho increment the
solver does not provide, and it is selected by (module, invoke, kernel) so
that the timed builds' caliper brackets exactly these forty executions and
not the diagnostic's two.

Unlike the other captured invokes this one has no kernel in its name.
PSyclone names an invoke after its kernel only when the invoke holds one;
invoke_11 holds a kernel between two built-ins, so it takes its ordinal. The
name is read from the generated PSy layer, not predicted.

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

# This algorithm's only invoke of dg_matrix_vector, and the one coded kernel
# within it. The name is the one the generated PSy layer uses, read from it
# rather than predicted.
TARGET_INVOKE = 'invoke_11'
TARGET_KERNEL = 'dg_matrix_vector_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
