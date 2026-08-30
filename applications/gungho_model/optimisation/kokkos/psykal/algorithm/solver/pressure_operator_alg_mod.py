##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The apply_helmholtz_operator kernel applies the Helmholtz operator that the
mixed-solver preconditioner inverts. For each column it forms a nine-point
sum: the central coefficient, the four horizontal neighbours reached through
a CROSS2D stencil, and the two levels above and two below reached through the
level index. The horizontal part is what makes this kernel unlike the five
captured before it -- it reads dofs belonging to other cells.

It is the second busiest loop C16_MG runs: 1755 entries over the run against
tri_solve's 1872, and every other loop is far behind both.

It is the first captured region to read a stencil dofmap. LFRic gives the
kernel three formals for that stencil,

    integer(kind=i_def), dimension(4), intent(in)                  :: smap_sizes
    integer(kind=i_def), intent(in)                                :: max_length
    integer(kind=i_def), dimension(ndf, max_length, 4), intent(in) :: smap

and the region needs no argument machinery of its own for them. The two
arrays are cell-sliced actuals in the PSy layer, so they become Views with
the cell index appended, exactly as the ordinary dofmap map_w3(:,cell) does;
max_length is a scalar and stays one. A 1-D stencil would not work this way,
which is why LFRicKokkosTrans accepts cross2d and refuses the other shapes by
name.

A stencil also makes the PSy layer emit a halo exchange in front of the loop.
That exchange computes its own depth by walking forward to the accesses that
read the field, so it is lowered before the loop is replaced rather than
after -- replacing the loop first removes the very accesses it walks to. This
is the first captured region where that ordering matters.

It is also the first captured region to take a logical argument. lam_mesh is
a logical(l_def), and it crosses the ABI by conversion rather than by width:
the interface declares logical(c_bool), value and the call site wraps the
actual in LOGICAL(..., c_bool). No width assertion is generated for it,
because there is none to make -- which matters here, since PSyclone's
precision map records l_def as 1 byte where LFRic defines it as kind(.false.)
and it measures 4. See PSyclone issue #1941.

The kernel holds one automatic local array,

    real(kind=r_solver), dimension(max_length,4) :: coeff

whose extent is a runtime value, so it is placed in Kokkos team scratch and
the region is launched over a TeamPolicy rather than a RangePolicy, as
tri_solve's is.

The kernel declares 'cell' as a local of its own and counts stencil branches
with it. That is the name the generated launch gives its own cell index, and
both land in the same C++ block, so the region renames its index rather than
leaving the collision to the compiler -- which would reject the translation
unit outright.

It is generated in single precision. R_SOLVER_PRECISION is not one of the
widths bin/build-gungho-model sets, so it takes lfric.mk's default of 32 and
r_solver reaches C++ as float. Unlike tri_solve this kernel is a single
procedure rather than a generic interface, so no member selection is
involved.

This algorithm holds one invoke of the kernel, the named apply_h at line 188
of the .x90, so there is no second call site left uncaptured for a comparison
to attribute a difference to.

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

# This algorithm's only invoke of apply_helmholtz_operator, and the one coded
# kernel within it. The invoke is named in the .x90, so the generated schedule
# takes that name rather than a numbered one.
TARGET_INVOKE = 'invoke_apply_h'
TARGET_KERNEL = 'apply_helmholtz_operator_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
