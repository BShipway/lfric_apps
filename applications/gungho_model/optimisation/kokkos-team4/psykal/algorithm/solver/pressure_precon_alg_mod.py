##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Local PSyclone transformation script capturing one LFRic cell loop as a
Kokkos region.

The tri_solve kernel applies the vertical part of the pressure
preconditioner: for each column it solves a tridiagonal system by the Thomas
algorithm, sweeping down the column to eliminate the sub-diagonal and back up
to substitute. Both sweeps are sequential in the level index, so the column is
the unit of parallelism and the levels within it are not.

It is the busiest loop C16_MG runs: 1872 entries over the run, 187 per
timestep, against 60 for the next busiest. That is why phase 3 stage 3 was
planned around it and why stage 3a returns to it -- see the divergence under
Task 3.7 of docs/plans/2026-08-29-phase-3-coverage.md in psy-ir-aidev, which
records finding mid-stage that it needed a pattern that stage did not model.

It is the first captured region to need scratch. The kernel holds two
automatic local arrays,

    real(kind=r_single), dimension(nlayers) :: x_new, tri_plus_new

whose extent is a runtime value, so neither can be a C++ local. They are
placed in Kokkos scratch, requested PerThread and private to the team rank
running the cell, which puts the region on a TeamPolicy rather than the
RangePolicy the other three regions captured by then use. Two cells sharing a
team therefore do not share a temporary, which is what makes the parallelism
over columns sound.

Stage 8 added a third launch shape and left this region the only one on the
second. Its two sweeps are recurrences, so
DependencyTools.can_loop_be_parallelised refuses its level loops and the
region keeps the flat team launch: a team of team_size_recommended ranks,
each rank taking one cell and its own team.thread_scratch(0) slice, with an
early return guarding the tail of the last league. The five other captured
regions moved to the hierarchical launch, TeamPolicy(ncells, Kokkos::AUTO)
with the level loops spread as TeamVectorRange. Being refused that shape by
measurement rather than left out of it is what makes this region the control
the launch-shape stage measured against: its generated C++ is byte-identical
across the change.

Like inject_wt_to_sh_w3 this kernel is kind-polymorphic.
sci_tri_solve_kernel_mod writes tri_solve_code as a generic interface over
tri_solve_code_r_single and tri_solve_code_r_double, declared in that order,
and this algorithm's actual arguments are r_solver_field_type.

Unlike inject_wt_to_sh_w3, the member selected here is the interface's
*first*. R_SOLVER_PRECISION is not one of the widths bin/build-gungho-model
sets, so it takes lfric.mk's default of 32; psyclone.cfg's precision map
agrees, recording r_solver as 4 bytes. LFRicKokkosTrans therefore resolves
the call to tri_solve_code_r_single and generates tri_solve_r_single_kokkos,
in float. This is the region that shows the selection rule is a width
comparison rather than a habit of taking the second procedure -- the two
captured polymorphic kernels now select opposite ends of the same two-member
interface. It is also the first captured region generated in single
precision; the other three are double.

The generated bind(C) interface asserts at compile time that r_single really
is 4 bytes, and names r_single -- the selected member's own kind -- rather
than the r_solver the algorithm layer passes. A build that made r_solver 8
bytes would select the r_double member and regenerate, so the assertion
guards the width the C++ was compiled for rather than the width the
algorithm happens to use.

This algorithm holds one invoke of the kernel, so there is no second call site
left uncaptured for a comparison to attribute a difference to. The invoke is
at line 123 of the .x90, inside apply_pressure_preconditioner.

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

# This algorithm's only invoke of tri_solve, and the one coded kernel within
# it. The name is the one the generated PSy layer uses, read from it rather
# than predicted.
TARGET_INVOKE = 'invoke_0_tri_solve_kernel_type'
TARGET_KERNEL = 'tri_solve_code'


def trans(psyir):
    '''
    Applies the global transformations, then captures the target loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    capture(psyir, __file__, TARGET_INVOKE, TARGET_KERNEL)
