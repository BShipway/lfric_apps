##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script for the LFRic API capturing every cell loop
LFRicKokkosTrans accepts as a Kokkos region, colouring first where a write is
shared.

This is 'kokkos-all' plus one flag: capture(psyir, coloured=True) applies
LFRicColourTrans to a loop whose kernel increments a field before capturing
it, so that the generated region takes the coloured answer to a write two
cells share rather than the atomic one. Cells of one colour meet at no dof,
so the update is a plain read-modify-write and the caller enters the region
once per colour.

WHAT THIS TREE CAPTURES THAT 'kokkos-all' DOES NOT

An atomic answers a read-modify-write of one element and nothing else, so
LFRicKokkosTrans refuses a kernel that updates a shared field with a
whole-array expression, or with any other statement no single atomic carries
out, and says 'Colour the loop instead'. Colouring is that instruction taken:
it puts the cells that meet at a dof in different launches, so the shape of
the update stops being a question. This tree therefore asks the coloured arm
first for a loop whose kernel has a shared write, and falls back to the atomic
arm where the coloured one is refused. That ordering is what lets it capture
those call sites; 'kokkos-all' asks the atomic arm only and leaves them as
Fortran. The two trees no longer capture the same sites, and the difference is
exactly the loops colouring alone can take.

'kokkos-all' is deliberately not changed to match. Its checksums are
bit-identical to the Fortran reference and it is the control every other
profile is measured against, so a capture that widened it would move the
control and the comparison with it. Colouring changes the order the
contributions to a shared dof are summed, which is a floating-point
difference: expected here, and a regression there.

Both answers are correct, so this tree exists to be compared: what
'kokkos-all' and 'kokkos-all-coloured' both capture they must capture to the
same checksums as 'minimum', by two different routes through the same regions.

The capture itself is capture_all.py at the root of optimisation/, shared
with 'kokkos-all' and 'kokkos-all-timed' so that the trees cannot drift, and
found by walking this file's parents in the same way the timed trees find
timed_region.py.

THE VIEW COMES BEFORE THE CAPTURE

psyclone_tools.view_transformed_schedule walks every loop of every invoke and
reads LFRicLoop attributes off it. Capturing a loop lowers the halo exchanges
of its invoke first, and lowering leaves sibling loops as plain PSyIR Loops
with no iteration_space, so the viewer run after a capture raises
AttributeError on the first invoke that captured anything beside a loop it did
not. The 'kokkos' profile never met this: its local scripts capture and do not
view, and its global script views algorithms that captured nothing. So the
view is taken first, of the schedule after the global transformations and
before the capture -- which is the same view 'kokkos' prints for the
algorithms it leaves alone -- and capture() prints one line per region after
it. The lowering is not a defect in the capture: it is what the whole
container undergoes at code generation anyway, and the eight regions the
'kokkos' profile captures reproduce the model's checksums with it.
'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)

# PSyclone puts only this script's own directory on sys.path. The capture
# module is shared and lives above the transformation trees.
_HELPER = [p for p in Path(__file__).resolve().parents
           if (p / 'capture_all.py').is_file()]
if not _HELPER:
    raise RuntimeError(f"{__file__} has no capture_all.py above it")
sys.path.insert(0, str(_HELPER[0]))

from capture_all import capture             # noqa: E402  needs the path above


def trans(psyir):
    '''
    Applies the global transformations, colours the shared writes, then
    captures every accepted loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    view_transformed_schedule(psyir)
    capture(psyir, coloured=True)
