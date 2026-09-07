##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script for the LFRic API capturing every cell loop
LFRicKokkosTrans accepts as a Kokkos region.

The 'kokkos-all' option applies the same minimum set of global transformations
as 'minimum' and 'kokkos', then captures every loop the transformation's
validate() passes -- in every invoke of every algorithm, because this is the
global script and there are no local ones beside it. 'kokkos' captures eight
named loops through eight local scripts; this captures whatever the
transformation can, so that a whole-model checksum comparison against
'minimum' measures every generated region at once.

The capture itself is capture_all.py at the root of optimisation/, shared
with 'kokkos-all-timed' so that the two trees cannot drift, and found by
walking this file's parents in the same way the timed trees find
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
    Applies the global transformations, then captures every accepted loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    view_transformed_schedule(psyir)
    capture(psyir)
