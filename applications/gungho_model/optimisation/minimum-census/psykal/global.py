##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
'minimum-timed', with a second caliper set covering every cell loop it leaves.

'minimum-timed' answers 'how long do the captured kernels take'. This answers
'how much of a timestep is cell loops at all', which is the denominator the
first number needs before it means anything. Both come out of one run, because
the two caliper sets partition the cell loops between them:

    cell_loops = region:cell_loops + the eight captured regions

The captured eight are timed under their own names by time_captured_loops, and
census_cell_loops covers everything else under one shared name. Nothing is
timed twice and nothing is missed, so the sum is the total and the ratio of
the parts is what the strategy asks for.

WHY THE CAPTURED CALIPERS GO FIRST

The order of the two calls is not free, and getting it wrong fails the build
rather than mismeasuring it -- but only because census_cell_loops was written
to refuse. time_captured_loops places its calipers on *named* loops: it reads
CAPTURED_REGIONS, finds the kernel each row names and walks up to the loop to
bracket. census_cell_loops places its caliper on whatever is left, and it works
out what is left by asking the same question of the same rows.

Placing the census first would leave _target_loop looking for a loop that is
now the child of a ProfileNode. Whether that resolves to the same node or to
nothing is a detail of the walk rather than a guarantee, and neither answer is
wanted: a captured kernel would end up either inside the shared caliper as well
as its own -- a nested toggle, which subtracts silently -- or outside both.
census_cell_loops raises RuntimeError on a candidate with a ProfileNode
ancestor precisely so that this order cannot be got wrong quietly.

There is no such constraint the other way round. time_captured_loops does not
look at what the census has placed because it never runs after it.

WHY THIS IS A SEPARATE TRANSFORMATION AND NOT A FLAG

'minimum-timed' is what a Kokkos region is compared against, and a comparison
is only worth anything if the thing compared has not moved. Adding 600-odd
calipers to it would change the number every earlier measurement was taken
against. This tree is therefore a sibling of 'minimum-timed', identical to it
apart from the census, so that the overhead of the census can be measured by
running both -- which is what task A6 does.

WHAT 'report' MEANS IN A CENSUS BUILD

timed_region.report walks every ProfileNode in the layer, so once the census
has run it lists the census calipers too and prefixes them all 'timed_region:'.
That prefix is the one place in the build log that does not name the module
which placed the caliper. It is left as it is because the full site list is
worth more than the label: report gives every caliper with its module and
region name, which is what a check of the placement reads, and census_report
gives the per-module totals underneath it.

'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)

# PSyclone puts only this script's own directory on sys.path. Both helpers are
# shared and so live above the transformation trees; the first parent holding
# *both* is used. Requiring them together rather than one at a time is not
# tidiness: cell_loop_census imports CAPTURED_REGIONS, INNER_COLOUR_LOOPS,
# _refuse_inside_openmp and _target_loop from timed_region, so a walk that
# found the two in different directories would have the census deciding what
# is already captured by reading a different timed_region from the one placing
# the captured calipers, and the partition above would silently stop holding.
_HELPER = [p for p in Path(__file__).resolve().parents
           if (p / 'timed_region.py').is_file()
           and (p / 'cell_loop_census.py').is_file()]
if not _HELPER:
    raise RuntimeError(
        f"{__file__} has no directory above it holding both timed_region.py "
        f"and cell_loop_census.py")
sys.path.insert(0, str(_HELPER[0]))

from timed_region import time_captured_loops, report      # noqa: E402
from cell_loop_census import census_cell_loops, census_report    # noqa: E402


def trans(psyir):
    '''
    Applies the minimum transformations, the captured calipers, then the
    census.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    time_captured_loops(psyir)
    census_cell_loops(psyir)
    report(psyir)
    census_report(psyir)
    view_transformed_schedule(psyir)
