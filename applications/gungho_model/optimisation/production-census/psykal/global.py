##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
'production-timed', with a second caliper set covering every cell loop it
leaves.

The coloured, threaded half of the breakdown. 'minimum-census' measures the
loop structure the Kokkos build has; this measures the loop structure a real
run has, and the two are a pair rather than one being a better version of the
other. As there, the two caliper sets partition the loops between them:

    cell_loops = region:cell_loops + the eight captured regions

THE ORDER OF THE FIVE CALLS IS NOT FREE

'production-timed' has four and this has five, and the new one is fixed at both
ends rather than only one.

colour_loops must come before the calipers. It walks each Routine's immediate
children and raises TransformationError -- "Must apply colour_loops BEFORE
profile_loops function in optimisation script" -- on finding a ProfileNode
among them, because a caliper already in place would hide the loop it was meant
to colour.

census_cell_loops must come after colour_loops, for a reason that is the same
constraint read from the other side. On a coloured schedule the loop to time is
the outer half of the colouring, and colour_loops is what creates it: run the
census first and it would bracket the single uncoloured loop, which
colour_loops would then refuse to touch. The census would have quietly moved
this build's loop structure instead of measuring it.

census_cell_loops must also come after time_captured_loops. Both work out which
loop to bracket from CAPTURED_REGIONS, and the census places its caliper on
what is left over; placing the census first would leave _target_loop looking
for a loop that is now inside a ProfileNode. census_cell_loops raises
RuntimeError on a candidate with a ProfileNode ancestor so that the wrong order
fails the build rather than nesting a toggle, which subtracts silently.

The calipers -- both sets -- must come before openmp_parallelise_loops.
timed_region refuses a caliper with an OpenMP ancestor for the reason
psyclone_tools does: PSyData is not thread-safe and timer_mod's table is module
state toggled by name, so a caliper entered by every thread corrupts every row.
census_cell_loops calls the same refusal on every candidate, so this half of
the constraint is enforced by the same code and with the same wording.

Every one of these refusals names the function to move, so a wrong order fails
the build rather than quietly mismeasuring it.

WHY THE CENSUS COUNT NEED NOT MATCH 'minimum-census's

The captured eight bracket the same work here as there, because _target_loop
returns the outer half of a colouring and skips the inner. The census does the
same -- INNER_COLOUR_LOOPS is what it filters on -- but that does not make the
two counts equal, because colouring does not merely rename loops. A loop over
cells becomes a loop over colours containing a loop over cells in colour, and
what the census counts here is the outer of those. Discontinuous spaces are not
coloured at all and are counted unchanged.

So a difference between this build's census total and 'minimum-census's is
information about how much of the model colouring reaches, and not a defect.
Both numbers are recorded with that said beside them.

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

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops,
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
    Applies redundant computation, colouring, both caliper sets and OpenMP.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    time_captured_loops(psyir)
    census_cell_loops(psyir)
    openmp_parallelise_loops(psyir)
    report(psyir)
    census_report(psyir)
    view_transformed_schedule(psyir)
