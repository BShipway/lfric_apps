##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
The coloured, OpenMP-threaded baseline, with a timing caliper round each of
the captured kernels.

This is 'meto-ex1a' -- redundant computation, colouring, OpenMP -- with the
calipers added. It answers the question 'minimum-timed' cannot: a Kokkos region
that beats an unthreaded Fortran loop has not yet been shown to beat what LFRic
would actually run on this host.

THE ORDER OF THE FOUR CALLS IS NOT FREE

colour_loops must come before the calipers. It walks each Routine's immediate
children and raises TransformationError -- "Must apply colour_loops BEFORE
profile_loops function in optimisation script" -- on finding a ProfileNode
among them, because a caliper already in place would hide the loop it was meant
to colour.

The calipers must come before openmp_parallelise_loops. timed_region refuses a
caliper with an OpenMP ancestor for the reason psyclone_tools does: PSyData is
not thread-safe and timer_mod's table is module state toggled by name, so a
caliper entered by every thread corrupts every row.

Both refusals name the function to move, so a wrong order fails the build
rather than quietly mismeasuring it.

WHY THE CALIPER SET IS THE SAME AS 'minimum-timed's

timed_region times the outermost loop of a colouring and skips its inner half,
so a coloured build and an uncoloured build bracket the same work under the
same seven names. That is what lets a row from this build sit beside a row from
'minimum-timed' and a row from 'kokkos-timed' in one table.

This is *not* what profile_loops does with its default colours_only=True: that
times a loop because colouring reached it, so a kernel on a discontinuous space
would be timed here and absent from the uncoloured baseline. The three timed
transformations have to report the same region set or the table has holes in
it, so the set is fixed by CAPTURED_REGIONS and not by what any one
transformation happens to produce.

'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops,
                            view_transformed_schedule)

# PSyclone puts only this script's own directory on sys.path. The caliper
# helper is shared by all three timed transformations and so lives above them;
# the first parent holding it is used, so that a copy placed inside a
# transformation tree would be found ahead of the shared one.
_HELPER = [p for p in Path(__file__).resolve().parents
           if (p / 'timed_region.py').is_file()]
if not _HELPER:
    raise RuntimeError(f"{__file__} has no timed_region.py above it")
sys.path.insert(0, str(_HELPER[0]))

from timed_region import time_captured_loops, report    # noqa: E402


def trans(psyir):
    '''
    Applies redundant computation, colouring, the calipers and OpenMP.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    time_captured_loops(psyir)
    openmp_parallelise_loops(psyir)
    report(psyir)
    view_transformed_schedule(psyir)
