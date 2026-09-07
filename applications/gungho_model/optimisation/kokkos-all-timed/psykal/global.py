##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
'kokkos-all', with the nine known regions calipered.

This is 'kokkos-all' plus one flag: capture(psyir, timed=True) places a
caliper on each captured call whose site is in timed_region.CAPTURED_REGIONS,
under the region name the 'minimum-timed' and 'kokkos-timed' builds use, so
the three timed builds join on those rows. Everything else captured here runs
uncalipered, for the reason timed_region.py's docstring gives: timer_mod holds
300 timers and a caliper on every kernel runs it out. What this build adds to
the comparison is the timestep itself, from LFRic's own timers, with the whole
model captured rather than eight loops of it.

The capture module is shared with 'kokkos-all' and found the same way, and
the schedule is viewed before it is captured for the reason 'kokkos-all's
global.py gives: the viewer reads LFRicLoop attributes, and a capture leaves
its invoke's sibling loops lowered.
'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)

# PSyclone puts only this script's own directory on sys.path. Both helpers are
# shared and live above the transformation trees; the first parent holding
# both is used, so the capture and the calipers come from one directory.
_HELPER = [p for p in Path(__file__).resolve().parents
           if (p / 'capture_all.py').is_file()
           and (p / 'timed_region.py').is_file()]
if not _HELPER:
    raise RuntimeError(
        f"{__file__} has no directory above it holding both capture_all.py "
        "and timed_region.py")
sys.path.insert(0, str(_HELPER[0]))

from capture_all import capture             # noqa: E402  needs the path above


def trans(psyir):
    '''
    Applies the global transformations, captures every accepted loop, and
    calipers the captured calls the other timed builds also bracket.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    view_transformed_schedule(psyir)
    capture(psyir, timed=True)
