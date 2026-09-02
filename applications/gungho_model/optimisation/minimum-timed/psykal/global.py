##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
'minimum', with a timing caliper round each of the captured call sites.

The Fortran baseline a Kokkos region is compared against. It is 'minimum' and
not 'meto-ex1a' on purpose: 'minimum' is what the Kokkos build is apart from
the capture, so a difference between the two is the capture rather than the
colouring and threading a production build would also bring. The coloured,
threaded comparison is 'production-timed', and it is a second measurement
rather than a better version of this one.

Unlike the Kokkos tree this has no per-algorithm scripts. It does not need
them: nothing here is transformed differently from anywhere else, and
time_captured_loops finds the seven call sites from CAPTURED_REGIONS whichever
algorithm PSyclone is generating when it meets one.

The calipers are placed by timed_region rather than by
psyclone_tools.profile_loops. profile_loops with its default colours_only=True
would place none at all here -- this build has no coloured loops -- and would
produce an empty timer report that looks exactly like a successful run.
colours_only=False would place them, but under names carrying a kernel ordinal
that shifts the moment a loop is captured, and on every coded kernel in the
application rather than the seven being compared. See timed_region's docstring
for both, and for why the second is fatal rather than merely wasteful.

'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval,
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
    Applies the minimum transformations, then the timing calipers.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    time_captured_loops(psyir)
    report(psyir)
    view_transformed_schedule(psyir)
