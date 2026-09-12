##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
The Fortran model with every cell loop calipered under one of two names:
'region:captured' for a loop the whole-model Kokkos capture takes, and
'region:residue' for one it leaves.

WHAT THIS MEASURES

Time-weighted coverage. Phase 5 and 6 counted coverage in loop ENTRIES --
94% of the entries C16_MG executes are captured -- and phase 7's first device
timings (2026-09-12) found the device model spending three quarters of a
Fortran step in the host's uncaptured work all the same. The two are
consistent only if the loops left behind are the expensive ones, and the
320-entry FFSL transport residue is the heaviest arithmetic in the model. This
tree measures it rather than infers it: run once, timer.txt holds two rows,
and

    captured / (captured + residue)

is the share of cell-loop time the capture can move to a device, which bounds
what any device can save (Amdahl's law with the residue as the serial part).

WHICH LOOPS ARE WHICH

The manifest beside this script, manifest.tsv, is bin/capture-manifest's
listing of the kokkos-all build's call sites -- psy module, invoke, kernel --
regenerated from a build of the capture profile:

    psy-ir-aidev/bin/capture-manifest scratch/work/kokkos-all-production \
        > lfric_apps/applications/gungho_model/optimisation/residue-census/manifest.tsv

so the census and the capture cannot disagree about a site except by the
manifest being stale, and its header names the revision it was taken from. A
loop is named by the module and invoke it sits in and by its first coded
kernel, which is how capture_all.py's fragments name a site.

WHAT IT DOES NOT MEASURE

Nothing here runs a Kokkos region: the loops are timed as the Fortran runs
them, on a build with the minimum transformations, so the split is what the
capture would take from a Fortran step, not what the device then costs. The
two SKIP sites (the XIOS-setup conversions) are in the manifest as captured
and are two of the cheapest loops in the model.

The census refuses to place a caliper inside another or inside an OpenMP
region, for the reasons cell_loop_census gives, so the two rows partition
the cell loops exactly or the build fails.
'''
import sys
from pathlib import Path

from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)

_HERE = Path(__file__).resolve().parent
_HELPER = [p for p in _HERE.parents
           if (p / 'timed_region.py').is_file()
           and (p / 'cell_loop_census.py').is_file()]
if not _HELPER:
    raise RuntimeError(
        f"{__file__} has no directory above it holding both timed_region.py "
        f"and cell_loop_census.py")
sys.path.insert(0, str(_HELPER[0]))
from cell_loop_census import (census_by_capture,           # noqa: E402
                              read_capture_manifest)

MANIFEST = _HERE.parent / 'manifest.tsv'


def trans(psyir):
    '''
    Applies the minimum transformations, then the two-name census.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    captured, residue = census_by_capture(psyir, read_capture_manifest(MANIFEST))
    print(f"residue-census: {captured} loop(s) under region:captured, "
          f"{residue} under region:residue")
    view_transformed_schedule(psyir)
