##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script for the LFRic API applying the same minimum
set of global transformations as the 'minimum' option.

The 'kokkos-timed' option differs from 'minimum-timed' only in the algorithms
that have a local transformation script beside this one. Everything else must
be generated identically, so that a numerical comparison between the two builds
measures the generated Kokkos regions and nothing else.

THIS SCRIPT PLACES NO CALIPERS, WHERE 'minimum-timed's DOES

Every entry in timed_region.CAPTURED_REGIONS names an algorithm that has a
local script in this tree, and psyclone_psykal.mk uses the local script or the
global one for a given algorithm, never both. So this script is never the one
generating a timed call site, and calling time_captured_loops here would search
for loops that by construction cannot be in front of it.

That asymmetry with 'minimum-timed' is the point rather than an oversight:
there the seven call sites are still loops and the global script is what meets
them; here they are Calls placed by the local scripts, which caliper them as
they go. Both trees end up bracketing the same seven pieces of work under the
same seven names, by different routes because the builds differ.

'''
from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)


def trans(psyir):
    '''
    Applies PSyclone redundant computation transformations on
    initialisation built-ins only.

    :param psyir: the PSyIR of the PSyIR-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    view_transformed_schedule(psyir)
