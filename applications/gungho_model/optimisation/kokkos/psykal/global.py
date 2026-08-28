##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script for the LFRic API applying the same minimum
set of global transformations as the 'minimum' option.

The 'kokkos' option differs from 'minimum' in one algorithm only, through the
local transformation script for 'moist_dyn_factors_alg_mod'. Everything else
must be generated identically, so that a numerical comparison between the two
builds measures the generated Kokkos region and nothing else.

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
