##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script applying colouring, OpenMP and redundant
computation, for the algorithms that no Kokkos region is captured from.

SHARED VERBATIM BY 'kokkos-omp' AND 'kokkos-omp-timed'

The two trees hold byte-identical copies of this file, because the calipers
that separate them go where the capture is and this script never meets a
captured algorithm -- see kokkos_region.py, and see 'kokkos-timed's global.py
for the same argument made about 'kokkos'. So

    diff -r optimisation/kokkos-omp/psykal optimisation/kokkos-omp-timed/psykal

reports kokkos_region.py alone, and a change made here belongs in both copies.

'kokkos-omp' IS THE EXPERIMENT AND 'kokkos-omp-timed' IS THE MEASUREMENT

Can a Kokkos region and LFRic's own colouring and OpenMP be applied to the
same build? Every production profile colours everything, and
LFRicKokkosTrans._validate_loop refuses a coloured loop outright, so if the
two cannot coexist then no production LFRic configuration can carry a Kokkos
region at all -- however good the region is.

The question is structural rather than numerical, and it is asked here, on a
host, deliberately: finding the answer during GPU bring-up would mean debugging
it with a device in the way.

'kokkos-omp' must stay untimed to answer it, because its checksums are the
evidence and they are compared against 'meto-ex1a' built untimed. The timed
copy stands to it exactly as 'production-timed' stands to 'meto-ex1a'.

WHY THIS FILE IS 'meto-ex1a's AND NOT 'kokkos's

'kokkos' applies redundant computation and nothing else, so that a comparison
against 'minimum' measures the generated regions alone. This transformation is
the opposite comparison: it is 'production-timed' with the Kokkos captures
added, so its uncaptured algorithms must be transformed exactly as that
build's are. The three calls below are 'meto-ex1a's, which is what
'production-timed' is built from.

The seven captured algorithms do not use this file. Each has a local script
beside it, and psyclone_psykal.mk uses one or the other and never both -- see
kokkos_region.py, which repeats these transformations after the capture.

'''
from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops,
                            view_transformed_schedule)


def trans(psyir):
    '''
    Applies redundant computation, colouring and OpenMP.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    openmp_parallelise_loops(psyir)
    view_transformed_schedule(psyir)
