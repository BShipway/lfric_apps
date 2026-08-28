##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
Diagnostic PSyclone transformation script: find out which coded-kernel loops a
configuration actually executes.

Reading the call graph is not enough. The moisture guard that stopped C16_MG
from ever entering the moist_dyn_gas loop sits three levels above the kernel,
in an algorithm that looks unconditional from the PSy layer down. So this
script measures instead, in two halves:

  * at generation time it writes one row per coded-kernel cell loop into
    WORKING_DIR/loop-inventory/, recording the structural facts that decide
    whether a loop could be captured as a Kokkos region;

  * at run time each of those loops announces itself on standard output with a
    'TRACE-LOOP' line, so the executed set is whatever the run prints.

Intersect the two and the candidates for capture are what is left:

    grep '^TRACE-LOOP' run.out | sort | uniq -c | sort -rn

This is instrumentation, not an optimisation. It exists to be built, run once
and read; nothing in the prototype depends on its output being reproducible
from a build that has it applied. Built-in loops are deliberately not traced:
they are neither capture candidates nor interesting here, and there are enough
of them to bury the output.

'''
import os
from pathlib import Path

from psyclone.domain.lfric import LFRicLoop
from psyclone.domain.lfric.lfric_builtins import LFRicBuiltIn
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.nodes import CodeBlock, Container

from psyclone_tools import redundant_computation_setval

# One file per PSy-layer module rather than one shared file: make runs several
# PSyclone processes at once, and appends to a shared file would interleave.
INVENTORY_NAME = 'loop-inventory'

# The row layout, recorded here rather than written into the files: they are
# meant to be concatenated, and a header per module would be a row per module
# that is not a loop.
COLUMNS = ('module', 'invoke', 'kernels', 'loop_type', 'iteration_space',
           'lower_bound', 'upper_bound', 'halo_depth', 'arguments',
           'quadrature', 'evaluator', 'cma', 'intergrid', 'codeblocks')


def inventory_path(module):
    '''
    :param str module: the PSy-layer module the rows belong to.

    :returns: where to write this module's share of the inventory.
    :rtype: :py:class:`pathlib.Path`

    :raises RuntimeError: if WORKING_DIR is unset.

    '''
    working_dir = os.environ.get('WORKING_DIR')
    if not working_dir:
        raise RuntimeError(
            "WORKING_DIR must name the PSyclone output directory: the loop "
            "inventory has nowhere else to go.")
    return Path(working_dir) / INVENTORY_NAME / f'{module}.tsv'


def describe_argument(argument):
    '''
    Summarises one kernel argument as the transformation would have to see it.

    :param argument: the argument to describe.
    :type argument: :py:class:`psyclone.psyGen.Argument`

    :returns: a compact description.
    :rtype: str

    '''
    space = getattr(argument, 'function_space', None)
    parts = [argument.argument_type, argument.intrinsic_type,
             str(argument.access).rsplit('.', maxsplit=1)[-1]]
    if space is not None:
        parts.append(space.orig_name)
    if getattr(argument, 'stencil', None):
        parts.append('stencil')
    return ':'.join(part for part in parts if part)


def has_codeblocks(kernel):
    '''
    Reports whether a kernel body is beyond PSyIR, which rules out capture.

    A kernel PSyclone cannot parse at all is reported as 'unparsed' rather
    than allowed to stop the build: this script is a survey, and a kernel it
    cannot read is a kernel the transformation could not capture either.

    :param kernel: the coded kernel to inspect.
    :type kernel: :py:class:`psyclone.psyGen.CodedKern`

    :returns: 'yes', 'no' or 'unparsed'.
    :rtype: str

    '''
    try:
        schedules = kernel.get_callees()
    except Exception:                       # pylint: disable=broad-except
        return 'unparsed'
    if len(schedules) != 1:
        return 'unparsed'
    return 'yes' if schedules[0].walk(CodeBlock) else 'no'


def loop_row(module, schedule, loop):
    '''
    Builds the inventory row for one loop.

    :returns: the tab-separated row, without a trailing newline.
    :rtype: str

    '''
    kernels = loop.kernels()
    names = '+'.join(kernel.name.lower() for kernel in kernels)
    arguments = ' '.join(
        describe_argument(argument)
        for kernel in kernels for argument in kernel.arguments.args)
    # pylint: disable=protected-access
    fields = (
        module,
        schedule.name,
        names,
        loop.loop_type or 'none',
        loop.iteration_space,
        loop._lower_bound_name,
        loop.upper_bound_name,
        str(loop.upper_bound_halo_depth),
        arguments,
        str(any(kernel.qr_required for kernel in kernels)),
        str(any(bool(kernel.eval_shapes) for kernel in kernels)),
        str(any(kernel.cma_operation is not None for kernel in kernels)),
        str(any(kernel.is_intergrid for kernel in kernels)),
        ' '.join(has_codeblocks(kernel) for kernel in kernels),
    )
    return '\t'.join(fields)


def trace_statement(schedule, tag):
    '''
    Builds the run-time announcement for one loop.

    :param schedule: the schedule whose symbols the statement is parsed in.
    :type schedule: :py:class:`psyclone.psyGen.InvokeSchedule`
    :param str tag: what the loop should call itself.

    :returns: a PSyIR statement that prints the tag.
    :rtype: :py:class:`psyclone.psyir.nodes.Node`

    '''
    return FortranReader().psyir_from_statement(
        f"write(*,'(A)') 'TRACE-LOOP {tag}'", schedule.symbol_table)


def trans(psyir):
    '''
    Applies the minimum transformations, then inventories and traces every
    coded-kernel loop.

    :param psyir: the PSyIR of the PSy layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)

    rows = []
    module = psyir.name
    for schedule in psyir.walk(InvokeSchedule):
        container = schedule.ancestor(Container)
        module = container.name if container else psyir.name
        for loop in schedule.walk(LFRicLoop):
            kernels = loop.kernels()
            if not kernels or any(isinstance(kernel, LFRicBuiltIn)
                                  for kernel in kernels):
                continue
            names = '+'.join(kernel.name.lower() for kernel in kernels)
            rows.append(loop_row(module, schedule, loop))
            tag = f'{module} {schedule.name} {names}'
            loop.parent.children.insert(
                loop.position, trace_statement(schedule, tag))

    if rows:
        path = inventory_path(module)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('\n'.join(rows) + '\n', encoding='utf-8')
