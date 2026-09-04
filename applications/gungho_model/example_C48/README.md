# `example_C48` — the same configuration at C48

This directory exists for one measurement: gate G0 asks whether the share of a
timestep that is cell loops holds as the mesh grows, so it needs a second
resolution that differs from `example/` in the mesh and in **nothing else**.

Select it with `LFRIC_EXAMPLE=example_C48`, which `lfric-run` already reads
and `psy-ir-aidev/bin/run-gungho-model` forwards into the container.

## What differs from `example/`

Seven lines across two files. Everything else is a byte-for-byte copy.

| File | Line | `example/` | here |
|---|---|---|---|
| `C48_MG.nml` | `mesh_file_prefix` | `'mesh_C16_MG'` | `'mesh_C48_MG'` |
| `C48_MG.nml` | `mesh_names` | `'C16','C8','C4','C2'` | `'C48','C24','C12','C6'` |
| `C48_MG.nml` | `mesh_maps` | `'C16:C8','C8:C4','C4:C2'` | `'C48:C24','C24:C12','C12:C6'` |
| `C48_MG.nml` | `edge_cells` | `16,8,4,2` | `48,24,12,6` |
| `configuration.nml` | `&base_mesh file_prefix` | `'mesh_C16_MG'` | `'mesh_C48_MG'` |
| `configuration.nml` | `&base_mesh prime_mesh_name` | `'C16'` | `'C48'` |
| `configuration.nml` | `&multigrid chain_mesh_tags` | `'C16','C8','C4','C2'` | `'C48','C24','C12','C6'` |

The file whose stem the mesh namelist is named for is the mesh it writes, which
is why this one is `C48_MG.nml` rather than `C16_MG.nml` with different
contents. `psy-ir-aidev/bin/generate-mesh` refuses the mismatch by name; its
header says why at length.

## What is deliberately the same

**`&timestepping dt = 3600` and `&time timestep_end = '10'` are unchanged**, and
that is a decision rather than an oversight. Shortening the run at C48 to keep
the wall time down would change the Courant number the model is integrating at,
and a model at a different CFL can take different branches — a different number
of solver iterations, a different limiter path. The two meshes' cell-loop shares
would then differ for two reasons at once and the comparison would say nothing.
Ten timesteps at `dt=3600` costs what it costs.

`&multigrid multigrid_chain_nitems=4` is unchanged because the chain is still
four meshes; only their names and sizes moved. `&partitioning` is unchanged
because `panel_decomposition='auto'` works out the decomposition from the rank
count, which is 1 for this measurement.

## The `.xml` files are copies

All **ten** of them, unchanged:

    file_def_ancil.xml           iodef.xml
    file_def_check_restart.xml   iodef_climate.xml
    file_def_diags_climate.xml   iodef_lam.xml
    file_def_diags_nwp.xml       iodef_nwp.xml
    file_def_diags_lam.xml
    file_def_initial.xml

They are XIOS I/O definitions and none of them names a mesh, so there is nothing
in them to change for a resolution. **An edit to `example/`'s copies needs
mirroring here**, and nothing enforces that — this paragraph is the whole of the
mechanism. It is stated rather than automated because a symlink would break the
property that makes this directory useful, which is that it is a plain copy a
reader can diff.

## The mesh file is generated, not committed

`mesh_C48_MG.nc` is not in the repository. `applications/gungho_model/*/*.nc` is
in `.gitignore`, and this directory's mesh is built by

    psy-ir-aidev/bin/generate-mesh \
        applications/gungho_model/example_C48/C48_MG.nml \
        applications/gungho_model/example_C48

`example/` differs here: its `mesh_C16_MG.nc` *is* committed, which is why the
`.gitignore` rule had to be written to reach this directory without reaching
that one.

**Nothing else the generator writes belongs here.** It finalises through LFRic's
`timer_mod` and so leaves a `timer.txt` in whatever directory it ran in.
`generate-mesh` deletes it, and must: `lfric-run` stages an example with `cp -r`,
so a `timer.txt` sitting here would be copied into every run tree before the
model runs, and a harness looking for the model's report could find the mesh
generator's instead. If one ever appears in this directory, it is a bug in
`generate-mesh` and not something to add to `.gitignore`.

## `qrparm.orog.ugrid.nc` is not copied

`example/` carries an orography ancillary and this directory does not. It is
not needed: `&orography orog_init_option='none'`, and `file_def_ancil.xml`
marks `orography_mean_ancil` `enabled=".FALSE."`, so nothing opens it. Copying
it would have been worse than omitting it — a C16 orography sitting in a C48
directory is a file that looks like it belongs and would be wrong the moment
someone switched the ancillary on.
