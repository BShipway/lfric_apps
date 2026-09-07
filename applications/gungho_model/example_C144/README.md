# `example_C144` — the same configuration at C144, at a fifteen-minute step

This directory exists for one measurement. Stage 13 put eight captured LFRic
regions on an H100 and found five of them launch-bound: most of their device
time was a fixed launch floor that does not grow with the problem. C16 is too
small to say whether the card is worth using, and C144 — 124416 columns, eighty-
one times C16 — is the mesh that settles it.

Select it with `LFRIC_EXAMPLE=example_C144`, which `lfric-run` already reads and
`psy-ir-aidev/bin/run-gungho-model` forwards into the container.

## What differs from `example_C48`

**Eight lines across two files.** `example_C48` is the template rather than
`example/`, because C48 already made every edit a resolution needs. Everything
else is a byte-for-byte copy.

| File | Line | `example_C48` | here |
|---|---|---|---|
| `C144_MG.nml` | `mesh_file_prefix` | `'mesh_C48_MG'` | `'mesh_C144_MG'` |
| `C144_MG.nml` | `mesh_names` | `'C48','C24','C12','C6'` | `'C144','C72','C36','C18'` |
| `C144_MG.nml` | `mesh_maps` | `'C48:C24','C24:C12','C12:C6'` | `'C144:C72','C72:C36','C36:C18'` |
| `C144_MG.nml` | `edge_cells` | `48,24,12,6` | `144,72,36,18` |
| `configuration.nml` | `&base_mesh file_prefix` | `'mesh_C48_MG'` | `'mesh_C144_MG'` |
| `configuration.nml` | `&base_mesh prime_mesh_name` | `'C48'` | `'C144'` |
| `configuration.nml` | `&multigrid chain_mesh_tags` | `'C48','C24','C12','C6'` | `'C144','C72','C36','C18'` |
| `configuration.nml` | `&timestepping dt` | `3600` | **`900`** |

The file whose stem the mesh namelist is named for is the mesh it writes, which
is why this one is `C144_MG.nml`. `psy-ir-aidev/bin/generate-mesh` refuses the
mismatch by name; its header says why at length.

## The eighth line is a decision, and `example_C48`'s README argues against it

**`dt = 900`, set by the owner on 2026-09-07.** Read `example_C48/README.md`
first and you will find the opposite case made at length: that `dt` is held at
3600 so two meshes differ only in the mesh, because a model at a different
Courant number can take different branches — a different number of solver
iterations, a different limiter path — and the two meshes' figures would then
differ for two reasons at once.

That argument is right, and it stops being decisive here. It buys comparability
between C16 and C48 at a Courant number both can integrate. C144's cells are a
third of C48's across, so `dt = 3600` would be integrating at three times C48's
Courant number, and comparability with a run that does not stay up is worth
nothing. A fifteen-minute step is the owner's judgement about what the model can
integrate at 69 km.

So this directory is deliberately not the minimal edit, and a reader who finds
only C48's argument should not read the difference as an oversight.

**`&time timestep_end = '10'` is unchanged.** Ten timesteps at `dt=900` is 2.5
hours of model time where C48's ten were ten hours. That is not a problem for
what this directory is for: the extraction captures each region's dataset on its
**first** invocation, and the ten steps exist to show the model integrating
stably and to store a checksum baseline, not to reach a particular model time.

`&multigrid multigrid_chain_nitems=4` is unchanged because the chain is still
four meshes; only their names and sizes moved. The chain is a fixed depth, so it
scales with the base mesh: C144,C72,C36,C18, whose coarsest level is 1944
columns where C48_MG's was 216 and C16_MG's 24. `&partitioning` is unchanged
because `panel_decomposition='auto'` works the decomposition out from the rank
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
in them to change for a resolution. **An edit to `example/`'s copies now needs
mirroring into two directories, this one and `example_C48`**, and nothing
enforces that — this paragraph is the whole of the mechanism. It is stated
rather than automated because a symlink would break the property that makes
these directories useful, which is that they are plain copies a reader can diff.

## The mesh file is generated, not committed

`mesh_C144_MG.nc` is not in the repository. `applications/gungho_model/*/*.nc`
is in `.gitignore` for this directory, and the mesh is built by

    psy-ir-aidev/bin/generate-mesh \
        applications/gungho_model/example_C144/C144_MG.nml \
        applications/gungho_model/example_C144

It is about 17 MB, nine times `example_C48`'s 1920772 bytes. `example/` differs
here: its `mesh_C16_MG.nc` *is* committed, which is why the `.gitignore` rules
had to be written to reach this directory and C48's without reaching that one.

**Nothing else the generator writes belongs here.** It finalises through LFRic's
`timer_mod` and so leaves a `timer.txt` in whatever directory it ran in.
`generate-mesh` deletes it, and must: `lfric-run` stages an example with `cp -r`,
so a `timer.txt` sitting here would be copied into every run tree before the
model runs, and a harness looking for the model's report could find the mesh
generator's instead. If one ever appears in this directory, it is a bug in
`generate-mesh` and not something to add to `.gitignore`.

## `qrparm.orog.ugrid.nc` is not copied

`example/` carries an orography ancillary and this directory does not, as
`example_C48` does not. It is not needed: `&orography orog_init_option='none'`,
and `file_def_ancil.xml` marks `orography_mean_ancil` `enabled=".FALSE."`, so
nothing opens it. Copying it would have been worse than omitting it — a C16
orography sitting in a C144 directory is a file that looks like it belongs and
would be wrong the moment someone switched the ancillary on.
