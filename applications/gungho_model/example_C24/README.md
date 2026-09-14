# `example_C24` — the same configuration at C24

The phase-7 campaign (Stream C of
`psy-ir-aidev/docs/plans/2026-09-12-phase-7-h100-host-campaign.md`) compares
three arms over a twelve-day deep baroclinic wave at C24, C48 and C144. C48
and C144 already have a directory each; this is the third, and it exists so
the cheapest rung of the ladder can be run in minutes rather than hours —
sizing runs, the comparator's own gate, and the full rank × thread grid that
would be unaffordable to debug at C48.

Select it with `LFRIC_EXAMPLE=example_C24`, which `lfric-run` already reads
and `psy-ir-aidev/bin/run-gungho-model` forwards into the container.

## What differs from `example_C48`

Four lines across two files. Everything else is a byte-for-byte copy, exactly
as `example_C48` is a byte-for-byte copy of `example/`.

| File | Line | `example_C48` | here |
|---|---|---|---|
| `C24_MG.nml` | `mesh_file_prefix` | `'mesh_C48_MG'` | `'mesh_C24_MG'` |
| `C24_MG.nml` | `mesh_names` | `'C48','C24','C12','C6'` | `'C24','C12','C6','C3'` |
| `C24_MG.nml` | `mesh_maps` | `'C48:C24','C24:C12','C12:C6'` | `'C24:C12','C12:C6','C6:C3'` |
| `C24_MG.nml` | `edge_cells` | `48,24,12,6` | `24,12,6,3` |
| `configuration.nml` | `&base_mesh file_prefix` | `'mesh_C48_MG'` | `'mesh_C24_MG'` |
| `configuration.nml` | `&base_mesh prime_mesh_name` | `'C48'` | `'C24'` |
| `configuration.nml` | `&multigrid chain_mesh_tags` | `'C48','C24','C12','C6'` | `'C24','C12','C6','C3'` |

The multigrid chain keeps its fixed depth of four and halves down to C3; the
`&multigrid multigrid_chain_nitems=4` line is therefore unchanged. The file
whose stem the mesh namelist is named for is the mesh it writes, which is why
this one is `C24_MG.nml`; `psy-ir-aidev/bin/generate-mesh` refuses the
mismatch by name.

## What is deliberately the same — including the length and the cadence

`&timestepping dt=3600`, `&time timestep_end='10'` and
`&io diagnostic_frequency=10` are unchanged from `example_C48`, and that is
the campaign's decision rather than an oversight. **A twelve- or fifteen-day
run is a run-time override, not a committed example.** The plan's Task C1
says so and `psy-ir-aidev/bin/run-gungho-model` implements it:

    LFRIC_EXAMPLE=example_C24 LFRIC_TIMESTEP_END=288 \
        LFRIC_DIAGNOSTIC_FREQUENCY=24 run-gungho-model minimum-production 6

is twelve days at `dt=3600` with daily output. Committing a second directory
per length would multiply the examples by the lengths and leave the reader no
way of telling which of them a figure came from; one directory per mesh, and
the length in the command, leaves the run's identity in the run's own record.

An untracked `example_c24mg_15day` existed in the canonical checkout before
this directory did, carrying `timestep_end='360'` and
`diagnostic_frequency=24` and no `C24_MG.nml` at all. This directory replaces
it: same mesh, the recipe that makes the mesh, and the two lengths back where
they belong.

## The `.xml` files are copies

All **ten** of them, unchanged, as `example_C48`'s are. None of them names a
mesh. **An edit to `example/`'s copies needs mirroring here and in the other
two**, and nothing enforces that — this paragraph is the whole of the
mechanism, for the reason `example_C48/README.md` gives.

## The mesh file is generated, not committed

`mesh_C24_MG.nc` is not in the repository;
`applications/gungho_model/example_C24/*.nc` is in `.gitignore` beside the
`example_C48` and `example_C144` rules. This directory's mesh is built by

    psy-ir-aidev/bin/generate-mesh \
        applications/gungho_model/example_C24/C24_MG.nml \
        applications/gungho_model/example_C24

run from `psy-ir-aidev` with the paths relative to `lfric_apps`. On
2026-09-14 that command wrote a 489 956-byte `mesh_C24_MG.nc` holding C24,
C12, C6 and C3 with the three intergrid maps, **byte-identical** to the mesh
the untracked `example_c24mg_15day` had been carrying — so the recipe here
reproduces the file the campaign's earlier hand-made directory ran on, and
the provenance of that file is no longer a matter of memory.

Nothing else the generator writes belongs here: it finalises through
`timer_mod` and `generate-mesh` deletes the `timer.txt` it leaves, because
`lfric-run` stages an example with `cp -r` and a stray report would be copied
into every run tree. If one appears here it is a bug in `generate-mesh`.

## `qrparm.orog.ugrid.nc` is not copied

As in `example_C48`: `&orography orog_init_option='none'` and
`file_def_ancil.xml` marks `orography_mean_ancil` `enabled=".FALSE."`, so
nothing opens it. A C16 orography sitting in a C24 directory would be a file
that looks like it belongs and would be wrong the moment the ancillary was
switched on.
