# Migrating from v1 to v2

## What remains compatible

- Existing Wizard-layout TOML files remain readable.
- `run_LQCD("parameters.toml")` and the `Params` conversion remain available
  as the legacy execution and debugging oracle.
- `run_wizard_legacy()` returns `Params` and preserves historical side effects.
- ILDG and Bridge text configuration I/O remain available; portable JLD2 is
  the default for the typed workflow and the required restart container.

## Preferred v2 workflow

```julia
using LatticeQCD

spec = load_simulation_spec("parameters.toml") # no Params construction
session = build_simulation(spec, GaugefieldsEnvironment())
summary = run!(session)                         # rank-zero progress
```

Use `run!(session; verbose=false)` for a GUI or library and attach an event
sink for structured progress. `show_config(spec)` and `show_config(session)`
provide compact human-readable summaries; ordinary REPL display does not dump
the allocated lattice arrays.

`run_wizard()` is now the Param-free Wizard and returns a typed
`SimulationSpec`. `run_wizardv2` is a compatibility alias. To normalize an old
file into the versioned schema:

```julia
spec = load_simulation_spec("old-wizard.toml")
write_simulation_spec("simulation-spec.toml", spec)
```

## Configuration and restart files

`save_configuration` writes gauge fields only. Pseudofermions, momenta, and
solver workspaces belong to one trajectory and are regenerated. A restart
checkpoint additionally stores trajectory/session counters, Metropolis RNG,
input fingerprint, and dependency versions. Save and restore only at completed
trajectory boundaries.

```julia
checkpoints = JLD2CheckpointOutput("restart"; every=10)
spec = SimulationSpec(config, schedule, OutputConfig(; checkpoints))
session = build_simulation(spec, environment)
run!(session)

continued = build_simulation(spec, environment)
load_checkpoint!(continued, "restart/restart_00000010.jld2")
run!(continued)
```

## MPI and GPU lifecycle

MPI is no longer loaded by importing LatticeQCD. Load and initialize it only
for a distributed run. Likewise, select and initialize the JACC accelerator
backend before loading LatticeQCD; a serial communicator is valid for one GPU
and Jupyter. See the [backend matrix](backends.md).

## Removed dead exports

The names `analyze`, `get_plaquette`, `get_plaquette_average`, `get_polyakov`,
`get_polyakov_average`, `get_trjs`, and `run_LQCD!` were exported by v1 even
though their implementation modules were not loaded. v2 removes those broken
exports. Use typed measurement plans/QCDMeasurements and `run_LQCD` or
`run!` instead.
