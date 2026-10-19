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

## Wizard corrections and text output in v2.0.2

The v2 wizard now treats thermalization and production inputs as separate counts.
For 50 warmup updates followed by 300 production updates starting at trajectory 1,
it writes legacy `Nthermalization = 51` and `Nsteps = 350`. Existing legacy files
retain their threshold interpretation; they are not silently reinterpreted.

Wilson pion kappa and Staggered pion mass/flavor settings now survive the wizard,
TOML, and typed conversion. Files produced by the old buggy wizard do not contain
the lost values: re-enter them or edit the nested `fermion_parameters` table.
Quenched inputs with fermion measurements now receive the required halo as well.
Both the legacy `Params` reader and the typed reader reject contradictory flat
and nested fermion values and unknown nested fermion keys. Explicit flat values
are preserved when the corresponding nested value is omitted.

Custom measurement paths also apply to flow-only runs and survive temporarily
deselecting all measurements. The wizard's verbosity setting controls the
dynamical fermion solver; session progress is controlled separately by
`run!(session; verbose=false)`, and backend diagnostics by
`GaugefieldsEnvironment(verbose=...)`.

The typed runner now honors wizard/legacy measurement and log destinations.
Files are appended on rank zero, including when `verbose=false` or when resuming
with `run!`. Parsing and building remain free of these file-creation side effects.
For explicit specifications, configure destinations with:

```julia
text = TextOutputConfig(
    logfile="logs/run.log",
    measurement_directory="measurements/run",
)
output = OutputConfig(; text)
```

A manually constructed `OutputConfig()` disables text output. To disable text
output after loading an existing file while retaining configuration/checkpoint
output, construct `OutputConfig(spec.output.configurations, spec.output.checkpoints)`
and use it in a new `SimulationSpec(spec.config, spec.schedule, output)`.

Measurement files use observable names such as `plaquette.txt` and
`pion_correlator.txt`; flowed measurements go under `flow/`. Each line contains
three tab-separated fields: trajectory, flow time (`nothing` for direct
measurements), and the complete Julia representation of the result. This format
preserves vectors and structured values and differs from the legacy file runner's
measurement-specific formats. Log files contain the same progress messages as
the console. Event sinks continue receiving the original structured values.

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
