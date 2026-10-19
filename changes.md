# Changes

## 2.0.2

- Preserve pion-measurement Wilson kappa and Staggered mass/flavor settings in
  both wizards, saved TOML, typed specifications, and gradient-flow measurements.
  Revisiting the v2 measurement form now retains previously entered values.
- Reject unknown nested fermion keys and conflicting flat/nested settings;
  explicit flat values no longer depend on dictionary iteration order.
- Include valence measurements when selecting gauge-field halo storage, including
  legacy quenched and file-loading inputs and measurements after gradient flow.
  Explicit typed specifications with insufficient measurement halos fail validation.
- Interpret v2 wizard thermalization and production inputs as counts and encode
  them correctly in legacy TOML. Existing legacy TOML thresholds remain unchanged.
- Add serializable `TextOutputConfig` for rank-zero log and measurement files.
  Wizard output paths now reach the typed runner, including quiet runs and
  stop/resume. Complete measurement values are appended separately from the
  compact console display; flowed measurements have a separate directory.
- Make the quick start run a simulation once and explain output files, optional
  plotting, and HMC step-size tuning.
- Add regression coverage using nondefault measurement parameters, direct and
  flowed measurements, canonical/legacy TOML round trips, and actual file output.
- Preserve custom measurement paths for flow-only runs and while measurements
  are temporarily deselected. Apply the same nested-parameter validation in the
  legacy `Params` reader and the typed reader.
- Reject zero starting trajectories, negative random seeds, and odd
  Sexton-Weingarten step counts at the wizard prompt. Clarify that its verbosity
  setting controls the dynamical fermion solver.
- Audit all wizard pages with nondefault values, page revisits, both TOML formats,
  all dynamical fermion/flavor choices, direct/flowed observables, input/output
  formats, and an executed flow-only run with independent save/checkpoint intervals.

## 2.0.0

### Simulation architecture

- Added concrete, fully typed configuration objects for lattices, gauge
  fields, gauge actions, fermion actions, solvers, integrators, updates,
  measurements, schedules, and configuration I/O. Abstract field types are
  not stored in the new runtime objects.
- Added `Simulation`, `SimulationSpec`, and `SimulationSession` APIs suitable
  for scripts, Jupyter notebooks, MPI/GPU applications, and future GUIs.
  Sessions expose typed events and cooperative step, run, stop, and resume
  operations. `run!` prints rank-zero progress and results by default, with
  `verbose=false` for event-driven GUI or library use.
- Added concise `text/plain` displays for typed specifications, simulations,
  sessions, schedules, output settings, and run summaries. Assigning a built
  session in a REPL or Jupyter cell therefore no longer dumps the underlying
  gauge-field storage.

### Parameter files and Wizard

- Added a Param-free TOML reader and a versioned canonical TOML schema.
  Existing Wizard TOML files remain readable, while the original `Params`
  conversion remains available as a compatibility and debugging oracle.
- Made `run_wizard()` the new Param-free typed Wizard. The historical Wizard
  is available as `run_wizard_legacy()`; `run_wizardv2` remains an alias of
  `run_wizard`. Simple mode now identifies its two-flavor Wilson HMC preset
  when selected and prints a human-readable summary of the actual fermion,
  smearing, update, measurement, and output choices before writing TOML.
  Back navigation now says that it returns to the previous section and
  discards that section's unfinished edits; text and number prompts accept
  the plain command `back` and report the destination after navigation. The
  opening v2 banner now explains both `back` and the clean `quit` action where
  the old Ctrl+C-only exit message appeared.

### Configurations, updates, and fermions

- Added gauge-only cold, hot, file, one-instanton, and embedded-instanton
  initialization through the Gaugefields v1 API.
- Use the Grid/Bridge++ momentum normalization in HMC configurations converted
  for the typed v2 workflow: Gaussian Lie-algebra coefficients have width
  `sqrt(2)`, the kinetic term is `p*p/4`, and momentum kicks carry a factor of
  two. An explicitly constructed `GaussianMomentumConfig(1.0, ...)` retains
  the historical LTK normalization, and the legacy `Params` execution path is
  unchanged. To compare the same trajectory with v1, convert the MD step size
  as `delta_tau_v2 = delta_tau_v1 / sqrt(2)`.
- Require Gaugefields v1.1.4 for its normalization-aware MD driver.
- Added typed HMC, heatbath, file-loading, and self-learning HMC execution.
  SLHMC keeps the exact target action separate from the action used for the
  MD proposal. Gauge-only Sexton-Weingarten input is warned about and ignored;
  gauge/fermion force splitting remains supported.
- Added Wilson, Wilson-clover, staggered, SU(N) HISQ, standard domain-wall,
  and Möbius domain-wall fermion configurations, including stout-smeared
  actions, through the LatticeDiracOperators MD-action interface.

### Measurements and gradient flow

- Added typed measurement plans, multiple measurements, and gradient-flow
  measurements. Fermionic observables remain implemented by
  QCDMeasurements.

### Configuration I/O and restart

- Added portable JLD2 configuration output. Rank zero assembles and
  writes one global gauge configuration, and pseudofermion workspaces are not
  serialized. Bridge text and ILDG remain supported; ILDG tests are skipped
  on Windows because the current c-lime binary is unavailable there.
- Added safe periodic HMC/SLHMC restart checkpoints. A checkpoint is written
  to a `.pending` JLD2 file and only renamed after both the global gauge
  configuration and trajectory-boundary state are complete. The interval is
  independent of ordinary configuration output, and restoring preserves the
  Metropolis RNG, counters, and session progress. The Wizard offers the
  interval in both simple and expert HMC modes. Dynamical-fermion trajectories
  clear chronological Krylov guesses at each boundary, so pseudofermions can
  be regenerated from their seed and trajectory after a portable restart.
  Checkpoints include a SHA-256 fingerprint of the typed physical input and
  numeric type plus Julia/package versions. Incompatible physics is rejected
  before loading; dependency differences warn or can be made strict.

Large, streaming JLD2 writes that avoid assembling a global configuration on
rank zero are intentionally deferred; v2.0 retains the portable rank-zero
checkpoint design.

### MPI, GPU, and package versions

- Made MPI an optional weak dependency loaded through `LatticeQCDMPIExt`.
  Serial and one-GPU notebook runs can explicitly use
  `Gaugefields.SerialCommunicator()` without initializing MPI.
- Updated the minimum supported package stack to Gaugefields 1.1.4,
  LatticeDiracOperators 1.1.2, and QCDMeasurements 1. LatticeMatrices is used
  through Gaugefields.

### Regression and release qualification

- Added generated regression matrices for all 464 current Wizard inputs and
  independent Wilson-clover, HISQ, and Möbius comparisons. Serial CPU, MPI
  one/two-rank, and NVIDIA H100 CUDA paths are covered.
- Verified the Grid/Bridge++ momentum normalization on threaded CPU and two
  MPI ranks, and on an NVIDIA H100 for gauge, staggered, and HISQ HMC plus an
  exact staggered checkpoint restart. See `comparison.md` for the numerical
  fingerprints and the old/new MD-time conversion comparison.
- Added bitwise restart regressions for Wilson, Wilson-clover, staggered RHMC,
  HISQ RHMC, standard/Möbius domain-wall, stout Wilson, and fermionic SLHMC.
- Added public-API and Aqua quality gates. Removed legacy names that were
  exported without loaded implementations, and fixed unbound type parameters
  in the lattice and legacy heatbath constructors.
- Declared the `Logging` standard library as a test dependency so isolated
  `Pkg.test` environments can run the Param-free TOML input regressions.
- Resolve legacy `Params` log and measurement directories with platform-native
  path handling, including absolute paths on a different Windows drive.
- Documented the v2 backend qualification matrix and v1 migration path.
- Verified the notebook workflow through an actual IJulia kernel on Julia
  1.11, including gauge-only heatbath, an existing Wilson-plus-stout TOML
  input, visible `run!` progress, and quiet `run!(...; verbose=false)` use.
