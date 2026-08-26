# Changes

## 2.0.0

- Added concrete, fully typed configuration objects for lattices, gauge
  fields, gauge actions, fermion actions, solvers, integrators, updates,
  measurements, schedules, and configuration I/O. Abstract field types are
  not stored in the new runtime objects.
- Added `Simulation`, `SimulationSpec`, and `SimulationSession` APIs suitable
  for scripts, Jupyter notebooks, MPI/GPU applications, and future GUIs.
  Sessions expose typed events and cooperative step, run, stop, and resume
  operations. `run!` prints rank-zero progress and results by default, with
  `verbose=false` for event-driven GUI or library use.
- Added a Param-free TOML reader and a versioned canonical TOML schema.
  Existing Wizard TOML files remain readable, while the original `Params`
  conversion remains available as a compatibility and debugging oracle.
- Made `run_wizard()` the new Param-free typed Wizard. The historical Wizard
  is available as `run_wizard_legacy()`; `run_wizardv2` remains an alias of
  `run_wizard`.
- Added gauge-only cold, hot, file, one-instanton, and embedded-instanton
  initialization through the Gaugefields v1 API.
- Added typed HMC, heatbath, file-loading, and self-learning HMC execution.
  SLHMC keeps the exact target action separate from the action used for the
  MD proposal. Gauge-only Sexton-Weingarten input is warned about and ignored;
  gauge/fermion force splitting remains supported.
- Added Wilson, Wilson-clover, staggered, SU(N) HISQ, standard domain-wall,
  and Möbius domain-wall fermion configurations, including stout-smeared
  actions, through the LatticeDiracOperators MD-action interface.
- Added typed measurement plans, multiple measurements, and gradient-flow
  measurements. Fermionic observables remain implemented by
  QCDMeasurements.
- Added portable JLD2 configuration checkpoints. Rank zero assembles and
  writes one global gauge configuration, and pseudofermion workspaces are not
  serialized. Bridge text and ILDG remain supported; ILDG tests are skipped
  on Windows because the current c-lime binary is unavailable there.
- Made MPI an optional weak dependency loaded through `LatticeQCDMPIExt`.
  Serial and one-GPU notebook runs can explicitly use
  `Gaugefields.SerialCommunicator()` without initializing MPI.
- Updated the supported package stack to Gaugefields 1.1.1,
  LatticeDiracOperators 1.1, LatticeMatrices 1.2.1 (through Gaugefields), and
  QCDMeasurements 1.
- Added generated regression matrices for all 464 current Wizard inputs and
  independent Wilson-clover, HISQ, and Möbius comparisons. Serial CPU, MPI
  one/two-rank, and NVIDIA H100 CUDA paths are covered.

Large, streaming JLD2 writes that avoid assembling a global configuration on
rank zero are intentionally deferred; v2.0 retains the portable rank-zero
checkpoint design.
