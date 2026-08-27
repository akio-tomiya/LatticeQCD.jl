# Backend support and qualification

LatticeQCD v2 does not branch its algorithms by accelerator vendor. Field
allocation, random streams, Dirac operators, MD, and I/O use the common
Gaugefields/LatticeMatrices/JACC interfaces. The table below distinguishes
implemented code paths from configurations actually qualified for v2.0.

| Execution mode | v2.0 status | Release qualification |
|---|---|---|
| Serial CPU | Supported | Full package tests; all 464 Wizard inputs; exact restart matrix |
| CPU MPI, 1 and 2 ranks | Supported | CI lifecycle, gauge HMC, Wilson HMC, JLD2 configuration and exact restart |
| One NVIDIA CUDA GPU | Supported for tested cases | H100: gauge, staggered and HISQ HMC; portable JLD2; fermion restart |
| MPI + CUDA multi-GPU | Experimental | Components implement MPI and CUDA, but their combination is not yet release-qualified |
| AMDGPU | Experimental | Backend-neutral path exists; no v2.0 hardware qualification |
| oneAPI | Experimental | Backend-neutral path exists; no v2.0 hardware qualification |

“Experimental” means that v2 deliberately contains no vendor-specific CPU
fallback, but the configuration is not part of the compatibility promise
until it passes hardware regression tests.

## Optional MPI

MPI is a weak dependency. A notebook or one-GPU process can use
`Gaugefields.SerialCommunicator()` without loading or initializing MPI. A
distributed process must load MPI explicitly before constructing the
environment, and every rank must call configuration/checkpoint save and load.

## Optional CUDA qualification

The hardware test is not part of ordinary hosted CI. On a CUDA runner with
LatticeQCD, JACC, and CUDA instantiated, select the backend in one Julia
process, restart Julia, and run:

```sh
julia --project=. -e 'using JACC; JACC.set_backend("cuda")'
julia --startup-file=no --project=. test/gpu_smoke.jl
```

`LQCD_CUDA_CASES=gauge,staggered,hisq` selects smoke cases, and
`LQCD_CUDA_RESTART_CASE=staggered` selects the exact restart case. The test
requires a real CUDA device and fails instead of silently falling back to CPU.
On a multi-GPU host, set `CUDA_DEVICE_ORDER=PCI_BUS_ID` and select the device
with `CUDA_TEST_DEVICE` before starting Julia.

## Portable files versus exact replay

Ordinary JLD2 configuration files and restart checkpoints contain one global
gauge configuration assembled by rank zero. They are intended to be readable
with a different CPU/GPU backend or MPI decomposition. Exact replay additionally
requires matching physical input, numeric type, package behavior, and
deterministic backend kernels. Package version differences therefore warn by
default, while `strict_versions=true` rejects them.

Large streaming writes that avoid rank-zero assembly remain outside v2.0.
