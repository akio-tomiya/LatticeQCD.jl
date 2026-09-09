# CPU update performance: LatticeMatrices and legacy backends

This document records the post-fix CPU benchmark run on 2026-08-25. It
compares the current typed `Simulation` path using
`Gaugefields.LatticeMatricesBackend()` with the historical LatticeQCD updater
using `Gaugefields.LegacyBackend()`. Compilation, configuration construction,
measurement, and file I/O are outside the timed region.

## Result

The LatticeMatrices backend is no longer slower than the legacy backend in any
of the six common gauge-update cases. Its median time is 4--38% lower. The
largest gains are in heatbath plus overrelaxation and gauge-only HMC.

| Update | Legacy backend | LatticeMatrices backend | LM time / legacy | Change |
|---|---:|---:|---:|---:|
| SU(2), `8^4` checkerboard heatbath | 56.336 ms | 52.824 ms | 0.938 | 6.2% faster |
| SU(3), `8^4` checkerboard heatbath | 101.200 ms | 96.833 ms | 0.957 | 4.3% faster |
| SU(3), `8^4` general-action heatbath | 94.754 ms | 88.574 ms | 0.935 | 6.5% faster |
| SU(3), `8^4` checkerboard heatbath + 3 OR | 531.040 ms | 365.717 ms | 0.689 | 31.1% faster |
| SU(2), `6^4` gauge HMC, 10 MD steps | 127.784 ms | 91.057 ms | 0.713 | 28.7% faster |
| SU(3), `6^4` gauge HMC, 10 MD steps | 326.675 ms | 203.562 ms | 0.623 | 37.7% faster |

These results supersede the earlier measurements that showed a 2.4--13.2x
LatticeMatrices slowdown. That regression was real, but it has been removed by
the changes below.

## Allocation per complete update

The values below are post-warm-up `@allocated` byte counts. Allocation is no
longer proportional to a boxed argument tuple at every lattice site. The
general-action and overrelaxation cases now allocate much less than the legacy
backend. The simple heatbath and HMC paths are faster but still allocate more
than legacy and remain candidates for later workspace-level optimization.

| Update | Legacy backend | LatticeMatrices backend |
|---|---:|---:|
| SU(2) checkerboard heatbath | 55,680 | 777,832 |
| SU(3) checkerboard heatbath | 582,480 | 779,400 |
| SU(3) general-action heatbath | 32,510,592 | 936,608 |
| SU(3) checkerboard heatbath + 3 OR | 319,261,968 | 3,111,264 |
| SU(2) gauge HMC | 293,112 | 6,644,592 |
| SU(3) gauge HMC | 500,376 | 7,955,248 |

## Cause and fix

There were two independent costs.

1. Mutating LatticeMatrices kernels passed a heterogeneous vararg tuple through
   the JACC Threads loop. In production-sized kernels, Julia stopped scalar-
   replacing parts of that boundary and allocated indices and argument tuples
   per site. LatticeMatrices v1.2.1 now binds the kernel and its arguments into
   one concrete, Adapt-compatible callable before `JACC.parallel_for`. An
   `8^4` SU(3) lattice `mul!` consequently falls from 254,544 bytes to zero
   heap allocation after warm-up on one thread.
2. The generic even/odd Wilson-line evaluator copied each accumulated product
   into another whole lattice and eagerly refreshed its halo after every link.
   Gaugefields v1.1.1 specializes the four-dimensional LatticeMatrices path to
   keep products in two reusable fields and leaves halo synchronization lazy
   until a shifted read requires it.

Both changes use one backend-neutral implementation. There is no conditional
Threads/CUDA/AMDGPU/oneAPI performance path.

## Numerical and accelerator validation

- Complete seeded SU(2) and SU(3) heatbath sweeps were compared element by
  element with the legacy storage implementation using identical site-based
  random streams. The SU(2) test passed 1,290 checks; the SU(3) test passed
  2,826 checks, plus 257 Float32 checks.
- The optimized Gaugefields heatbath path was run on an NVIDIA H100 NVL. After
  13 seeded SU(3) sweeps, both CPU and CUDA produced the normalized plaquette
  `0.5589393945177883`.
- LatticeMatrices HISQ CPU tests passed 193 checks, and the same 193 checks
  passed on the H100 after removing device-side dynamic matrix allocation from
  the pullback path.

## Conditions

| Item | Value |
|---|---|
| CPU | Intel Xeon Gold 6526Y |
| Julia | 1.11.8 |
| Julia threads | 1 |
| BLAS threads | 1 |
| Initial field | cold, followed by 3 untimed updates |
| Samples | 7 |
| Minimum sample duration | approximately 0.15 s using repeated updates |
| Statistic | median seconds per complete update |
| LatticeQCD | 2.0.0 working tree |
| Gaugefields | 1.1.1 working tree |
| LatticeMatrices | 1.2.1 working tree |
| LatticeDiracOperators | 1.1.0 working tree |

The runs shared a host with unrelated jobs. Legacy and LatticeMatrices cases
were measured in the same execution context, and repeated interleaved SU(3)
heatbath measurements gave the same ordering. These medians are suitable for
detecting the original large regression; an exclusively reserved node should
still be used for publication-quality absolute numbers.

The old global random stream and the new site-based stream do not generate the
same trajectory across the two high-level fixtures. All final plaquettes were
finite and physically normal. Exact equality is checked separately by the
shared-random-stream link-level tests described above.

## Reproduction

The driver is `benchmark/cpu_common_updates.jl`. Select `legacy`,
`typed_legacy`, or `typed_lm` with `LQCD_BENCH_MODE`. For example:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
LQCD_BENCH_MODE=typed_lm LQCD_BENCH_SAMPLES=7 \
LQCD_BENCH_WARMUP=3 LQCD_BENCH_SAMPLE_SECONDS=0.15 \
julia --startup-file=no --project=. benchmark/cpu_common_updates.jl
```

Run `legacy` and `typed_lm` in the same host/cgroup context. The
`typed_legacy` mode isolates the typed LatticeQCD orchestration overhead for
heatbath; the legacy momentum implementation cannot honor the explicit seeded
momentum contract required by typed HMC and is therefore not used for those two
HMC rows.
