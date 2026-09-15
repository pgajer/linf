# Measured resource use

These fixed constructed workloads compare linf before and after the September
2026 transfer/landmark allocation changes. They are scale examples, not a
biological benchmark or a maximum supported dataset size.

Measurements: Apple M4 Max, 64 GiB RAM, macOS arm64, R-devel 2026-06-24
r90190, three fresh processes per case/operation. Input generation uses seed
240915, all cases fit three levels with support 2 and absorb policy. The sparse
case uses a compressed-column Matrix object; the other cases use dense R matrices.
All 36 complete result objects are identical before and after the changes.

| Case | Samples | Features | Nonzero fraction | Depth-3 groups |
|---|---:|---:|---:|---:|
| small | 1000 | 30 | 0.990 | 273 |
| sparse | 5000 | 200 | 0.035 | 323 |
| branching | 3000 | 300 | 0.033 | 301 |

## Median elapsed seconds and whole-process memory

The operation timer excludes input generation and the initial reference fit.
Peak resident memory (RSS) is measured by `/usr/bin/time` across the **whole
fresh R process**, including package loading, data generation and reference
fitting. It is not the operation's incremental memory. Small RSS differences
should not be interpreted as demonstrated memory savings.

| Case | Operation | Before (s) | After (s) | Before RSS (MiB) | After RSS (MiB) |
|---|---|---:|---:|---:|---:|
| small | fit | 0.202 | 0.185 | 300.6 | 297.7 |
| small | transfer | 0.265 | 0.163 | 305.6 | 293.3 |
| small | landmarks | 0.442 | 0.285 | 352.7 | 294.4 |
| small | embedding | 0.002 | 0.001 | 297.9 | 292.4 |
| sparse | fit | 0.844 | 0.794 | 405.0 | 406.6 |
| sparse | transfer | 3.607 | 1.242 | 380.5 | 400.8 |
| sparse | landmarks | 0.488 | 0.363 | 437.0 | 395.2 |
| sparse | embedding | 0.065 | 0.064 | 402.4 | 411.9 |
| branching | fit | 0.676 | 0.671 | 391.6 | 392.5 |
| branching | transfer | 2.921 | 0.484 | 376.2 | 366.7 |
| branching | landmarks | 0.435 | 0.339 | 402.8 | 365.4 |
| branching | embedding | 0.019 | 0.019 | 374.8 | 377.0 |

Transfer profiling identified per-candidate work and indexing as the dominant
cost. Candidate feature positions are now prepared once and scored as vectors;
parent/child counts include observed pairs only. Landmarks split membership
once and assemble result rows once. Fit and embedding algorithms are unchanged
in this comparison; variation in those timings reflects measurement noise.

The embedding deliberately returns a dense n by (p - 1) numeric matrix: the
output alone needs approximately `8 * n * (p - 1)` bytes, excluding dimnames,
attributes and temporary matrices. Sparse input does not remove that cost.
Use the benchmark runner on your intended dimensions before budgeting a large
embedding; preserve sparse matrices for supported hierarchy operations.

The reproducible runner is `tools/benchmarks/run.py`, with workload definitions
in `tools/benchmarks/measure.R`; see the repository's `CONTRIBUTING.md`. Timings
are specific to these sizes, sparsities, group structures and this machine.
There is no general linear-scaling or cross-platform speed guarantee.
