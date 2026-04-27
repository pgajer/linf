# n0 Schedule Comparison: `50,25,12,10` vs `50,25,25,25`

Generated: 2026-04-26

## Goal

Evaluate whether a stricter paper-facing AGP absorb hierarchy should replace the
current steadily relaxing local-QZA SILVA hierarchy.

The motivating concern is that the current schedule

- `50,25,12,10` for the depth-4 local-QZA AGP run, and
- `50,25,12,10,8,5` in the older 6-depth run

allows deep retained lineage sets with very small sample support. If `n0` is
intended as a minimum support threshold for statistically viable lineage sets,
then values such as `12`, `10`, `8`, and especially `5` are hard to defend as
the primary paper-facing default.

## Runs

### Baseline local-QZA SILVA AGP hierarchy

- AGP run:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4`
- `n0` schedule: `50,25,12,10`
- Frozen transfer root:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26`

### Stricter rerun

- AGP run:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4-n0_50_25_25_25`
- `n0` schedule: `50,25,25,25`
- Frozen transfer root:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26_n0_50_25_25_25`

## Hierarchy-Level Effect

Depths 1 and 2 are unchanged, because both schedules use the same thresholds
there. The effect is entirely in depths 3 and 4.

| depth | old n0 | new n0 | old labels | new labels | old median size | new median size | old min | new min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 50 | 50 | 50 | 50 | 132.5 | 132.5 | 56 | 56 |
| 2 | 25 | 25 | 158 | 158 | 93.5 | 93.5 | 26 | 26 |
| 3 | 12 | 25 | 402 | 219 | 44.0 | 92.0 | 12 | 26 |
| 4 | 10 | 25 | 464 | 221 | 44.0 | 92.0 | 12 | 26 |

Interpretation:

- `50,25,25,25` roughly halves the number of retained depth-3/4 labels.
- It also roughly doubles the typical depth-3/4 lineage-set support.
- Most importantly for methods defensibility, it raises the minimum retained
  lineage-set size from `12` to `26` at depths 3 and 4.

## Transfer-Level Effect

Aggregate mapped-sample totals across the same 15 comparison families:

| depth | old mapped total | new mapped total | delta |
|---|---:|---:|---:|
| 1 | 3346 | 3346 | 0 |
| 2 | 2533 | 2533 | 0 |
| 3 | 2249 | 2463 | +214 |
| 4 | 2131 | 2407 | +276 |

Depth-loss rows in the transfer audit:

- old: `44`
- new: `43`

Interpretation:

- The stricter hierarchy does **not** reduce root overlap.
- It does **increase** deeper nominal-depth assignment coverage, because more
  samples can legitimately terminate in stable retained lineage sets that remain
  valid at later depths.

Largest mapping gains occurred in:

- `gevers_2014` depth 3: `249 -> 313`
- `gevers_2014` depth 4: `222 -> 307`
- `jacobs_2023_ibs_150bp sample` depth 3: `407 -> 448`
- `jacobs_2023_ibs_150bp sample` depth 4: `390 -> 436`
- `jacobs_2023_ibs_150bp subject` depth 3: `300 -> 329`
- `jacobs_2023_ibs_150bp subject` depth 4: `286 -> 319`

Halfvarson changed only modestly at depth 3/4, and HMP2 remained depth-1-only
under both schedules.

## Association-Level Effect

Top-10 rows with `q <= 0.10` across all comparison-depth summaries:

- old: `48`
- new: `63`

This increase should **not** be interpreted as stronger biology across the
board. It mainly reflects coarser depth-3/4 state spaces and smaller multiple-
testing burdens.

### IBD-facing takeaways

#### Halfvarson

Depth 1 and depth 2 are unchanged and remain the strongest transfer layers.

At depths 3 and 4, the stricter hierarchy makes the pooled and Crohn signals
less extreme but still significant, while the UC depth-3/4 signal drops below
`q <= 0.10`.

Examples:

- `IBD_vs_Healthy`
  - depth 3 best `q`: `1.58e-04 -> 2.26e-02`
  - depth 4 best `q`: `2.17e-04 -> 2.54e-02`
- `Crohn_vs_Healthy`
  - depth 3 best `q`: `9.37e-03 -> 3.93e-02`
  - depth 4 best `q`: `1.32e-02 -> 4.23e-02`
- `UC_vs_Healthy`
  - depth 3 best `q`: `1.97e-03 -> 0.163`
  - depth 4 best `q`: `2.46e-03 -> 0.186`

So the stricter hierarchy is more conservative for the main IBD transfer signal
at deeper depths.

#### Gevers

Still null. The deeper mapping gain does not recover corrected Crohn-vs-healthy
transfer.

#### HMP2

Unchanged. The current mucosal-only slice remains effectively depth-1-only.

### IBS / OFGCD side effects

The stricter hierarchy does help some non-IBD exploratory comparisons:

- Jacobs 150bp sample-level IBS gains depth-3/4 support
  - depth 3 best `q`: `2.05e-02 -> 3.33e-03`
  - depth 4 best `q`: `0.133 -> 3.85e-03`
- Jacobs 150bp subject-level IBS gains depth-3/4 suggestive support
  - depth 3 best `q`: `0.387 -> 8.33e-02`
  - depth 4 best `q`: `1.00 -> 9.21e-02`
- PRJEB84421 deeper nominal depths become much more stable, but this mostly
  reflects stable lineage-set carry-through rather than the emergence of new
  strong disease-vs-control signals.

## Recommendation

For the methods paper, I would use:

- **primary hierarchy:** `50,25,25,25`
- **sensitivity / finer-grained exploratory hierarchy:** `50,25,12,10`

Why:

1. `50,25,25,25` is much easier to justify as a minimum-support rule for
   retained lineage sets.
2. It still preserves the major depth-1/2 IBD transfer story.
3. It improves deeper assignment coverage.
4. It avoids paper-facing dependence on lineage sets with support near `10` or
   `12`.

Why not use it as the only story:

1. It coarsens depth-3/4 structure substantially.
2. It weakens the strongest deeper Halfvarson signals, especially UC.
3. It changes the identity of some “best” deeper transferred labels simply by
   merging more fine-grained descendants.

So the cleanest paper-facing position is:

- `50,25,25,25` as the **default, defensible minimum-support hierarchy**
- `50,25,12,10` as a **resolution-increasing sensitivity analysis** that shows
  what is gained, and what becomes less statistically stable, when smaller deep
  lineage sets are allowed.
