# Action Plan: Absorb dCST Rerun for the Gut Paper

Date: 2026-04-10

Consulted comparison report:
- [`merged_ct_dcst_manuscript__2026_04_10.tex`](/Users/pgajer/current_projects/CT_clearance/docs/merged_ct_dcst_manuscript__2026_04_10.tex)

Target manuscript:
- [`gut_application_paper.tex`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex)

## Purpose

Rerun the gut-paper DCST analyses using the absorb low-frequency policy and
test whether the hierarchy supports scientifically useful refinement beyond the
currently reported depth 2.

This is not a cosmetic relabeling exercise. The absorb dCST hierarchy must be built
from the start and then refined recursively, because deeper absorb-only splits
are not recoverable by simply switching a rare-built hierarchy to the absorb
view afterward.

## Core Decision

Use absorb as the primary hierarchy for the rerun and treat the current
rare-policy paper results as the baseline comparator.

The main idea borrowed from the CT-clearance report is:

- build one absorb dCST hierarchy recursively,
- summarize structure by depth before making phenotype claims,
- use omnibus depth-level tests as the gatekeeper for going deeper,
- then do focused label-level and biological follow-up only where the deeper
  structure is supported and interpretable.

## Recommended Implementation Path

Use the `linf` R package directly rather than extending the older Python
prototype in
[`agp_dcst_analysis.py`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/agp_dcst_analysis.py).

Why:

- the package already supports absorb recursion cleanly through
  [`linf.csts()`](/Users/pgajer/current_projects/linf/R/linf.R),
  [`refine.linf.csts()`](/Users/pgajer/current_projects/linf/R/linf.R), and
  [`refine.linf.csts.iter()`](/Users/pgajer/current_projects/linf/R/linf.R),
- the old Python script is a 5,000-sample rare-policy prototype and is no
  longer a good source of truth for the manuscript-scale analysis,
- the CT-clearance absorb workflow is conceptually much closer to the package
  implementation than to the older Python prototype.

## Run Roots And Output Convention

Create a new absorb-specific run rather than overwriting current rare-policy
artifacts.

Recommended discovery run root:
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-10-absorb-depthscan/`

Recommended validation run root:
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/absorb_2026-04-10/`

Recommended paper-facing staging note:
- [`papers/gut-dcst/notes/action_plans/ACTION_PLAN_absorb_dcst_rerun_2026-04-10.md`](/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/action_plans/ACTION_PLAN_absorb_dcst_rerun_2026-04-10.md)

Do not replace the current manuscript-facing CSVs and figures until the absorb
branch passes the stop-go checks below.

## Depth Schedule

Use a fixed absorb discovery schedule modeled on the CT-clearance report:

- depth 1: `n0 = 50`
- depth 2: `n0 = 25`
- depth 3: `n0 = 12`
- depth 4: `n0 = 10`
- depth 5: `n0 = 8`
- depth 6: `n0 = 5`

Interpretation rule:

- depths 1 to 4 are inferential candidates,
- depths 5 and 6 are structural and exploratory unless they remain well
  populated and phenotype-supported.

## Immediate Priorities

Prioritize phenotype branches in this order based on the current rare-policy
outputs:

1. IBD
2. Autoimmune
3. IBS
4. Seasonal allergies
5. Acid reflux

IBD is the primary deep branch because it is currently the strongest and most
coherent phenotype family in the paper and already carries the external
validation story.

## Work Packages

## Package 0: Freeze The Baseline

Goal:

- preserve the current rare-policy outputs as the baseline comparator.

Tasks:

- record the exact current discovery files used by the manuscript:
  - [`full_cohort_adjusted_results.csv`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_adjusted_results.csv)
  - [`depth2_adjusted_results.csv`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/depth2_adjusted_results.csv)
  - [`sensitivity_clean_adjusted_results.csv`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/sensitivity_clean_adjusted_results.csv)
  - [`full_cohort_dcst_assignments.csv`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_dcst_assignments.csv)
- copy or symlink the baseline files into the absorb run root under
  `baseline_rare_reference/`.
- save a short machine-readable manifest:
  `baseline_manifest.tsv`

Deliverables:

- `baseline_rare_reference/`
- `baseline_manifest.tsv`

## Package 1: Rebuild The Discovery Hierarchy In Absorb Mode

Goal:

- generate one recursive absorb dCST hierarchy on the AGP discovery table using the
  same sample and feature filtering as the current paper.

Tasks:

- load the same PRIME-derived SILVA species table and metadata as the current
  paper.
- apply the current discovery filters:
  - minimum library size `>= 1000`
  - feature retention at count `>= 2` in at least `5%` of retained samples
- build the absorb dCST hierarchy from depth 1 to depth 6 using the fixed schedule
  above.
- save both active labels and policy-specific views at every depth.

Recommended outputs:

- `agp_absorb_assignments.tsv.gz`
- `agp_absorb_label_summary_by_depth.tsv`
- `agp_absorb_reassignment_summary.tsv`
- `agp_absorb_object.rds`

Required columns in `agp_absorb_assignments.tsv.gz`:

- `Run`
- `depth1_absorb`
- `depth2_absorb`
- `depth3_absorb`
- `depth4_absorb`
- `depth5_absorb`
- `depth6_absorb`
- `dominant_taxon_depth1`
- `dominant_pair_depth2`

## Package 2: Structural Depth Scan

Goal:

- determine whether deeper absorb refinement remains meaningful before doing
  phenotype-heavy analyses.

Tasks:

- summarize label counts by depth.
- summarize largest-label sizes by depth.
- quantify how many samples are absorbed at each depth.
- generate depth-transition cross-tabs and heatmaps:
  - `depth1 -> depth2`
  - `depth2 -> depth3`
  - `depth3 -> depth4`
  - optionally `depth4 -> depth5`
- record whether deeper levels continue to branch nontrivially or mostly
  fragment.

Recommended outputs:

- `depth_summary.tsv`
- `largest_labels_by_depth.tsv`
- `depth12_crosstab.tsv`
- `depth23_crosstab.tsv`
- `depth34_crosstab.tsv`
- `depth12_heatmap.png`
- `depth23_heatmap.png`
- `depth34_heatmap.png`

Decision checkpoint:

- if depth 3 already collapses into too many tiny labels with weak parent-child
  structure, stop the inferential branch at depth 2 and keep deeper levels
  supplementary only.

## Package 3: Omnibus Phenotype Scan By Depth

Goal:

- use depth-level tests to decide which phenotype branches are worth taking
  deeper.

Tasks:

- for each phenotype and each depth, compute an omnibus label-by-phenotype
  association test.
- summarize effect size with Cramer's V.
- adjust across the full phenotype-by-depth panel with BH.

Recommended phenotype set:

- IBD
- IBS
- Autoimmune
- Seasonal allergies
- Acid reflux
- Kidney disease
- CDI
- Lung disease
- Migraine
- Diabetes
- Cardiovascular disease
- Obesity

Recommended outputs:

- `omnibus_by_depth.tsv`
- `omnibus_q_heatmap.png`
- `omnibus_v_heatmap.png`

Minimum go criterion for deeper follow-up:

- omnibus `q <= 0.10` at the candidate depth
- and Cramer's V does not collapse relative to the immediately shallower depth

Main-text criterion:

- omnibus support at the candidate depth plus at least one interpretable parent
  lineage with adequate sample size for follow-up

## Package 4: Label-Level Association Reruns

Goal:

- rerun the paper's label-level analyses on absorb labels for all supported
  phenotype-depth panels.

Tasks:

- rerun univariate 2x2 association tests at every supported depth.
- rerun covariate-adjusted models analogous to the current paper.
- keep the same covariate logic as the rare-policy paper unless a real defect is
  found.
- produce depth-specific association tables instead of forcing everything into
  one depth-1-centric table layout.

Recommended outputs:

- `depth1_univariate_results.tsv`
- `depth2_univariate_results.tsv`
- `depth3_univariate_results.tsv`
- `depth4_univariate_results.tsv`
- `depth1_adjusted_results.tsv`
- `depth2_adjusted_results.tsv`
- `depth3_adjusted_results.tsv`
- `depth4_adjusted_results.tsv`
- `headline_absorb_results.tsv`

Interpretation rule:

- do not treat a deeper level as useful merely because it produces more small
  labels; it must either sharpen an existing phenotype signal or reveal a
  coherent lineage that the shallower analysis hid.

## Package 5: Contamination-Aware Sensitivity Rerun

Goal:

- test whether the absorb results remain supported after the manuscript's
  contamination-aware filtering.

Tasks:

- apply the same contamination-aware filtering logic used in the current paper.
- rerun the adjusted association analyses at the supported absorb depths.
- compare retained versus dropped signals against the baseline rare-policy
  results.

Recommended outputs:

- `clean_depth1_adjusted_results.tsv`
- `clean_depth2_adjusted_results.tsv`
- `clean_depth3_adjusted_results.tsv`
- `clean_depth4_adjusted_results.tsv`
- `absorb_vs_clean_retention_summary.tsv`

Key question:

- does absorb mainly rescue interpretable IBD-like lineages, or does it mostly
  increase oral/skin-sensitive instability?

## Package 6: Parent-Child Lineage Follow-Up

Goal:

- replace the current one-off depth-2 case study with lineage analysis modeled
  on the CT report.

Tasks:

- identify the best supported absorb lineage, almost certainly starting with the
  IBD branch.
- compare the selected child dCST state against the remainder of its parent dominance-lineage.
- if depth 4 is supported, repeat one step deeper only for that selected dominance-lineage.
- keep this as a phenotype-specific lineage analysis, not a full feature
  discovery project.

Recommended outputs:

- `lineage_priority_table.tsv`
- `ibd_parent_child_summary.tsv`
- `ibd_parent_child_plot.png`
- optionally `ibd_depth4_descendants.tsv`

Decision rule:

- only take a lineage one step deeper if the child-versus-parent contrast adds a
  clearer biological story than the shallower label already provided.

## Package 7: Taxonomy Concordance Rerun

Goal:

- test whether the absorb dCST hierarchy is stable across SILVA and GG2 after
  harmonization.

Tasks:

- rebuild the absorb dCST hierarchy on the overlapping GG2 table.
- compare SILVA and GG2 absorb labels at depths 1 to 3.
- summarize exact-label agreement, harmonized agreement, and large-state
  relabeling patterns.

Recommended outputs:

- `silva_gg2_absorb_concordance_depth1.tsv`
- `silva_gg2_absorb_concordance_depth2.tsv`
- `silva_gg2_absorb_concordance_depth3.tsv`
- `silva_gg2_absorb_heatmap_depth1.png`
- `silva_gg2_absorb_heatmap_depth2.png`

Interpretation rule:

- taxonomy concordance should be used to calibrate naming and stability, not as
  a requirement that labels match exactly character-for-character.

## Package 8: External Validation Rerun

Goal:

- rerun completed validation cohorts under absorb and allow selected cohorts to
  extend to depth 3 or 4 when sample size permits.

Tasks:

- build an absorb-capable validation runner rather than reusing the current
  rare-only, depth-1/2-only validator in
  [`run_external_case_control_validation.py`](/Users/pgajer/current_projects/gut_microbiome/scripts/run_external_case_control_validation.py)
  unchanged.
- keep cohort-scaled thresholds for smaller validation cohorts, but extend the
  threshold chooser to depth 3 and depth 4.
- rerun first:
  - HMP2 / IBDMDB
  - Halfvarson 2017
  - Gevers 2014
  - PRJEB84421
- add IBS-focused external cohorts afterward, especially the in-progress Jacobs
  branch, once those matrices are ready.

Recommended outputs:

- `validation_absorb_summary.tsv`
- one absorb assignment file per cohort and comparison
- one depth-specific result table per cohort and comparison
- one compact cohort-level comparison figure to discovery

Decision rule:

- validation should test directional and structural generalization, not require
  exact recovery of every AGP label.

## Package 9: Rare Versus Absorb Comparison

Goal:

- quantify what actually changes when low-frequency samples are absorbed instead
  of sent to `RARE_DOMINANT`.

Tasks:

- compare label counts, omnibus support, strongest phenotype hits, and the fate
  of the current `RARE_DOMINANT` signals.
- track where the rare-policy tail is redistributed in the absorb dCST hierarchy.

Recommended outputs:

- `rare_vs_absorb_overview.tsv`
- `rare_vs_absorb_omnibus_comparison.tsv`
- `rare_vs_absorb_lineage_transfer.tsv`
- `rare_vs_absorb_summary_figure.png`

This comparison should likely become a supplement section unless absorb becomes
the new primary analysis for the manuscript.

## Recommended Rerun Order

1. Package 0: baseline freeze
2. Package 1: absorb discovery hierarchy
3. Package 2: structural depth scan
4. Package 3: omnibus phenotype scan
5. Stop-go decision on whether depth 3 and depth 4 deserve inferential follow-up
6. Package 4: label-level association reruns
7. Package 5: contamination-aware rerun
8. Package 6: selected dominance-lineage follow-up
9. Package 7: taxonomy concordance rerun
10. Package 8: external validation rerun
11. Package 9: rare-versus-absorb comparison
12. Only then update manuscript-facing figure/table builders

## Stop-Go Rules

Stop at depth 2 for the main paper if any of the following occurs:

- depth 3 becomes structurally fragmented without clear parent-child
  refinement,
- omnibus support weakens materially at depth 3 for the best phenotype branches,
- deeper labels become too small to support stable adjusted analyses,
- validation does not support the deeper branch even directionally.

Promote depth 3 to the main paper if most of the following occur:

- the hierarchy remains well populated at depth 3,
- at least one phenotype branch, ideally IBD, shows omnibus support at depth 3,
- the parent-child lineage interpretation is clearer than the current depth-2
  story,
- contamination-aware reruns retain the main signal,
- validation is at least directionally supportive.

Promote depth 4 only if:

- depth-4 support is reproducible at the omnibus level,
- the additional split is lineage-coherent rather than just sparse,
- the biology becomes clearer rather than noisier.

## Manuscript Integration Plan

If absorb outperforms rare:

- make absorb the primary discovery analysis in the paper,
- keep rare versus absorb as a supplement or brief robustness subsection,
- rewrite the current depth-2 refinement section into a lineage-refinement
  section,
- update the IBD validation section to match the absorb branch.

If absorb is informative but not clearly better:

- keep the current rare-policy paper as the main text,
- report absorb-policy dCSTs as a structural supplement with selected dominance-lineage follow-up.

If absorb fails to add stable signal:

- retain the current rare-policy paper architecture and close the absorb branch
  as a negative but useful calibration exercise.

## First Concrete Deliverables To Build

The first implementation sprint should produce only these files:

- `agp_absorb_object.rds`
- `agp_absorb_assignments.tsv.gz`
- `depth_summary.tsv`
- `largest_labels_by_depth.tsv`
- `depth12_heatmap.png`
- `depth23_heatmap.png`
- `depth34_heatmap.png`
- `omnibus_by_depth.tsv`
- `omnibus_q_heatmap.png`
- `omnibus_v_heatmap.png`
- `GO_NO_GO_depth3_depth4.md`

Only after those exist should the deeper phenotype-specific and validation work
begin.
