# Relative Abundance Matrix Inventory

Date: 2026-04-29

This is a first-pass inventory for building combined AGP plus validation
relative-abundance tables for graph modeling. The source files below are
count-like taxon matrices unless noted otherwise; the relative-abundance table
for phase 1 should be created by row-normalizing each matrix after dropping any
zero-library rows.

## Candidate Matrices

| table_id | cohort | project/accession | source group | body/site | amplicon / read-family evidence | source matrix rows | usable relative rows | taxon features | taxa with max rel >= 1% | taxa with rel >= 1% in >= 1% samples | notes |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---|
| agp_local_qza_silva | American Gut Project | PRJEB11419 / Qiita 10317 | AGP local QZA | stool/gut | V4, 100 nt Deblur context | 24566 | 24566 | 3620 | 1181 | 168 | SILVA-assigned from local AGP md5 ASVs and collapsed to taxon labels. |
| hmp2_prime | HMP2 / IBDMDB | PRJNA398089 | PRIME processed | rectum/ileum/colon/cecum; not stool | V4, paired 2x250 in PRIME metadata | 176 | 176 | 4593 | 187 | 115 | Clinically strong IBD cohort, mucosal and longitudinal. |
| prjeb84421_prime | OFGCD-FI-2025 | PRJEB84421 | PRIME processed | stool/feces | V3-V4, single 300 in PRIME metadata | 73 | 73 | 4593 | 92 | 92 | Clean stool inflammatory complement, smaller and pediatric. |
| halfvarson_2017 | Halfvarson 2017 | PRJEB18471 / ERP020401 / Qiita 1629 | external reprocessed | stool/fecal | region not locally confirmed; single 100 bp family used | 591 | 591 | 806 | 249 | 158 | Longitudinal fecal IBD cohort. |
| gevers_2014 | Gevers 2014 | PRJEB13679 + PRJEB13680 | external reprocessed | final branch is stool-filtered | region not locally confirmed; single 175 bp families used | 396 | 396 | 1019 | 375 | 209 | Pediatric Crohn/control validation branch; source cohort is mixed-site. |
| jacobs_2023_ibs_250bp | Jacobs 2023 IBS | PRJNA812699 | external reprocessed | stool/fecal | region not locally confirmed; paired 250 bp family used | 239 | 239 | 710 | 214 | 155 | Baseline IBS/control branch. |
| jacobs_2023_ibs_150bp_sample | Jacobs 2023 IBS | PRJNA812699 | external reprocessed | stool/fecal | region not locally confirmed; paired 150 bp family used | 503 | 502 | 1465 | 275 | 144 | One zero-library row, `SRR18588904`, must be dropped before row normalization. |

Jacobs 150 bp will use the sample-level matrix only; the subject-level
sensitivity matrix is intentionally excluded from the combined graph-model
source set.

## Interpretation

The 1% diagnostic columns are computed after row-normalization of each source
count table:

- `taxa with max rel >= 1%`: a permissive feature screen; the feature reaches
  at least 1% relative abundance in at least one usable sample.
- `taxa with rel >= 1% in >= 1% samples`: a stricter abundance-plus-prevalence
  screen; the feature reaches 1% relative abundance in at least 1% of usable
  samples.

For a first combined graph-model pass, the stricter column is a good starting
point because it reduces extreme one-sample features while keeping the rule
simple and reproducible.

## Generated Stats

Machine-readable diagnostics were written to:

`notes/project/relative_abundance_matrix_inventory_stats_2026-04-29.tsv`

## Implemented Phase-1 Bundles

The phase-1 build script lives in:

- `/Users/pgajer/current_projects/gut_microbiome/scripts/build_phase1_relative_abundance_graph_inputs.R`

The generated output root is:

`/Users/pgajer/current_projects/gut_microbiome/outputs/graph_models/phase1_relative_abundance_2026-04-29`

The implemented phase-1 bundle is the terminal compact-label table:

`/Users/pgajer/current_projects/gut_microbiome/outputs/graph_models/phase1_relative_abundance_2026-04-29/combined/combined_relative_abundance_union_1pct.tsv.gz`

This is the recommended first `iknn_3x3` input. It combines only:

- `agp_local_qza_silva`
- `halfvarson_2017`
- `gevers_2014`
- `jacobs_2023_ibs_250bp`
- `jacobs_2023_ibs_150bp_sample`

The PRIME-derived validation tables, `hmp2_prime` and `prjeb84421_prime`, are
excluded from the primary phase-1 table because their retained features are
mostly genus-or-higher labels and are not commensurate with the AGP/external
species-like terminal labels. The primary table has 26,294 samples and 313
union features after the 1% within-dataset filter.

The table preserves terminal assigned taxon labels and does not collapse to
genus-or-higher labels. The local AGP matrix is the local-QZA SILVA assignment
built from AGP md5 ASVs using the same taxonomy-assignment path used for the
validation datasets; PRIME AGP taxonomy is not used.

## Feature Overlap Report

Before running `iknn_3x3`, feature overlap by dataset, source group, and
amplicon/read-family group was summarized in:

`/Users/pgajer/current_projects/gut_microbiome/outputs/graph_models/phase1_relative_abundance_2026-04-29/feature_overlap_report/feature_overlap_report.html`

The report uses exact terminal compact taxon labels from the phase-1 1%
relative-abundance feature set. It should therefore be interpreted as an
overlap check for the current graph-model feature representation, not as a
general biological estimate of shared taxa across differently ranked taxonomy
tables.

The report was regenerated after removing the PRIME-derived validation tables
from the primary bundle. Under exact terminal labels, AGP and the external
reprocessed stool cohorts share 151 retained features.

## PCA Dimension-Selection Report

The first implemented `iknn_3x3` preflight branch uses the primary relative
abundance table, the 1% feature filter, and centered PCA without feature
standardization (`center = TRUE`, `scale = FALSE`). The report is:

`/Users/pgajer/current_projects/gut_microbiome/outputs/graph_models/phase1_relative_abundance_2026-04-29/pca_dimension_selection/pca_dimension_selection_report.html`

The PCA thresholds are:

- 90% total explained variance: 53 PCs
- 95% total explained variance: 80 PCs
- 99% total explained variance: 140 PCs

The first 100 PCs explain 97.02% of total variance. Reusable PC scores and
loadings through PC140 were written under the same `pca_dimension_selection`
output directory so that all three candidate thresholds are covered.
