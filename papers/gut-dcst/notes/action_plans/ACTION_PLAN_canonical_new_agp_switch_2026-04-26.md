# Action Plan: Canonical New-AGP Switch

Date: 2026-04-26

## Goal

Switch the gut-dCST manuscript’s analysis layer to the taxonomy-harmonized AGP
reference and the corrected frozen-transfer algorithm, then rebuild the
manuscript-facing downstream assets from that canonical source.

Canonical roots:

- AGP hierarchy:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4-n0_50_25_25_25`
- Frozen transfer:
  `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26_n0_50_25_25_25`

## Step 1. Canonical source switch

1. Patch manuscript-facing builder scripts so they no longer read:
   - `outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive`
   - `outputs/dcst_validation/frozen_agp`
2. Replace those roots with the canonical local-QZA AGP hierarchy and corrected
   transfer root above.
3. Where practical, centralize Python-side constants in `paper_paths.py` or in
   the builder itself to reduce future drift.

## Step 2. Rebuild AGP follow-up layer on the canonical hierarchy

1. Run `run_absorb_label_followup.R` on the canonical AGP run for at least:
   - `IBD`
   - `IBS`
   - `Autoimmune`
2. Confirm the `ibd_label_followup` outputs needed by:
   - Figure 2A
   - Supplementary Figure S1
   - AGP follow-up tables
   are present under the canonical run.

## Step 3. Rebuild manuscript-facing validation assets

1. Regenerate the main external-validation asset table/figure:
   - `build_validation_assets.py`
2. Regenerate the IBD overview figure:
   - `build_ibd_results_overview_figure.py`
3. Regenerate supplement reviewer-sensitivity assets:
   - `build_ibd_reviewer_sensitivity_assets.R`
4. Regenerate the calprotectin threshold-sensitivity assets:
   - `build_calprotectin_context_models.R`
5. Regenerate the validation-dataset appendix tables:
   - `build_external_validation_association_tables.py`

## Step 4. Update AGP discovery-side retained-feature summaries

1. Recompute the retained-feature lineage profile table from the canonical AGP
   counts and assignments.
2. Update the Supplementary Table S7 values in the supplement source if they
   change.

## Step 5. Rebuild manuscript PDFs and inspect drift

1. Rebuild:
   - main manuscript PDF
   - supplement PDF
2. Check for:
   - figure drift
   - table drift
   - changed quantitative claims in main text
3. Summarize any narrative updates still needed after the asset switch.

## Expected deliverables

- canonical-source decision note
- patched builder scripts
- canonical AGP follow-up outputs
- regenerated manuscript figures/tables based on the new AGP hierarchy and
  corrected transfer algorithm
- rebuilt PDFs

## Execution Status

Completed on 2026-04-26.

### Completed

1. Patched manuscript-facing builders to use:
   - canonical AGP hierarchy:
     `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-26-agp-silva-local-qza-absorb-depth4-n0_50_25_25_25`
   - canonical frozen transfer:
     `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/frozen_agp_silva_local_qza_2026-04-26_n0_50_25_25_25`
2. Reran AGP absorb-label follow-up for:
   - `IBD`
   - `IBS`
   - `Autoimmune`
3. Regenerated manuscript-facing validation assets:
   - `TABLE_V1_external_validation_summary.tsv`
   - `FIGURE_2_ibd_results_overview.png`
   - `FIGURE_4_external_validation_modes.png`
   - `FIGURE_S1_agp_ibd_followup_states.png`
4. Regenerated reviewer-sensitivity and subject-aware tables from the canonical
   roots, including updated `S8`–`S13` assets.
5. Regenerated calprotectin threshold-sensitivity assets (`S13H`, `S13I`,
   `Figure S2`) from the canonical transfer layer.
6. Regenerated the external-validation appendix tables.
7. Rebuilt:
   - `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
   - `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper_supplement.pdf`

### Important Outcome

The switch is technically successful, but it changes the AGP discovery-side IBD
story:

- the canonical local-QZA AGP `IBD` follow-up no longer yields frequentist
  depth-1 to depth-4 `q < 0.05` labels;
- the strongest retained AGP IBD labels shift toward states such as
  `Segatella sp.`, `Faecalibacterium prausnitzii`,
  `Segatella sp.__Faecalibacterium prausnitzii`, and
  `Agathobacter faecis__Bacteroides dorei__Faecalibacterium prausnitzii`;
- transfer-side results remain regenerated and coherent under the corrected
  shared-taxonomy workflow.

### Residual Editorial Work

The main remaining task is narrative rather than computational:

- the hand-written AGP discovery discussion and any hard-coded
  `Bacteroides/Lachnospiraceae/Alistipes` prose in the manuscript and supplement
  should be revised to match the canonical local-QZA AGP follow-up outputs.
