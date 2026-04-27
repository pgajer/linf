# Action Plan: Threshold-Free Halfvarson Calprotectin Refactor (2026-04-25)

## Goal
Replace the current prototype script [build_calprotectin_context_models.R](/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_calprotectin_context_models.R) with a manuscript-ready workflow that:

- uses the real Halfvarson validation matrix and the same sample universe as the current manuscript,
- fits a threshold-free mixed-effects dCST analysis as the primary calprotectin refinement,
- keeps the current `250 µg/g` calprotectin split as a supplemental sensitivity layer,
- and leaves the mechanistic context-aware *Faecalibacterium prausnitzii* modeling as an optional second-stage extension rather than part of the primary manuscript pipeline.

## Why this refactor is needed
The current prototype is not yet safe to use for manuscript outputs because it has four structural problems:

1. It reads from the PRIME AGP species table rather than the Halfvarson validation matrix used elsewhere in the paper.
2. It points to a nonexistent generic Halfvarson assignments file instead of the comparison-specific validation outputs.
3. It depends on undeclared gz-reading behavior (`R.utils`) and fails before building the modeling frame in the current environment.
4. It mixes two different scientific questions:
   - threshold-free dCST association with continuous calprotectin,
   - mechanistic context dependence of *F. prausnitzii*.

The manuscript needs the first of those now. The second is interesting, but it should be staged separately so the primary paper stays coherent and verifiable.

## Recommended replacement architecture

### Stage 1: Manuscript-ready threshold-free dCST analysis
Build a new source-of-truth script focused on the current manuscript question:

- **Primary question**: which rebuilt or AGP-transferred dCST labels are associated with higher continuous calprotectin burden in Halfvarson Crohn disease when repeated samples per subject are handled correctly?
- **Primary model family**: one mixed-effects logistic model per dCST label
  - `glmer(I(label == L) ~ scale(log10(calprotectin)) + (1 | subject_id), family = binomial)`
- **Primary outputs**:
  - odds ratio per +1 SD or per +1 log10 unit calprotectin increase,
  - 95% confidence interval,
  - Wald or likelihood-based p-value,
  - BH-adjusted q-value within each depth/mode comparison panel.

### Stage 2: Keep current thresholded analysis as a sensitivity layer
Retain the current `250 µg/g` split, but demote it to a supporting role:

- **Question**: do the same dCST labels separate higher inflammatory burden visits when calprotectin is dichotomized at `250 µg/g`?
- **Role in manuscript**: supplemental clinical sensitivity, not primary evidence.
- **Outputs**:
  - the current Fisher-based or count-based OR table can remain,
  - optionally add a threshold-sweep appendix later if useful.

### Stage 3: Optional mechanistic context model
After Stage 1 is stable, optionally build a separate exploratory follow-up for:

- continuous calprotectin vs CLR(*F. prausnitzii*),
- with context dependence modeled either by depth-2 dCST or by biologically motivated community scores.

This should be implemented as a separate script or a clearly isolated second mode, not intertwined with the primary manuscript generator.

## Source-of-truth inputs for the replacement script
The new script should use only the validated Halfvarson assets already aligned with the paper:

- Counts matrix:
  - [validation_counts.tsv.gz](/Users/pgajer/current_projects/gut_microbiome/outputs/validation_ingest/external_harmonized/halfvarson_2017/final_matrix/validation_counts.tsv.gz)
- Validation metadata:
  - [validation_metadata.tsv](/Users/pgajer/current_projects/gut_microbiome/outputs/validation_ingest/external_harmonized/halfvarson_2017/final_matrix/validation_metadata.tsv)
  - [validation_metadata_table.tsv](/Users/pgajer/current_projects/gut_microbiome/outputs/validation_ingest/external_harmonized/halfvarson_2017/validation_metadata_table.tsv)
- Clinical metadata with calprotectin and disease descriptors:
  - [halfvarson_supplementary_dataset1_parsed.tsv](/Users/pgajer/current_projects/gut_microbiome/data/validation_cohorts/external/halfvarson_2017/metadata/halfvarson_supplementary_dataset1_parsed.tsv)
- Current rebuilt dCST assignments:
  - [halfvarson_2017__Crohn_vs_Healthy_dcst_assignments.csv](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/absorb/halfvarson_2017/halfvarson_2017__Crohn_vs_Healthy_dcst_assignments.csv)
- If AGP-transfer calprotectin analyses are also needed, the script should explicitly identify and use the corresponding AGP-transfer assignment source already used by the manuscript rather than inferring labels ad hoc.

## Proposed script behavior

### A. Data assembly
1. Read the Halfvarson final validation matrix directly from the harmonized validation tree.
2. Read the current Crohn-comparison dCST assignment file directly from the validation outputs.
3. Join the counts, assignments, and clinical metadata by `run_accession`.
4. Restrict to Crohn samples with finite positive calprotectin.
5. Preserve the exact sample universe used for the manuscript wherever possible.
6. Emit a build summary stating:
   - total Crohn runs with calprotectin,
   - number of unique subjects,
   - number of repeated-sample subjects,
   - counts per depth/mode panel,
   - and any dropped samples with explicit reasons.

### B. Primary threshold-free dCST analysis
For each mode (`rebuilt`, `AGP-transfer`) and each depth performed in the paper:

1. Enumerate eligible labels with sufficient case count.
2. Fit one mixed-effects logistic model per label:
   - outcome: membership in the label,
   - predictor: `scale(log10(calprotectin))` or a clearly reported unscaled `log10(calprotectin)`,
   - random intercept: `subject_id`.
3. Report for each label:
   - label,
   - depth,
   - mode,
   - `n_label`,
   - `n_nonlabel`,
   - effect size as OR,
   - CI,
   - p-value,
   - q-value.
4. Apply BH correction within each depth-by-mode panel.

### C. Subject-aware sensitivity layer
Because Halfvarson is longitudinal, the new script should also include one subject-aware confirmation layer for the continuous model:

1. One-sample-per-subject resampling:
   - sample one visit per subject,
   - refit the threshold-free per-label models,
   - repeat `B` times (for example `B = 500` or `1000`),
   - report recurrence of `q < 0.05`, median OR, and median q.
2. Keep the interpretation scale consistent with the current bootstrap language already adopted in the manuscript.

This is not strictly required for the first pass, but it is strongly recommended if the threshold-free results are going to be discussed in the main text.

### D. Thresholded `250 µg/g` sensitivity layer
The new script should also regenerate the current dichotomized analysis on the same sample universe so it stays synchronized with the threshold-free results:

1. Define `high_calprotectin = calprotectin >= 250`.
2. Rebuild the current label-level association tables.
3. Emit a summary that makes clear this is:
   - a pragmatic high-inflammatory-burden contrast,
   - not the primary analysis,
   - and not a remission/activity classifier.

### E. Optional mechanistic second stage
Do not block manuscript integration on this stage.

If added later, it should live as:

- either a separate script such as `build_calprotectin_context_mechanistic_models.R`,
- or a clearly optional `--mode mechanistic` branch in the refactored script.

Preferred contents of the optional stage:

1. One continuous-outcome mixed model:
   - `log10(calprotectin) ~ Fp_abundance + (1 | subject_id)`
2. Optional context interaction model:
   - `log10(calprotectin) ~ Fp_abundance * context + (1 | subject_id)`
3. Any CLR-based compositional modeling should be explicitly labeled exploratory and should not replace the dCST-centered primary analysis.

## Deliverables

### Script deliverable
- Replace or substantially rewrite:
  - [build_calprotectin_context_models.R](/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_calprotectin_context_models.R)

### Table deliverables
Suggested outputs:

- threshold-free primary dCST table:
  - `TABLE_S13F_halfvarson_calprotectin_threshold_free.tsv`
  - `TABLE_S13F_halfvarson_calprotectin_threshold_free.tex`
- threshold-free subject-aware sensitivity:
  - `TABLE_S13G_halfvarson_calprotectin_threshold_free_subject_sensitivity.tsv`
  - `TABLE_S13G_halfvarson_calprotectin_threshold_free_subject_sensitivity.tex`
- retained thresholded `250 µg/g` sensitivity:
  - `TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tsv`
  - `TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tex`

Exact numbering can be adjusted to fit the current supplement, but the threshold-free primary and the thresholded sensitivity should be clearly separated.

### Figure deliverables
Only add figures if they materially help interpretation. The priority order is:

1. forest or dot plot of threshold-free ORs by depth/mode,
2. optional subject-aware recurrence plot,
3. mechanistic context figures only if Stage 2 is pursued.

## Manuscript integration plan

### Main text
1. Update the Halfvarson within-disease paragraph to say that calprotectin follow-up was revisited with a threshold-free mixed-effects analysis.
2. Keep the `250 µg/g` result in the text only if it adds a simple clinical gloss; otherwise move it fully to supplement mention.
3. Make the threshold-free result the preferred interpretation if it is stable and readable.

### Methods
1. Add a short subsection explaining:
   - continuous `log10(calprotectin)` modeling,
   - random intercept for repeated samples,
   - BH correction within each label panel,
   - and the role of the `250 µg/g` thresholded analysis as a sensitivity layer.

### Supplement
1. Insert the new threshold-free tables near the current Halfvarson within-disease follow-up block.
2. Keep the mechanistic context model out of the main supplemental flow unless it is mature enough to interpret responsibly.

## Acceptance criteria
The refactor is ready for manuscript use only if all of the following are true:

1. The script runs end-to-end in the current environment without relying on undeclared packages.
2. The modeling frame is built from the same Halfvarson validation branch as the manuscript.
3. The script explicitly uses the comparison-specific dCST assignments that already exist locally.
4. Threshold-free results and `250 µg/g` sensitivity results are generated from the same joined sample universe.
5. The outputs are easy to map into the supplement without ambiguous term naming.
6. The script emits enough build diagnostics that sample drift can be audited later.

## Recommended execution order
1. Replace the data-ingestion layer so the script uses the real Halfvarson validation matrix and assignment files.
2. Implement the primary threshold-free mixed-effects dCST analysis.
3. Rebuild the existing `250 µg/g` table on the same sample universe as a sensitivity layer.
4. Add subject-aware resampling for the threshold-free model if the first-pass results look promising.
5. Integrate tables into the supplement and update the main text.
6. Only then decide whether the mechanistic *F. prausnitzii* context model adds enough value to justify a second-stage analysis.

## Recommendation
For this manuscript, the safest and strongest path is:

- **Primary**: threshold-free mixed-effects dCST analysis on continuous calprotectin,
- **Secondary**: current `250 µg/g` threshold as a clinical sensitivity,
- **Optional future extension**: context-aware *F. prausnitzii* mechanism model.

That gives the paper a cleaner statistical story without giving up the clinically intuitive thresholded result, and it avoids overcommitting the manuscript to a prototype mechanistic model before the core validation pipeline is solid.
