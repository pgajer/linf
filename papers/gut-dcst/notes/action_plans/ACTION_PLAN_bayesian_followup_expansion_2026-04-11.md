# Action Plan: Expand Absorb Bayesian Follow-up Beyond IBD, IBS, and Autoimmune

Date: 2026-04-11

Primary manuscript:
- [`gut_application_paper.tex`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex)

Current absorb discovery backbone:
- [`2026-04-11-absorb-depthscan-adaptive`](/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive)

Current follow-up runner:
- [`run_absorb_label_followup.R`](/Users/pgajer/current_projects/gut_microbiome/scripts/run_absorb_label_followup.R)

## Status Update

Batch 1 was completed on 2026-04-11 using:

- `Acid_reflux`
- `Seasonal_allergies`
- `Cardiovascular_disease`
- `Obesity`

Current interpretation:

- `Acid_reflux` is the clearest new main-text candidate.
- `Seasonal_allergies` is the strongest secondary candidate and can support a
  lighter main-text treatment or a high-priority supplement section.
- `Cardiovascular_disease` is strong statistically but is better kept
  supplement-first because the lineage-level interpretation is more vulnerable
  to confounding.
- `Obesity` does not currently justify manuscript expansion.

## Purpose

Extend the absorb-policy dCST dominance-lineage frequentist/Bayesian follow-up beyond IBD, IBS, and
autoimmune disease in a controlled way so that:

- the downstream pipeline becomes broadly reusable across projects,
- the gut paper can evaluate additional phenotype branches without becoming
  diffuse,
- Bayesian modeling is used where it adds value most clearly: deeper, sparser,
  lineage-level states.

## Current Starting Point

Already completed:

- IBD
- IBS
- Autoimmune

These three branches are already developed enough to anchor the current
manuscript rewrite.

The next branches should be chosen from the absorb omnibus results rather than
from intuition alone. Across depths 1 to 4, the strongest additional candidates
by depthwise Cramer's V growth and q-value support are:

- Cardiovascular disease: depth-1 to depth-4 Cramer's V `0.138 -> 0.315`
- Acid reflux: `0.100 -> 0.207`
- Seasonal allergies: `0.096 -> 0.204`
- Obesity: `0.079 -> 0.197`
- Lung disease: `0.058 -> 0.184`
- Migraine: `0.056 -> 0.182`
- CDI: `0.047 -> 0.177`
- Diabetes: `0.050 -> 0.173`

Not all of these are equally attractive for the manuscript. Some are strong
statistically but biologically heterogeneous or heavily confounded.

## Main Recommendation

Separate the next phase into analysis order and manuscript priority.

Analysis order:

1. Acid reflux
2. Seasonal allergies
3. Cardiovascular disease
4. Obesity
5. Migraine
6. Lung disease
7. Diabetes
8. CDI
9. Kidney disease only if needed later

Why this order:

- Acid reflux and seasonal allergies combine strong omnibus signal with large
  case counts and straightforward clinical framing.
- Cardiovascular disease has very strong absorb structure, but its phenotype is
  likely to carry heavier medication and comorbidity confounding, so it is best
  handled after the easier gastrointestinal and immune branches.
- Obesity is worth running early because it is common and historically central
  to gut-microbiome work, but interpretation is likely to be more diffuse.
- Migraine and lung disease are useful second-wave tests of generality.
- Diabetes and CDI are still interesting, but smaller case counts make them
  less efficient as immediate manuscript-expansion targets.

## Promotion Rules

### Promote to main text

A new phenotype branch should be promoted into the main Results only if all of
the following are true:

- the omnibus absorb signal remains strong across multiple depths,
- at least one lineage remains interpretable biologically at depths 2 to 4,
- frequentist and Bayesian directions broadly agree,
- the branch is not driven entirely by a few ultra-rare labels,
- there is a reasonable literature anchor for the ecological story,
- the phenotype is not so heterogeneous that the discussion becomes mostly
  caveats.

### Keep in supplement

Keep a branch in the supplement if:

- the omnibus signal is clearly real,
- the Bayesian layer is informative,
- but the lineage-level story is diffuse, heavily confounded, or not yet
  compelling enough for the main narrative.

This is the most likely destination for:

- cardiovascular disease
- obesity
- migraine
- lung disease

### Hold out of the manuscript for now

Do not write up a branch yet if:

- the signal is strong only at the omnibus level,
- the adjusted label-level results are unstable,
- the Bayesian layer points to tiny labels without a coherent parent lineage,
- or the phenotype definition is too ambiguous to support a biologically honest
  interpretation.

## Recommended Run Batches

Use the existing absorb run root:

- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive`

### Batch 1: highest manuscript-value expansion

Outcomes:

- `Acid_reflux`
- `Seasonal_allergies`
- `Cardiovascular_disease`
- `Obesity`

Suggested command:

```bash
Rscript /Users/pgajer/current_projects/gut_microbiome/scripts/run_absorb_label_followup.R \
  --run-dir /Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive \
  --output-dir /Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive/round2_label_followup_batch1 \
  --outcomes Acid_reflux,Seasonal_allergies,Cardiovascular_disease,Obesity \
  --depths 1,2,3,4 \
  --min-label-n 50 \
  --min-case-in-label 5 \
  --posterior-draws 4000
```

Expected decision after Batch 1:

- choose zero to two of these branches for main-text promotion,
- keep the rest as supplement candidates.

### Batch 2: generality and second-wave phenotypes

Outcomes:

- `Migraine`
- `Lung_disease`
- `Diabetes`
- `CDI`

Suggested command:

```bash
Rscript /Users/pgajer/current_projects/gut_microbiome/scripts/run_absorb_label_followup.R \
  --run-dir /Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive \
  --output-dir /Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive/round2_label_followup_batch2 \
  --outcomes Migraine,Lung_disease,Diabetes,CDI \
  --depths 1,2,3,4 \
  --min-label-n 50 \
  --min-case-in-label 5 \
  --posterior-draws 4000
```

Expected decision after Batch 2:

- retain the clearest branch as a supplement case study if useful,
- otherwise summarize these branches in one supplement overview table.

### Batch 3: conditional and low-priority extension

Only run if a specific clinical reason emerges:

- `Kidney_disease`
- any other phenotype that becomes especially interesting after Batch 1 or 2

This batch should be driven by scientific need, not by completeness for its
own sake.

## Standard Outputs To Review For Each New Branch

For every batch, review:

- `batch_summary.md`
- per-outcome `*_followup_summary.md`
- per-outcome `*_top_hits.tsv`
- per-outcome `*_depth_summary.tsv`
- per-outcome `*_adjusted_or_concordance.png`

Minimum questions to answer before manuscript promotion:

1. Does the branch have a coherent parent-to-child lineage story?
2. Does the Bayesian layer clarify deeper states rather than simply multiplying
   weak hits?
3. Are the strongest signals biologically plausible and literature-consistent?
4. Would the branch improve the paper's scientific focus, or only increase its
   breadth?

## Manuscript Strategy

Preferred structure after expansion:

- main text:
  - IBD
  - IBS
  - Autoimmune
  - optionally one or two additional branches from Batch 1
- supplement:
  - remaining Batch 1 branches
  - selected Batch 2 branches
  - cross-outcome Bayesian/frequentist concordance overview
  - complete label-level result tables

The paper should stay selective. The goal is not to narrate every phenotype
with a significant omnibus test. The goal is to show that the absorb dCST hierarchy
and Bayesian follow-up form a reusable, biologically interpretable downstream
analysis framework.

## Immediate Next Step

Run Batch 1 first.

If Batch 1 yields:

- one clearly interpretable branch:
  add it to the main manuscript and keep the rest supplementary.
- two clearly interpretable branches:
  choose the one that best broadens the paper beyond inflammatory disease for
  the main text and keep the other in supplement unless the manuscript still
  reads tightly.
- no clearly interpretable branches:
  keep the current main text focused on IBD, IBS, and autoimmune disease and
  present the broader Bayesian expansion as a supplement-only framework result.
