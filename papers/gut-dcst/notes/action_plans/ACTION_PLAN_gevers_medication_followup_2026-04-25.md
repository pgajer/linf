# Action Plan: Gevers Medication Follow-up (2026-04-25)

## Goal
Add an exploratory medication-stratified Gevers analysis that uses currently available sample-level ENA medication flags without overclaiming adjustment quality.

## Rationale
Gevers 2014 is the only currently analyzed IBD validation cohort in the local project state with usable sample-level medication metadata joined to the analyzed stool subset. Halfvarson and HMP2 do not currently have populated medication covariates in their harmonized analysis tables.

## Available flags in the analyzed Gevers stool subset
- `antibiotics`
- `mesalamine`
- `steroids`
- `biologics`

`biologics` has no positive samples in the analyzed subset and will therefore be excluded from the exploratory follow-up.

## Planned contrasts
1. `Crohn vs healthy` within explicitly medication-negative samples for:
   - antibiotics
   - mesalamine
   - steroids
2. `Medication-positive vs medication-negative` within Crohn cases for:
   - antibiotics
   - mesalamine
   - steroids

## Planned outputs
- Extend `build_ibd_reviewer_sensitivity_assets.R` with a Gevers medication follow-up builder.
- Generate manuscript-facing table assets:
  - `TABLE_S13E_gevers_medication_followup.tsv`
  - `TABLE_S13E_gevers_medication_followup.tex`
- Add `Supplementary Table S13E` to the supplement.
- Update the Gevers Results/Discussion text in the main manuscript to summarize:
  - medication-negative subset transfer remains null;
  - antibiotic exposure reshapes case-side transferred labels;
  - steroid follow-up shows a rebuilt-only exploratory signal;
  - mesalamine follow-up is weaker.

## Interpretation guardrails
- Treat the analysis as exploratory because medication fields are incomplete and not jointly modeled with disease activity or location.
- Do not describe this as a formal medication-adjusted model.
- Emphasize that the key question is whether the null Gevers transfer is reversed in medication-negative samples; if not, medication alone is unlikely to explain the portability boundary.
