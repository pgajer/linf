# Action Plan: Pre-Submission Revision for mSystems

Date: 2026-04-11

## Goal

Bring the gut absorb dCST manuscript materially closer to a submission-ready
state by fixing the highest-impact scientific, structural, and presentation
issues identified in the pre-submission review.

## Execution plan

1. Checklist and scope lock
   - Status: completed
   - Write down the exact blockers and choose a discovery-focused absorb
     framing until absorb external validation is complete.

2. Methods and reproducibility upgrade
   - Status: completed
   - Add phenotype definition rules and sample counts for focal outcomes.
   - Add a reproducible description of the frequentist and Bayesian follow-up
     models, including covariates and decision summaries.

3. Figure and table cleanup
   - Status: completed
   - Build a real absorb overview figure from the depth-scan outputs.
   - Build a real IBD lineage figure from the absorb follow-up outputs.
   - Remove the external-validation placeholder from the main text.
   - Fix table numbering and PDF presentation issues.

4. Supplement structure
   - Status: completed
   - Create a real supplementary-materials source file.
   - Populate it first with phenotype definitions/sample counts and a compact
     pure-vs-absorb depth-change summary.

5. Submission-facing manuscript pass
   - Status: completed
   - Remove project-tracking language.
   - Align claims with evidence actually shown in the paper.
   - Keep the main story centered on the strongest absorb findings.

6. Verification
   - Status: completed
   - Rebuild the manuscript PDF.
   - Inspect the rendered pages for title-page, figure, and table correctness.

## Success criteria for this pass

- No main-text placeholders remain.
- No visible build stamp remains in the submission PDF by default.
- The manuscript no longer promises validation/robustness results that are not
  actually shown.
- Phenotype definitions, sample counts, and Bayesian methods are explicitly
  documented.
- A supplementary-materials file exists in the paper workspace.

## Expected remaining limitation after this pass

- Absorb-specific external validation may still remain an outstanding scientific
  strengthening step if it is not rerun in this revision cycle.

## Completed outputs

- Main manuscript revised at
  `papers/gut-dcst/manuscript/gut_application_paper.tex`.
- Real absorb overview and IBD-lineage figures added under
  `papers/gut-dcst/assets/figures/`.
- Supplement source added at
  `papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`.
- Main manuscript and supplement PDFs rebuilt and visually checked.
