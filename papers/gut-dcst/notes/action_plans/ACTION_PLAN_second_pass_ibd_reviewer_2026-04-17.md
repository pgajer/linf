# Second-Pass IBD Reviewer Action Plan

Date: 2026-04-17

Scope:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/tables/`

Goal:

Address the second-pass IBD/gut-biologist reviewer recommendations before
submission-facing polishing. The emphasis is reviewer-proofing rather than
changing the central story: preserve the IBD-centered absorb-policy dCST
manuscript while making feature filtering, AGP phenotype limitations,
artifact-prone taxa, lineages containing `Faecalibacterium`, and external
validation boundaries harder to misread.

Actions:

- [x] Clarify that dCST dominance is determined after the stated retained-feature
  filtering, so labels are retained-feature dominance states rather than claims
  about every raw feature in the original table.
- [x] Bring the self-reported and non-mutually exclusive AGP phenotype caveat
  closer to the IBD discovery result, including limited control for disease
  subtype, activity, medication, and antibiotic exposure.
- [x] Add a short artifact/non-gut annotation caveat for labels such as
  `Unassigned`, `Eukaryota`, `Mitochondria`, `Pseudomonas`, and `Staphylococcus`,
  and make clear that these labels are retained for transparent unsupervised
  accounting but are not used as main biological anchors.
- [x] Add a lineage-profile table for the key IBD dCSTs so readers understand
  that a depth label records ordered retained-feature dominance context; a
  label containing `Faecalibacterium` does not imply that
  `Faecalibacterium` is increased overall in IBD.
- [x] De-emphasize rare `Morganella` and `Proteus` in the abstract while keeping
  them as explicitly exploratory rare-state findings in the Results.
- [x] Explain mixed taxonomic-rank labels as deepest resolved SILVA annotations
  available for retained features.
- [x] Calibrate "inflammatory" language so it is clear that the manuscript refers
  to disease-contextual inflammatory phenotypes and literature-supported
  dysbiosis patterns, not direct host inflammatory biomarker measurements.
- [x] Strengthen external-validation reporting by adding cohort context to the
  supplement and simplifying the main-text validation table so it is readable.
- [x] Rebuild and visually inspect the manuscript and supplement PDFs.
- [x] Commit and push the full staged change set.

Completion note:

The manuscript and supplement were rebuilt successfully after these edits. The
only remaining LaTeX diagnostics were minor underfull-box warnings in main-text
paragraphs, with no missing citations, overfull boxes, or build failures.

Expected outcome:

The paper should remain centered on IBD as the strongest absorb-dCST branch,
but it should read as more cautious and more biologically transparent to an IBD
reviewer. The main text should not overclaim clinical mechanism from
self-reported AGP metadata, and the supplement should carry the extra context
needed to make validation and taxonomic interpretation auditable.
