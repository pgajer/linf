# Pre-Submission Checklist for mSystems

Date: 2026-04-11

Purpose: convert the reviewer-style assessment into a concrete pre-submission
checklist ordered first by impact on editorial/reviewer acceptance and then by
execution effort.

## Tier 1: Submission blockers

- Replace all main-text placeholder figures with final figures or remove the
  sections that still depend on them.
- Remove visible draft-only build metadata from the submission PDF.
- Eliminate in-progress language such as "being extended", "in parallel",
  "Batch 1", and similar project-tracking phrasing from the manuscript body.
- Resolve the mismatch between the claims made in the Introduction/Abstract and
  the evidence actually shown in the paper, especially around external
  validation, taxonomy robustness, and contamination-aware filtering.
- Fix duplicate/confusing table numbering in the rendered PDF.

## Tier 2: High-impact scientific fixes

- Add explicit phenotype-definition rules for the AGP outcomes analyzed in the manuscript.
- Add per-phenotype sample counts:
  - total filtered samples
  - case counts
  - control counts
  - complete-case counts for adjusted models
- Expand the Bayesian methods description so another group could reproduce the
  analysis without reading the code first.
- Decide on the validation stance:
  - either rerun external validation in absorb mode and report it
  - or narrow the manuscript honestly to a discovery-focused absorb analysis
- Provide a real supplementary-materials structure rather than referring to a
  supplement that does not yet exist.

## Tier 3: Narrative tightening

- Keep the main paper centered on the strongest and most interpretable
  phenotype branches.
- Consider whether acid reflux and seasonal allergies should remain in the main
  text or be softened / shifted to supplementary analyses.
- Keep the pure dCST branch visible only as baseline context and
  supplementary-reference material.

## Tier 4: Presentation polish

- Ensure figures and tables appear in a clean journal-ready sequence.
- Check caption wording so captions explain the figure rather than duplicating
  panel text.
- Rebuild the PDF and visually inspect the first page, figure pages, and all
  table pages before submission.

## Recommended immediate strategy

- Make the paper internally honest and submission-like now by centering it as a
  discovery-focused absorb dCST manuscript.
- Replace the absorb overview and IBD lineage placeholders with real figures
  based on current absorb outputs.
- Remove the external-validation placeholder from the main text until absorb
  validation is complete.
- Add phenotype definitions, sample counts, and Bayesian model details in the
  manuscript and supplement.
