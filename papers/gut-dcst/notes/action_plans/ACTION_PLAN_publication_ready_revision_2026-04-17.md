# Publication-Readiness Revision Plan

Date: 2026-04-17

Target files:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/`

Goal:

Move the gut-dCST paper closer to an mSystems-ready submission draft by
addressing the reviewer-style concerns from the IBD/gut microbiome review pass.

Execution plan and status:

- [x] Change the title to `Dominant Community State Types Reveal Inflammatory Dominance Patterns in the Human Gut Microbiome`.
  Updated the main manuscript title, PDF metadata title, and supplement title.

- [x] Tighten the abstract so it reads as a scientific summary rather than an editorial inventory.
  Recentered the abstract on the biological problem, absorb-policy dCSTs, the IBD result, calibrated external validation, and the role of Bayesian follow-up.

- [x] Standardize terminology around `dCST`, `absorb-policy dCST`, `pure-policy dCST`, and `dominance-lineage`.
  Removed remaining all-caps `DCST` usage from the manuscript and supplement, defined dominance-lineage once in the Introduction, and used that wording consistently for deeper labels.

- [x] Clarify the operational absorb rule in Methods.
  Rewrote the construction description so the algorithmic idea is separated from implementation details and explains how low-frequency dominance sample sets are reassigned to retained named states.

- [x] Add case/control support for the most important rare IBD states.
  Added explicit IBD/control support for the rare enterobacterial examples and framed them as suggestive sparse-state signals rather than stand-alone biomarkers.

- [x] Run and report an IBD-focused control-definition sensitivity check addressing AGP self-report and non-mutually exclusive phenotype controls.
  Re-ran the IBD absorb-policy dCST dominance-lineage follow-up after excluding non-IBD controls carrying any other AGP phenotype string analyzed in this study (IBS, autoimmune disease, acid reflux, or seasonal allergies). Added the result to the IBD Results text and Supplementary Table S7.

- [x] Simplify Figure 1 using the CT-paper logic without importing the dense CT schematic wholesale.
  Rebuilt Figure 1 as a gut-facing dCST schematic covering a toy abundance table, dominance sample sets, pure versus absorb handling of small states, and depth-2 dominance-lineages.

- [x] Simplify Figure 2C so main-text phenotypes remain visually primary.
  Removed acid reflux and seasonal allergy curves from the main panel and retained IBD, IBS, and autoimmune disease as the main-text phenotype branches.

- [x] Improve Figure 3 so rare-state counts are visible and the IBD result reads as the flagship biological figure.
  Added IBD/control counts to the representative IBD-state estimates and sharpened the caption around common versus rare-state interpretation.

- [x] Add a compact main-text external-validation table and calibrate Figure 4 as supporting visualization rather than the only validation summary.
  Added a main-text validation table distinguishing rebuilt-cohort validation from direct AGP-derived label transfer, then revised the Figure 4 caption so it supports rather than overstates label portability.

- [x] Split the cramped supplement validation table into rebuilt-cohort and direct AGP-derived label-transfer summaries.
  Reorganized the supplement validation material into separate rebuilt and direct AGP-derived label-transfer tables for readability.

- [x] Clean Data and Code Availability formatting.
  Reworked the section into a compact reproducibility-facing list and used path formatting for long script names.

- [x] Rebuild manuscript and supplement PDFs.
  Rebuilt both PDFs with the manuscript build scripts after figure and source edits.

- [x] Visually inspect key rendered pages.
  Checked title/abstract, Figure 1, Figure 2, Figure 3, main validation table, Figure 4, and Supplementary Table S7 in rendered PDF pages.

Items not fully resolvable by manuscript editing alone:

- Future external-validation datasets may strengthen or weaken the story. The current manuscript should remain calibrated as strong for IBD framework reproducibility, moderate for selected AGP-derived label portability, and limited for broad universal dictionary claims.
- The HMP2 pooled IBD comparison remains weaker than Halfvarson in the current processed subset. The manuscript treats it as useful but not decisive transfer support.
- IBS and autoimmune disease remain secondary branches until outcome-specific external validation is available.

Completion notes:

- The main manuscript should remain centered on IBD as the anchor branch.
- IBS and autoimmune disease remain secondary main-text branches.
- Acid reflux and seasonal allergies remain supplementary exploratory examples.
- The CT-paper nomenclature should be borrowed selectively; the gut paper should not adopt the full CT Figure 4 notation load in the main text.
