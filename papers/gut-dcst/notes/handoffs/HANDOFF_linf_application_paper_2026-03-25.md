# Handoff: linf Gut Application Paper

Date: 2026-03-25

This document is the current handoff for the `linf` gut application-paper
thread. It complements, and partly supersedes, the earlier seed-analysis
handoffs created on 2026-03-24.

## Purpose

Move the gut DCST work from "promising exploratory analysis" to a coherent
application-paper project that supports:

1. the `linf` CRAN submission,
2. a mathematically readable gut-application preprint, and
3. a disease-review library that can be reused in the paper discussion.

## High-Level Decision

Treat the project as three linked but distinct products:

1. **CRAN/package product**
   - stable, compact, reproducible, conservative in claims
   - uses bundled example data and short didactic narratives
2. **Application paper**
   - full-cohort AGP analysis with adjusted models, sensitivity analysis,
     depth-2 findings, and a clear claim hierarchy
3. **Disease-review library**
   - structured literature assets used to contextualize the paper
   - should support the paper, not block it

## Current State

### Package/docs

- `README.Rmd` already presents a compelling story, but the gut example carries
  more inferential weight than is ideal for a CRAN-facing README.
- `vignettes/linf-intro.Rmd` currently does too much. It starts as a methods
  tutorial and then becomes a gut disease-analysis vignette.
- `vignettes/linf-vaginal.Rmd` is the cleanest biologically grounded teaching
  document in the package.
- `data/agp_gut.rda` is a stratified, enriched demonstration dataset, not a
  population-representative cohort.

### Application-paper assets

- `dev/Phase1_Full_Cohort_DCST_Analysis_Report.pdf` already reports a
  30,290-sample full-cohort analysis with:
  - adjusted logistic regression
  - contamination sensitivity analysis
  - depth-2 associations
  - GreenGenes2 comparison
- `vignettes/articles/gut-dcst-disease-analysis.Rmd` is still a prototype
  article describing a 5,000-sample exploratory analysis with `eval = FALSE`.
- `dev/ACTION_PLAN_gut_dcst_paper.docx` still reads mostly as a "what to do
  next" plan, even though some of that work now exists in report form.

### Disease reviews

- The bibliography effort is now substantially cleaner and more useful than the
  original umbrella-bucket version.
- The current source of truth for literature assets is under
  `~/current_projects/gut_microbiome/outputs/dcst_analysis/disease_reviews/`.
- The most important files there are:
  - `CLEAN_ALL_BIBLIOGRAPHY_CANDIDATES.tsv`
  - `CLEAN_EXACT_SHARED_COMMON_ASSETS.tsv`
  - `CLEAN_SPLIT_SUMMARY.md`

## Strategic Decisions To Keep

1. Keep the vaginal example as the clean validation/example ecosystem.
2. Keep the gut example as the main application and public-interest driver.
3. Keep the disease-review effort, but tie it to paper claims and mechanisms.
4. Write for a mathematician:
   - define microbiome terms on first use,
   - explain what odds ratios and BH correction mean,
   - separate exploratory signal from causal or mechanistic claims,
   - explain why depth-2 DCSTs matter geometrically and biologically.

## Strategic Decisions To Change

1. **Split the package intro from the gut application story.**
   Recommended package structure:
   - `linf-intro`: method-first tutorial with toy data and a very light real-data
     touch
   - `linf-vaginal`: main biological validation vignette
   - gut article: package website / pkgdown article, not the pedagogical core

2. **Be explicit that `agp_gut` is illustrative.**
   The bundled gut dataset is stratified and enriched. It is appropriate for
   demonstrating the pipeline, but not for making strong prevalence claims.

3. **Bring the gut article into sync with the report.**
   The current article should either:
   - be rewritten as a "prototype/exploratory analysis" article, or
   - be replaced by a report-aligned article built around the 30,290-sample
     results.

4. **Do not let the full review program delay the paper.**
   Review work should first prioritize the phenotypes and disease families that
   appear in the adjusted AGP results.

## Immediate Objectives For The New Thread

### Objective A: Reconcile package, report, and paper story

Produce one unified statement of:

- what the package demonstrates,
- what the full AGP analysis establishes,
- what remains exploratory,
- what still needs validation.

### Objective B: Build a manuscript scaffold

Turn the current assets into a manuscript-ready structure with:

- claim hierarchy,
- figure plan,
- table plan,
- limits and caveats stated early and clearly.

### Objective C: Define package-paper boundaries

Decide exactly what lives in:

- `README`
- `linf-intro`
- `linf-vaginal`
- gut pkgdown article
- application paper

## Recommended Deliverables From The Next linf Thread

1. `papers/gut-dcst/notes/MANUSCRIPT_SCAFFOLD_gut_application_paper.md`
   - section-by-section paper outline
   - one-sentence claim for each section

2. `papers/gut-dcst/notes/CLAIM_INVENTORY_gut_application_paper.tsv`
   - claim
   - evidence source
   - exploratory vs adjusted vs sensitivity-supported
   - literature support needed

3. `papers/gut-dcst/notes/FIGURE_PLAN_gut_application_paper.md`
   - must-have figures
   - optional figures
   - which existing report figures can be reused
   - which should be redrawn

4. Updated package-doc plan
   - either a revised markdown action plan or direct edits to README/vignettes

## Recommended First Tasks In The New linf Thread

1. Read:
   - `papers/gut-dcst/archive/2026-03-24-phase1/Phase1_DCST_Report_by_Disease.pdf`
   - `vignettes/articles/gut-dcst-disease-analysis.Rmd`
   - `vignettes/linf-intro.Rmd`
   - `README.Rmd`
2. Compare the report claims against the current vignette/article language.
3. Draft a paper claim hierarchy:
   - methodological claim
   - biological association claim
   - robustness claim
   - limit/uncertainty claim
4. Propose the final package vignette split before editing prose.
5. Only then start rewriting README/vignettes/manuscript scaffold.

## Things The New Thread Should Avoid

- Do not present the bundled `agp_gut` subset as representative of AGP.
- Do not merge disease-review bucket names directly into AGP phenotype names
  without explanation.
- Do not make mechanistic claims in the paper without a clear literature anchor.
- Do not let the article continue to describe the project as a 5,000-sample
  "full analysis" if the real report already goes beyond that.

## Key Open Decisions For Pawel

These are the few decisions worth explicit confirmation:

1. Should `linf-intro` remain a combined tutorial, or become a mostly methods
   vignette?
   Recommendation: make it mostly methods.
2. Should the gut application paper go to arXiv before the CRAN submission, or
   in parallel?
   Recommendation: parallel is fine if the manuscript is already coherent and
   caveated.
3. How much validation is required before the first preprint?
   Recommendation: enough internal robustness and honest caveats for arXiv; an
   external cohort can be Phase 2 rather than a hard precondition.

## Inputs For Literature Context

When the manuscript needs disease context, the thread should rely on the clean
review assets under:

`~/current_projects/gut_microbiome/outputs/dcst_analysis/disease_reviews/`

Use those assets as structured background, not as a substitute for checking
which specific claims are needed by the AGP results.
