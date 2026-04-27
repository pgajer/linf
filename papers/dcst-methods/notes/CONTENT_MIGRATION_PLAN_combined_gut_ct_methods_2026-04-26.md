# Content Migration Plan: Combined Gut+CT dCST Methods Paper

## Purpose

This document translates the high-level outline into a practical migration plan.

The goal is to build a new methods-forward manuscript in:

`/Users/pgajer/current_projects/linf/papers/dcst-methods`

using three upstream sources:

- gut application manuscript:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- CT manuscript:
  `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`
- original Linf theory paper:
  `/Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf`

This is a **selective migration plan**, not a merge-everything plan.

## Governing editorial rule

Every migrated block must help support the new umbrella claim:

> dCSTs are a deterministic, interpretable microbiome-state representation that is useful across more than one biological niche and more than one downstream analysis mode.

If a section does not support that claim, it should stay out of the new paper or move to supplement.

## Source roles

### Gut manuscript role

The gut paper contributes:

- the practical motivation from clustering-defined enterotypes
- the deterministic taxonomic dCST framework in a large ecological background cohort
- the external transfer and portability logic
- the shared-taxonomy harmonization story
- the honest repeated-sample / subject-aware limitation story

### CT manuscript role

The CT paper contributes:

- a second microbiome niche
- a second data regime
- a stronger biologically structured application
- functional follow-up that demonstrates interpretability beyond taxonomic labels

### Original Linf paper role

The Linf paper contributes only:

- the minimal methods/philosophy framing needed to explain why the framework is deterministic and stable
- possibly one small conceptual figure idea if useful

It does **not** contribute full theory sections.

## Migration categories

## Category A: Move with moderate rewriting

These are high-value sections that should migrate into the new main paper after rewriting.

### A1. Gut methods core

Source:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`

Move:

- dCST definition at the microbiome-facing level
- absorb-policy explanation
- AGP local-QZA SILVA harmonization
- external cohort and transfer-mode design
- subject-aware sensitivity description

Rewrite goals:

- strip out gut-only framing from the method definition
- make the method sound general before the application-specific details begin

### A2. Gut main results skeleton

Source:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`

Move:

- AGP exploratory prioritization of IBD
- external transfer design logic
- Halfvarson as the clearest stool-based support
- HMP2 and Gevers as portability boundaries

Rewrite goals:

- shorten this relative to the current gut manuscript
- preserve the strongest transferable result
- keep the subject-aware limitation visible and honest

### A3. CT methods and cohort design

Source:

- `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

Move:

- labeled CT subset versus broader CT+reference cohort logic
- VOG-cluster and full-VOG analysis layers
- restriction of CT outcome tests to the labeled subset
- functional annotation follow-up logic

Rewrite goals:

- harmonize wording with the gut application methods section
- keep the dCST core method shared, while letting the CT application have its own feature universe and follow-up layer

### A4. CT main results skeleton

Source:

- `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

Move:

- lower-resolution VOG-cluster signal
- primary full-VOG outcome-associated branches
- clearance-oriented versus persistence-oriented interpretation
- taxonomic and functional follow-up

Rewrite goals:

- present CT as a second proof point for the framework, not as a competing full paper
- preserve biological richness without importing every CT-specific detail

### A5. Cross-application synthesis language

Source:

- not directly migratable as-is; must be newly written

Create:

- a new synthesis subsection in Results
- a new synthesis layer in Discussion

Purpose:

- make the paper read as one methods paper rather than two adjacent applications

## Category B: Move in reduced form

These components are worth keeping, but only in compressed form.

### B1. Linf conceptual framing

Source:

- `/Users/pgajer/current_projects/Linf_paper_arXiv_submission/Linf_paper.pdf`

Keep only:

- deterministic sample-wise state assignment
- within-sample rank-based invariance
- stability relative to adding/removing other samples
- contrast with clustering-defined labels

Do not carry over:

- long geometry exposition
- cube embedding sections
- broad non-microbiome generalization

### B2. Gut within-disease follow-up

Source:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`

Recommendation:

- keep only the strongest or most informative within-disease follow-up in main text
- likely keep Crohn location
- move UC extent and calprotectin details to supplement unless they are needed for one specific portability-limit point

### B3. CT robustness layers

Source:

- `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

Recommendation:

- keep capped-feature robustness in reduced form
- avoid letting the main paper become dominated by CT sensitivity variants

## Category C: Supplement-only material

These are useful, but should not burden the main paper.

### C1. Additional gut exploratory phenotype screens

Move to supplement:

- IBS
- autoimmune disease
- acid reflux
- seasonal allergies

### C2. Extended gut audit and sensitivity detail

Move to supplement:

- full overlap / assignment audit tables
- extended reviewer-sensitivity tables
- thresholded calprotectin details

### C3. Extended CT technical tables

Move to supplement:

- long annotation tables
- extra branch-level listings
- extended Bayesian outputs

### C4. Any residual Linf-theory note

If needed at all, put in supplement:

- a short methods note clarifying why deterministic rank-based assignment differs from clustering

## Category D: Leave out entirely

These components should not migrate unless a concrete later need emerges.

### D1. Large portions of Linf theory

Do not import:

- detailed geometric constructions
- topological/computational geometry ambitions
- wide generalization claims far beyond microbiome state labeling

### D2. Reviewer-driven local detail that no longer serves the new umbrella paper

Do not import into main text:

- highly localized gut follow-up detours
- narrow troubleshooting narratives
- manuscript-history rationale that only made sense during gut-paper iteration

### D3. CT-paper detail that belongs in a future standalone biology paper

Do not import:

- very fine-grained branch-by-branch storytelling that is more useful in a CT-specific biology framing

## Figure migration plan

## Figure class 1: Reuse with adaptation

### Gut overview figure

Current source:

- gut Figure 2 logic

Use in new paper as:

- the portability / transfer application figure

Needed changes:

- likely simplify and remove gut-specific local history
- make visual style compatible with the CT figure system

### Conceptual dCST schematic

Current source:

- gut Figure 1 concept

Use in new paper as:

- general methods figure

Needed changes:

- make it explicitly cross-application rather than gut-specific

## Figure class 2: Import and redesign

### CT application figure

Current source:

- main CT result figures in the CT manuscript

Use in new paper as:

- one or two principal CT results figures

Needed changes:

- align design language with the new methods paper
- emphasize framework output rather than CT-paper local narrative

## Figure class 3: New figure required

### Cross-application synthesis figure

Need:

- one figure that makes the combined paper feel unified

Possible content:

- deterministic framework at left
- gut portability proof point in the middle
- CT biological/functional proof point at right

or

- a property matrix:
  - deterministic assignment
  - portability
  - subject-aware robustness
  - functional interpretability

## Table migration plan

## Main-text tables

Recommended:

- one compact gut results table
- one compact CT results table
- possibly one high-level framework-properties table if needed

## Supplement tables

Recommended:

- retain most long-form association tables in supplement
- avoid repeating long tables already better communicated by figures

## Writing order

## Phase 1: Lock the content architecture

Do first:

1. finalize section map
2. decide what stays in main text for gut
3. decide what stays in main text for CT
4. identify which figures are reused, redesigned, or new

## Phase 2: Build shared methods core

Do second:

1. write a fresh dCST framework section
2. write application-specific methods subsections
3. define shared terms and notation

## Phase 3: Migrate and rewrite results

Do third:

1. migrate gut results in shortened form
2. migrate CT results in methods-paper form
3. write cross-application synthesis section

## Phase 4: Finish discussion and opening sections

Do last:

1. write Discussion after the main results are settled
2. write Introduction after the results architecture is stable
3. write Abstract last

## Concrete migration checklist

### Gut

- [ ] identify which current gut results paragraphs remain main-text worthy
- [ ] decide whether Crohn location is the only within-disease follow-up kept in main text
- [ ] mark which gut supplement items remain relevant to the umbrella paper

### CT

- [ ] identify the one or two CT results figures that best support the framework claim
- [ ] decide whether the VOG-cluster layer stays in main text or supplement
- [ ] decide how much functional follow-up belongs in main text

### Linf

- [ ] extract only the small conceptual statements worth reusing
- [ ] avoid importing any long-form theory section

### New paper

- [ ] create the main manuscript source
- [ ] create a figure inventory note
- [ ] create a source-to-section mapping note if needed

## Bottom line

The best version of this paper is not:

- gut plus CT plus Linf pasted together

It is:

- a new methods manuscript with a shared dCST core,
- one transfer/portability application in gut,
- one biologically richer functional application in CT,
- and only a minimal Linf-derived conceptual layer.
