# Reviewer-Style Journal-Prep Review for the Gut dCST Manuscript

Date: 2026-04-11

Scope reviewed:

- Main manuscript:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- Supplement:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`
- Rendered PDF:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/build/gut_application_paper.pdf`
- Generic operating framework:
  `/Users/pgajer/current_projects/CT_clearance/docs/incremental_submission_style_paper_development_strategy.md`

Overall judgment:

The manuscript is scientifically promising and already far more journal-facing than an internal analysis memo, but it is still not in the cleanest submission-oriented state. The main remaining issues are editorial rather than computational: the narrative center is broader than the strongest evidence base, the paper still narrates some of its own project history, and several main-text sections would be stronger if moved to the supplement.

## Findings Ordered by Severity

### 1. Major: the main narrative is still broader than the strongest evidence base

The manuscript currently presents five phenotype branches as part of one main-text story: IBD, IBS, autoimmune disease, acid reflux, and seasonal allergies. That breadth makes the paper read more like a portfolio of completed analyses than a tightly centered submission draft. The strongest external and biological support is clearly in the inflammatory branch, especially IBD, with more limited support for the rest. A submission-oriented paper should either narrow the main text to IBD plus one or two secondary branches or explicitly demote the other branches to supplement.

Key locations:

- [`gut_application_paper.tex#L186`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L186)
- [`gut_application_paper.tex#L370`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L370)
- [`gut_application_paper.tex#L468`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L468)
- [`gut_application_paper.tex#L644`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L644)
- [`gut_application_paper.tex#L694`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L694)
- [`gut_application_paper.tex#L744`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L744)

### 2. Major: development-history framing remains in the scientific narrative

The paper still repeatedly tells the reader how manuscript branches were organized rather than simply presenting the science. Phrases about pure dCST analyses being “retained” as supplementary reference results, not being “repeated here,” or not being “carried as a second full results stream” read like editorial history from an evolving project rather than submission prose. This makes the manuscript feel merged from multiple draft logics.

Key locations:

- [`gut_application_paper.tex#L184`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L184)
- [`gut_application_paper.tex#L239`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L239)
- [`gut_application_paper.tex#L340`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L340)
- [`gut_application_paper.tex#L403`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L403)
- [`gut_application_paper.tex#L473`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L473)
- [`gut_application_paper.tex#L819`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L819)

### 3. Major: title, abstract, and conclusion still imply a broader paper than the central evidence supports

The title and abstract still read like a general disease-wide gut dominance paper, while the manuscript’s best-supported story is more specific: absorb dCSTs are most convincing in the inflammatory branch, with IBD as the anchor phenotype and broader generalization still more selective. The conclusion moves in the right direction, but the front end of the paper still creates a larger umbrella than the current reviewer-facing center should carry.

Key locations:

- [`gut_application_paper.tex#L144`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L144)
- [`gut_application_paper.tex#L173`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L173)
- [`gut_application_paper.tex#L204`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L204)
- [`gut_application_paper.tex#L871`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L871)

### 4. Moderate: methods are reproducible, but the manuscript still mixes scientific method with project framing

The core methods are now reasonably reproducible, especially for phenotype definition, absorb dCST hierarchy construction, Bayesian follow-up, and external validation. However, parts of Methods still read as manuscript-management explanation rather than pure scientific method, especially around pure-versus-absorb placement and local script ecosystem framing. The scientific object of each test is clearer than before, but the section could still be cleaner and more journal-efficient.

Key locations:

- [`gut_application_paper.tex#L325`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L325)
- [`gut_application_paper.tex#L340`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L340)
- [`gut_application_paper.tex#L346`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L346)
- [`gut_application_paper.tex#L383`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L383)
- [`gut_application_paper.tex#L891`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L891)

### 5. Moderate: figure and table economy can be tightened

The paper currently uses both Table 1 and Figure 2 to communicate the high-level cohort overview, and Table 2 functions largely as a manuscript-internal summary of results that are then narrated again phenotype by phenotype. For a submission draft, one of the overview display items should probably move to the supplement, and Table 2 is a strong candidate for demotion because it mainly summarizes material already described elsewhere.

Key locations:

- [`gut_application_paper.tex#L421`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L421)
- [`gut_application_paper.tex#L445`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L445)
- [`gut_application_paper.tex#L451`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L451)

### 6. Moderate: the Discussion still partly reads as a project summary rather than a submission discussion

The Discussion has improved substantially, but some language still reads like editorial triage of the project rather than journal-facing scientific interpretation. Phrases such as “stories are coherent but not equally mature” and the broad comparative sweep across multiple secondary phenotype branches keep the section closer to an internal synthesis memo than to a crisp discussion built around one central result and its implications.

Key locations:

- [`gut_application_paper.tex#L825`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L825)
- [`gut_application_paper.tex#L840`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L840)
- [`gut_application_paper.tex#L847`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L847)

### 7. Minor: the “Remaining phenotypes” section belongs outside the main text

The section naming cardiovascular disease, migraine, lung disease, and other branches does useful project-management work, but it is not helping the paper read like a clean submission manuscript. If those branches are not developed in this paper, they should live in supplement or working notes rather than in a main-text subsection.

Key locations:

- [`gut_application_paper.tex#L789`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex#L789)

### 8. Minor: the supplement is useful but should become more curated and less summary-like over time

The supplement is now carrying appropriate material, especially phenotype definitions, pure-versus-absorb divergence, selected IBD lineage detail, and the rebuilt-cohort versus direct AGP-derived label-transfer validation comparison. The current S4 interpretation column is helpful, but as the manuscript tightens it should evolve toward a stable reviewer-facing support document rather than a running summary of current project interpretations.

Key locations:

- [`gut_application_paper_supplement.tex#L114`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex#L114)
- [`gut_application_paper_supplement.tex#L124`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex#L124)

## Placement Decision at the Current Manuscript State

### Main text

At the manuscript’s current state, the following belong in the main text:

- the deterministic/absorb dCST concept and rationale;
- the AGP gut cohort definition and primary methods;
- the absorb dCST hierarchy overview;
- IBD as the anchor phenotype;
- the external-validation story for the inflammatory branch;
- at most one or two secondary phenotype branches, with IBS and autoimmune disease the strongest candidates if breadth is retained.

### Supplement / appendix

The following belong in the supplement by default:

- pure dCST reference analyses;
- detailed pure-versus-absorb comparison material;
- dense Bayesian summary tables;
- the full rebuilt-cohort versus direct AGP-derived label-transfer external-validation comparison;
- additional phenotype branches not chosen as core stories, especially acid reflux and seasonal allergies if the paper is narrowed;
- any retained high-level summary table that duplicates the results narrative.

### Working notes only

The following should remain in notes unless the paper’s scope changes:

- cardiovascular disease and other side branches not chosen for the main paper;
- tuning history and branch-selection logic;
- incomplete validation cohorts or speculative next-step analyses;
- manuscript-organization notes that explain why material was moved between sections.

## Journal-Prep Bottom Line

The manuscript is close enough to be worth tightening now, not later. The next pass should not primarily add more science. It should narrow and stabilize the story so that the paper reads like one reviewer-facing submission draft rather than a well-organized synthesis of several promising analysis branches.
