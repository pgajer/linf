# Action Plan: External IBD Transfer Refactor

Date: 2026-04-24

This plan implements the rewrite schema in:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/MANUSCRIPT_REWRITE_SCHEMA_external_ibd_transfer_2026-04-24.md`

## Goal

Refactor the gut dCST manuscript so that:

- AGP is used primarily for hierarchy derivation and ecological structure learning;
- clinically annotated external IBD cohorts carry the main disease-association narrative;
- CD and UC are split explicitly where the external cohorts allow it;
- AGP self-reported disease associations become exploratory context rather than the primary biological evidence.

## Scope for this pass

This pass focuses on the main manuscript only:

- `manuscript/gut_application_paper.tex`

It may adjust captions and section ordering, but it does not require a new figure-construction round unless strictly necessary for coherence.

## Stage 1. Freeze and document the pre-refactor state

Tasks:

- archive the current manuscript `.tex` and rebuilt `.pdf` as a milestone snapshot;
- write a rewrite-schema note;
- write this action plan.

Done criteria:

- the pre-refactor manuscript is preserved and easily retrievable;
- the refactor logic is documented before editing the live source.

## Stage 2. Rebuild the front matter around the new paper

Tasks:

- rewrite the Abstract around AGP-derived structure plus external clinical transfer;
- rewrite the closing paragraphs of the Introduction so the paper’s main question is external IBD transfer rather than AGP phenotype discovery;
- revise any preview language that still centers AGP self-reported IBD as the primary result.

Done criteria:

- title, abstract, and introduction all describe the same paper;
- the reader expects an external-IBD-centered results section.

## Stage 3. Reorder the Results around external IBD cohorts

Tasks:

- shrink the current AGP outcome-wide discovery overview into a short hierarchy / prioritization section;
- convert the current AGP IBD subsection into contextual support rather than the main evidence block;
- elevate AGP-derived external transfer results into the central Results narrative;
- explicitly organize the external IBD story around:
  - Halfvarson 2017,
  - HMP2 / IBDMDB,
  - Gevers 2014,
  - PRJEB84421 only if needed as a secondary inflammatory complement.

Done criteria:

- the main Results read as an external IBD transfer paper;
- AGP disease associations no longer dominate the narrative.

## Stage 4. Split IBD by subtype where admissible

Tasks:

- make CD and UC explicit in the Halfvarson discussion;
- make pooled IBD versus subtype-specific support explicit in HMP2;
- treat Gevers as the Crohn-focused anchor;
- ensure wording distinguishes what was shown for pooled IBD versus CD or UC separately.

Done criteria:

- subtype-specific claims are stated exactly and only where supported by the current analyses.

## Stage 5. Remove non-central phenotype branches from the main Results

Tasks:

- remove or sharply compress the IBS subsection;
- remove or sharply compress the autoimmune subsection;
- replace them with one sentence directing readers to exploratory supplementary material.

Done criteria:

- the main Results no longer read as a multi-phenotype survey;
- the IBD transfer story is not diluted.

## Stage 6. Reframe the Methods

Tasks:

- make AGP hierarchy construction the first methodological role of AGP;
- describe AGP phenotype coding as exploratory and self-reported;
- present frozen AGP-derived transfer into external cohorts as the main validation strategy;
- retain rebuilt-cohort validation as a secondary comparison mode.

Done criteria:

- Methods clearly distinguish structure learning from disease inference;
- a reader can understand why external cohorts are central to admissible IBD interpretation.

## Stage 7. Tighten Discussion and Conclusion

Tasks:

- center the Discussion on transferable representation of known IBD dysbiosis;
- explicitly state that the strongest disease-facing evidence comes from clinically annotated external cohorts;
- remove broad non-IBD framing from the main interpretive close;
- tighten the Conclusion around IBD-focused portability and interpretability.

Done criteria:

- the paper closes on a narrower but stronger claim;
- the Discussion reads like a submission draft rather than a project portfolio summary.

## Stage 8. Rebuild and verify

Tasks:

- rebuild the manuscript PDF;
- confirm the build succeeds and the updated structure is coherent;
- note any remaining issues that should be deferred to a later pass.

Done criteria:

- the refactored manuscript compiles cleanly;
- the result is reviewer-readable and structurally consistent.

## Execution order for this pass

1. Complete Stage 1.
2. Edit front matter.
3. Reorder and rewrite Results.
4. Reframe Methods.
5. Tighten Discussion and Conclusion.
6. Rebuild and verify.

## Non-goals for this pass

- no new cohort downloads;
- no new pooled cross-cohort meta-analysis;
- no major figure rebuild unless the manuscript becomes incoherent without it;
- no supplement rewrite beyond wording references from the main text.
