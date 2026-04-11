# Gut dCST Manuscript Operating Note

Date: 2026-04-11

This note adapts the generic manuscript-development framework in
`/Users/pgajer/current_projects/CT_clearance/docs/incremental_submission_style_paper_development_strategy.md`
to the gut dCST paper centered on
`/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`.

## Reusable operating principles extracted from the generic strategy

1. Maintain one canonical manuscript.
The paper should live in one primary source file, with one curated supplement, rather than multiple near-final branches.

2. Separate stable science from exploratory work.
Every completed analysis must be placed into one of three buckets before it is written up: main text, supplement, or working notes only.

3. Calibrate claims before integration.
No result should enter the manuscript until its strongest justified claim and its strongest caveat are both known.

4. Keep the manuscript reviewer-facing at all times.
The paper should read as if an external reviewer could read it tomorrow, which means avoiding internal project history, placement narration, and local workflow framing in the scientific narrative.

5. Update the paper in controlled increments.
Each round of manuscript changes should be small enough that the narrative, claim strength, and figure/table economy can be reviewed immediately.

6. Give every figure and table a job.
Main-text display items must earn space by carrying a central scientific point; otherwise they belong in the supplement or should be deleted.

7. Rebuild and visually inspect after each meaningful change.
The rendered PDF is part of the manuscript state, not a later packaging step.

8. Keep side branches out of the main manuscript until they are genuinely central.
Exploratory phenotype branches, method tuning, and project-management summaries should live in notes, not in the paper.

## Gut-paper-specific translation

### 1. Main biological question

The gut paper should be centered on one question:

Can absorb dCSTs provide a biologically interpretable and externally supportable hierarchy of gut microbiome dominance states, with the clearest evidence in inflammatory phenotypes?

Corollary questions such as Bayesian stabilization of rare states, pure-versus-absorb comparison, and external label transfer are important, but they should support that central question rather than compete with it.

### 2. Canonical manuscript files

- Main manuscript:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- Canonical supplement:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`
- Main build script:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_gut_application_paper_pdf.sh`
- Supplement build script:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts/build_gut_application_paper_supplement_pdf.sh`

### 3. What belongs in the main text by default

The main text should retain only the following:

- the biological and conceptual motivation for deterministic gut community states;
- the AGP discovery cohort and phenotype definitions needed to understand the paper;
- the absorb dCST hierarchy as the primary analytical framework;
- the main cohort overview that establishes why the absorb hierarchy matters;
- IBD as the anchor phenotype branch;
- at most one or two secondary phenotype branches that materially strengthen the paper rather than broaden it for its own sake;
- external validation only to the extent that it calibrates the central inflammatory claim;
- the main limitations, especially around portability, self-reported phenotypes, and cross-cohort generalization.

### 4. What belongs in the supplement by default

The supplement should hold the following unless they become central:

- pure dCST association results and pure-versus-absorb comparison tables;
- expanded phenotype-screen summaries;
- dense Bayesian follow-up tables;
- full rebuilt-versus-frozen external-validation matrices;
- deeper lineage inventories that support but do not anchor the main argument;
- additional phenotype branches such as acid reflux, seasonal allergies, and other non-IBD examples unless they are explicitly retained as one of the few main-text secondary stories;
- expanded methods details that aid reproducibility but would interrupt the main narrative.

### 5. What should remain in working notes only

The following should stay out of the manuscript unless their role changes materially:

- cardiovascular, migraine, lung disease, and other side branches not chosen for the paper's final story;
- method-tuning history;
- partial external-validation runs;
- abandoned framing options;
- manuscript-organization notes;
- analysis branches that are still changing in interpretation.

### 6. Claim-calibration rules for the gut paper

- IBD may be described as the clearest and best-supported phenotype branch.
- IBS and autoimmune disease can be described as secondary or supportive branches only if they remain in the main text.
- Acid reflux and seasonal allergies should be described as exploratory or hypothesis-generating unless they become central and externally supported.
- External validation should be described as strongest for inflammatory framework reproducibility, not as universal portability of one fixed label dictionary.
- Bayesian results should be described as stabilizing or extending the readout for rarer states, not as replacing the main inferential framework.

### 7. Language rules specific to this manuscript

Do not use the following styles in the main text:

- editorial placement narration such as “retained as supplementary reference analyses” when the sentence can simply state what is analyzed in the paper;
- merged-project language such as “developed in parallel,” “carried forward,” or “the current branch”;
- notebook-style transitions such as “we also ran,” “for completeness,” or “the script generated”;
- local project framing that tells the reader how the manuscript evolved instead of what the science shows.

### 8. Immediate operating rule for future additions

Before any new gut-paper result is integrated, decide all three of the following first:

1. Is it central to the current submission story?
2. If not central, is it a supplement item or only a working-note item?
3. Can the result be summarized in one paragraph and supported by one figure or table or fewer?

If the answer to the third question is no, the result is not ready for main-text integration.

### 9. Recommended current story center

At the manuscript's current stage, the default submission-oriented center should be:

- absorb dCST framework in the AGP gut cohort;
- IBD as the anchor biological example;
- external validation of the inflammatory branch;
- one or two secondary branches at most.

This means breadth should be added only when it clearly strengthens the paper rather than merely documents completed analyses.
