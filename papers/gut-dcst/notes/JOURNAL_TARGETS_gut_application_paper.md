# Journal Targets: Gut DCST Application Paper

Date: 2026-04-03

## Why This Note Exists

The current planning stack records preprint timing and leaves `journal priority`
as an open decision, but it does not yet contain a concrete journal shortlist
for the gut DCST manuscript.

Relevant earlier notes:

- `/Users/pgajer/current_projects/linf/notes/project/CODEX_TASK_OWNERSHIP_2026-03-25.md`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/handoffs/HANDOFF_linf_application_paper_2026-03-25.md`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/MANUSCRIPT_SCAFFOLD_gut_application_paper.md`

## Current Manuscript Profile

At the moment, the paper reads most naturally as:

- a methods-application paper
- a gut-microbiome use case for the `linf` / DCST framework
- a statistically careful association paper rather than a causal paper
- a manuscript for mathematically literate and computational microbiome readers
- a paper whose fit improves substantially if validation and reproducibility are
  made prominent

This matters because the best venue depends on whether we want to present the
paper primarily as:

1. a microbiome paper with a strong new analysis framework, or
2. a computational/methods paper with a strong gut application

## Decision Criteria

When choosing a target journal, the main questions should be:

1. Does the journal welcome methods-plus-application papers, or does it expect
   deeper mechanistic biology?
2. How strong does external validation need to be before submission?
3. Will a mathematically explicit framing help or hurt fit with the readership?
4. How much benchmarking against alternative community-typing approaches do we
   need?
5. Are strong reproducibility expectations compatible with our current code/data
   release plan?
6. Do we want the first submission to optimize for prestige, fit, or speed?

## Recommended First-Pass Shortlist

### 1. `mSystems`

Why it fits:

- broad systems microbiology scope
- explicitly hospitable to tools and techniques
- likely to tolerate a balanced methods-plus-application story better than a
  mechanism-heavy microbiome journal
- good match if we want the paper to stay readable to both microbiome and
  quantitative audiences

What would strengthen the submission:

- a crisp statement of method novelty
- clean reproducibility story
- at least one coherent validation section

Main risk:

- if the manuscript reads as a large exploratory association screen without a
  sufficiently strong general methodological message

Working recommendation:

- best current first-choice venue if we keep the present balance between method,
  microbiome application, and interpretability

### 2. `PLOS Computational Biology`

Why it fits:

- good home for computational methods that generate real biological insight
- attractive if we want to foreground the deterministic DCST framework rather
  than the gut-disease application alone
- strong fit if the paper clearly explains why this framework improves on
  clustering-based typing

What would strengthen the submission:

- broader methodological framing beyond AGP
- explicit comparison against alternative typing approaches
- strong code/software visibility

Main risk:

- if the paper still reads primarily as a microbiome results paper rather than a
  computational-biology contribution

Working recommendation:

- strongest alternative if we deliberately make this a method-first paper

### 3. `Microbiome`

Why it fits:

- directly aligned readership
- high visibility for host-associated microbiome work
- plausible home if the manuscript matures into a stronger validation-backed gut
  microbiome story

What would strengthen the submission:

- stronger external validation
- tighter biological interpretation of the most credible signals
- careful explanation of why DCSTs add more than descriptive association mining

Main risk:

- the official scope explicitly says the journal is especially interested in
  work that goes beyond descriptive omics surveys and includes experimental or
  theoretical support for proposed microbiome functions, so a mostly
  association-driven paper could be judged too observational

Working recommendation:

- ambitious stretch target after stronger validation and a tighter biologically
  oriented framing

### 4. `The ISME Journal`

Why it fits:

- strong microbial ecology readership
- scope includes community ecology, theoretical advances, and bioinformatics
  linked to ecological questions
- possible fit if we frame DCSTs as a general microbial community-typing and
  community-ecology contribution

What would strengthen the submission:

- broader ecological significance beyond one disease-screening use case
- clearer general lessons about microbial community structure
- stronger emphasis on generalizable theory, not just AGP results

Main risk:

- high bar for broad ecological interest and impact
- weaker fit if the paper remains centered on one gut-disease application story

Working recommendation:

- realistic only if we intentionally broaden the paper beyond an application
  manuscript

### 5. `NAR Genomics and Bioinformatics`

Why it fits:

- strong methods and reproducibility orientation
- good option if we emphasize software, benchmarking, and reusable workflow
- attractive fallback if we want the paper judged as computational infrastructure
  rather than primarily as microbiome biology

What would strengthen the submission:

- explicit benchmarking or comparison framework
- clearly deposited code/data
- sharper description of the reusable computational deliverable

Main risk:

- less microbiome-specific readership
- weaker fit if we do not lean into method, FAIRness, and software-style
  reproducibility

Working recommendation:

- best fallback if we choose a reproducibility-forward methods identity

## Working Recommendation

If we had to choose a direction now, the most sensible order is:

1. `mSystems` as the best fit for the current manuscript shape
2. `PLOS Computational Biology` if we decide to push the method-first framing

Hold as stretch or conditional options:

- `Microbiome` after stronger external validation and stronger biological
  framing
- `The ISME Journal` only if we broaden ecological significance
- `NAR Genomics and Bioinformatics` if we deliberately package this as a more
  explicit methods/reproducibility paper

## Framing Changes By Target

If targeting `mSystems`:

- lead with deterministic community typing
- emphasize interpretability and reproducibility
- keep disease associations central but not overstated

If targeting `PLOS Computational Biology`:

- lead with computational novelty
- compare against clustering-based alternatives
- make generalizability beyond AGP more explicit

If targeting `Microbiome`:

- lead with biologically interpretable gut-microbiome findings
- foreground validation and careful disease-specific interpretation
- avoid sounding like a broad exploratory screen

If targeting `The ISME Journal`:

- frame the work as a microbial community ecology contribution
- highlight general principles rather than only disease results

If targeting `NAR Genomics and Bioinformatics`:

- foreground FAIRness, reusable code, benchmarking, and workflow design

## Open Decisions For Pawel

1. Do we want the first journal round to optimize for fit or for stretch?
2. Do we want the paper to read primarily as microbiome biology or as
   computational methodology?
3. How much external validation do we want before the first journal submission?
4. Do we want a preprint before the first journal submission, or in parallel?

## Sources Checked On 2026-04-03

- [Microbiome aims and scope](https://link.springer.com/journal/40168/aims-and-scope)
- [PLOS Computational Biology journal information](https://journals.plos.org/ploscompbiol/s/journal-information)
- [The ISME Journal guide to authors](https://www.nature.com/documents/ismej-gta.pdf)
- [ASM `mSystems` journal description](https://asm.org/webinars/msystems-webinar-series)
- [NAR Genomics and Bioinformatics scope and criteria](https://academic.oup.com/nargab/pages/scope_and_criteria)
