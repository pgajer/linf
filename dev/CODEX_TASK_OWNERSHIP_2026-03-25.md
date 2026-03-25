# Codex Task Ownership and Delegation

Date: 2026-03-25

This document defines which parts of the `linf` gut-application and
disease-review program I should own directly, which parts are good subagent
work, and where Pawel's decisions are most valuable.

## Work I Should Own Directly

These tasks benefit from central judgment and continuity, so I should do them
myself rather than fully delegate them.

### 1. Package-paper boundary management

- review `README`, vignettes, and paper assets together
- keep the package story conservative and CRAN-friendly
- prevent the stratified `agp_gut` demo from being overstated

### 2. Manuscript claim shaping

- convert report results into a clear claim hierarchy
- separate:
  - descriptive results
  - adjusted associations
  - sensitivity-supported findings
  - speculative mechanisms
- make sure the discussion is honest and mathematically readable

### 3. Cross-disease synthesis

- review disease-specific literature outputs together
- identify shared mechanisms and shared papers
- decide what belongs in the paper discussion versus supplementary review assets

### 4. Handoff and action-plan maintenance

- keep the planning documents current
- make sure future Codex threads inherit the latest project state

### 5. Final review at each major phase

- review package-doc changes before close-out
- review literature-review outputs after subagents return
- review manuscript sections as they are drafted

## Work That Is Good To Delegate

These tasks are concrete, parallelizable, and easy to review centrally.

### Literature and review tasks

- disease-specific top-paper extraction
- DOI/link/full-text collection
- shared-paper overlap detection
- first-pass disease summaries
- figure inventory from review papers
- table extraction from papers and supplements

### Analysis support tasks

- locating specific metadata fields and code paths
- building candidate figure lists from report outputs
- validation-cohort scouting
- bibliography normalization and citation cleanup

### Documentation support tasks

- first draft of disease briefs
- first draft of glossary/explainer sidebars
- first pass at manuscript section summaries

## Work That Should Stay With Pawel

These are judgment calls where Pawel's preference or scientific positioning is
the main constraint.

- final target audience emphasis:
  - more mathematical
  - more microbiome-facing
  - balanced
- arXiv timing relative to CRAN
- journal priority
- appetite for releasing a first paper before external validation is added
- package-level decisions that materially affect public presentation

## Recommended Collaboration Model

### Phase 1: Central planning

I own:

- revised handoff documents
- package/paper split recommendation
- manuscript scaffold

### Phase 2: Parallel collection

Subagents own:

- disease-specific evidence gathering
- figure inventory
- citation and link collection

I own:

- prompt design
- merge/review/correction

### Phase 3: Writing

I own:

- section architecture
- discussion framing
- mathematical/contextual explanation

Subagents can support:

- first-pass section drafts
- table assembly
- source checking

### Phase 4: Review

I own:

- final pass on claims, caveats, and consistency
- integration across package docs, application paper, and review assets

## Concrete Tasks I Can Start Taking On Now

1. Rewrite the planning stack so it matches the current project state.
2. Prepare the handoff for a dedicated `linf` application-paper thread.
3. Update the gut disease-review action plan around the clean split
   bibliography.
4. Review the paper as it progresses through phases and flag overclaims,
   missing context, and places where mathematician-oriented explanation is
   needed.
5. Design and supervise the subagent pass for disease reviews and common assets.
6. Help shape figure strategy:
   - what can be reused from current outputs
   - what should be redrawn
   - what should come from literature only when licensing permits
7. Help convert the literature review into paper discussion assets:
   - mechanism summaries
   - alternative hypotheses
   - future-research paragraphs

## Suggested Near-Term Sequence

1. Finish the new handoff/action-plan docs.
2. Launch the dedicated `linf` paper thread from the new handoff.
3. Update the disease-review workflow and run the next structured review pass.
4. Bring the review outputs back into the manuscript discussion plan.

## Definition Of Success

Success is not just "more documents" or "more papers collected." Success means:

- the package docs are clearer,
- the paper claims are better calibrated,
- the review assets are reusable,
- and the project stays easy to hand off without losing context.
