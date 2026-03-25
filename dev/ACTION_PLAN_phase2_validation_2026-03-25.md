# Action Plan: Phase 2 Validation Track

Date: 2026-03-25

## Purpose

Run Phase 2 as a **parallel validation track** for the gut application paper.

The validation track should strengthen the paper without freezing the main
writing effort. The manuscript can continue to develop from the completed AGP
Phase 1 analysis while Phase 2 tests whether the strongest and most credible
claims replicate in an external cohort.

## Why This Runs In Parallel

Phase 1 and Phase 2 answer different questions:

- **Phase 1** establishes the discovery landscape in AGP.
- **Phase 2** asks which parts of that landscape generalize.

These should run in parallel because:

1. the paper architecture, claim hierarchy, and figures can already be built
   from Phase 1,
2. validation should target specific claims rather than block all writing,
3. the right validation cohort may support only a subset of AGP phenotypes.

So the correct workflow is:

- write the paper from Phase 1 now,
- validate the highest-value claims in parallel,
- upgrade, narrow, or caveat claims when Phase 2 results come in.

## Local Starting Point

Based on the current local PRIME assets:

- `PRJNA398089` is already present in the PRIME manifest as `LIBD-US-2018`
  and corresponds to the iHMP / HMP2 IBDMDB cohort.
- In the local manifest it has:
  - `176` gut samples,
  - rich longitudinal and multi-omics context,
  - clear IBD relevance,
  - existing local metadata rows.
- In the local sample metadata file, `PRJNA398089` currently resolves to:
  - `85` Crohn's disease rows,
  - `46` ulcerative colitis rows,
  - `45` healthy rows.

This makes HMP2 the best immediate Phase 2 candidate.

## Critical Scope Caveat

HMP2 is a very good **IBD validation** cohort, but it is not a perfect
one-to-one ecological match to AGP stool data.

In the current local PRIME manifest, the recorded HMP2 body sites are:

- `Rectum`
- `Colon`
- `Ileum`
- `Cecum`

That means HMP2 should be used first to validate:

- IBD-related DCST structure,
- Crohn vs UC contrasts,
- depth-2 heterogeneity,
- inflammation-associated dominance patterns.

It should **not** be sold as a perfect external validation for every AGP
phenotype such as IBS, reflux, allergies, migraine, or obesity.

## Validation Goals

Phase 2 should answer these questions, in this order:

### Goal 1: structural validation

Do comparable depth-1 and depth-2 DCST patterns emerge in the external cohort?

### Goal 2: disease-focused validation

Do the strongest AGP findings that are relevant to the external cohort
replicate in direction and broad magnitude?

### Goal 3: technical validation

Are the key findings robust to a different cohort, metadata structure, and
sample context?

## Cohort Priority

### Primary: HMP2 / IBDMDB

Use first because it is already discoverable locally and it matches the
highest-value disease family for validation: IBD.

What HMP2 can validate well:

- Crohn-related DCST signals
- UC-related DCST signals
- IBD-enriched Enterobacteriaceae / inflammatory patterns
- depth-2 subtype heterogeneity

What HMP2 cannot validate well:

- non-IBD AGP disease families
- stool-specific population-prevalence claims
- broad citizen-science phenotype screens

### Secondary: MetaHIT

Use if accessible and harmonizable.

Best role:

- stool-based technical/ecological comparison
- metabolic and broader gut-structure comparison

### Tertiary: COMBO or another clinically curated stool cohort

Use if HMP2 proves too mismatched in body site or if a better stool-based IBD
cohort becomes available.

## Parallel Workstreams

## Track A: central manuscript track

This continues regardless of validation progress.

I should own:

- manuscript scaffold refinement
- claim hierarchy
- figure plan
- disease-review integration
- marking which claims are provisional pending validation

Deliverables:

- manuscript draft
- claim inventory updates
- figure/table drafts

## Track B: validation preflight

This should start immediately.

Tasks:

1. confirm exact HMP2 rows in the PRIME abundance and metadata files
2. confirm what body sites are represented in the actual sample-level metadata
3. define the phenotype map:
   - healthy
   - Crohn's disease
   - ulcerative colitis
   - combined IBD
4. define comparability rules with AGP
5. define the replication targets

Deliverable:

- `dev/VALIDATION_PREFLIGHT_HMP2.md`

## Track C: validation execution

This begins once Track B confirms HMP2 is workable.

Tasks:

1. subset `PRJNA398089` from the PRIME abundance table
2. run the same QC and DCST pipeline as closely as possible
3. generate:
   - depth-1 DCST assignments
   - depth-2 DCST assignments
   - size distributions
   - disease association tables
   - concordance summaries against AGP discovery findings
4. produce a compact validation report

Suggested output location:

- `~/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/`

## Track D: alternative-cohort scouting

This is a useful sidecar task and can run in parallel with HMP2 execution.

Tasks:

1. assess whether MetaHIT is locally available or readily ingestible
2. assess whether COMBO or another stool-based IBD cohort is a better ecological
   match
3. compare:
   - body site
   - sequencing region
   - disease labels
   - control availability
   - sample size

Deliverable:

- `dev/VALIDATION_COHORT_COMPARISON.md`

## Replication Targets

Phase 2 should not try to replicate everything. It should target the strongest
and most appropriate claims.

### Tier 1 replication targets

- broad IBD enrichment/depletion patterns at depth 1
- Crohn vs UC contrasts where metadata allows
- depth-2 heterogeneity within major parent DCSTs

### Tier 2 replication targets

- specific high-signal taxa such as Morganella, Bacteroides, Prevotella,
  Faecalibacterium, and Enterobacteriaceae-related types

### Out of scope for HMP2-first validation

- IBS
- acid reflux
- seasonal allergies
- migraine
- obesity
- other non-IBD AGP phenotype bins

## What I Should Own Directly

I should keep these tasks central:

1. choosing the validation target claims
2. deciding what counts as replication
3. writing the interpretation of partial vs failed replication
4. integrating validation results into the manuscript

These are not just data-processing tasks. They shape the scientific argument.

## What Is Good To Delegate

### Subagent candidate 1: validation preflight explorer

Best role:

- inspect local data and metadata for HMP2
- confirm sample counts and body-site caveats
- identify exact phenotype labels and usable columns
- report blockers before any heavy pipeline work starts

Why delegate:

- bounded, read-heavy, easy to run in parallel

### Subagent candidate 2: validation pipeline worker

Best role:

- implement the HMP2 extraction and pipeline run
- write outputs into a dedicated validation directory
- produce machine-readable result tables and a concise report

Why delegate:

- concrete production work with a clear write scope

### Optional subagent candidate 3: alternative cohort scout

Best role:

- compare MetaHIT, COMBO, and any other local candidates
- produce a ranked recommendation memo

Why delegate:

- separate from the immediate HMP2 critical path

## My Recommended Delegation Strategy

Yes, I am proposing delegation, but **not** of the entire Phase 2.

The clean split is:

1. **I keep the scientific control layer.**
   - replication criteria
   - target claims
   - interpretation
   - manuscript integration

2. **Subagents handle bounded sidecar work.**
   - HMP2 preflight
   - HMP2 pipeline execution
   - alternative cohort scouting

## Recommended Sequence

### Step 1: run one preflight subagent now

Purpose:

- verify the HMP2 candidate fully before committing execution work

### Step 2: if HMP2 checks out, run one worker subagent

Purpose:

- execute the validation pipeline and generate outputs

### Step 3: in parallel, optional alternative-cohort scout

Purpose:

- identify whether a stool-based backup cohort would strengthen the paper more
  than HMP2 alone

## Success Criteria

Phase 2 is successful if it produces:

1. a clear statement of what HMP2 validates and what it does not validate,
2. a replication table for the top IBD-related claims,
3. at least one validation figure suitable for the paper,
4. a principled decision on whether HMP2 alone is enough for the first preprint
   or whether a second cohort is worth pursuing.

## Immediate Next Deliverables

1. `dev/VALIDATION_PREFLIGHT_HMP2.md`
2. `dev/VALIDATION_REPLICATION_TARGETS.tsv`
3. validation outputs under
   `~/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/`
4. `dev/VALIDATION_COHORT_COMPARISON.md` if alternative-cohort scouting is run
