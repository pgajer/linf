# Validation Preflight: HMP2 / IBDMDB / PRJNA398089

Date: 2026-03-25

## Findings-First Summary

1. **Yes, PRJNA398089 can be cleanly identified and extracted from the current
   PRIME assets.**
   - The project is present in the local PRIME manifest.
   - It appears in the PRIME species download log.
   - All `176` metadata runs for `PRJNA398089` were found in the current SILVA
     species abundance table.

2. **The local sample-level metadata is strong enough for an HMP2-first
   validation pass.**
   - Phenotype labels map cleanly to:
     - `Healthy` = `45`
     - `Crohn's disease (CD), IBD` = `85`
     - `Ulcerative colitis (UC), IBD` = `46`
   - This is sufficient for:
     - Healthy vs combined IBD
     - CD vs UC
     - Healthy vs CD
     - Healthy vs UC

3. **The biggest scientific caveat is body site mismatch, not data
   availability.**
   - The current local body-site labels are:
     - `Rectum` = `88`
     - `Ileum` = `76`
     - `Sigmoid Colon` = `5`
     - `Transverse colon` = `2`
     - `Descending (left-sided) colon` = `2`
     - `Ascending (right-sided) colon` = `1`
     - `Cecum` = `1`
     - `Unknown` = `1`
   - This means HMP2 is not a direct stool-to-stool ecological match to AGP.

4. **There is also a repeated-measures caveat.**
   - The project is marked as `Time_series = yes`.
   - Sample names suggest repeated measures:
     - `176` runs
     - about `81` subject-like prefixes from `Sample_Name`
   - So a naive sample-level replication is possible for a first pass, but a
     stronger validation should eventually account for within-subject
     repetition.

## Recommended Decision

**Proceed with HMP2 as the first validation cohort: yes.**

But proceed with the correct scope:

- use HMP2 as an **IBD-focused validation cohort**
- do **not** use it as a universal validator for the full AGP disease screen

## What HMP2 Can Validate Well

1. Whether the DCST framework still yields interpretable depth-1 and depth-2
   structure in an external, clinically grounded gut cohort.
2. Whether IBD-related DCST enrichments and depletions replicate directionally.
3. Whether CD and UC differ in a way that depth-1 or depth-2 DCSTs can detect.
4. Whether depth-2 refinement adds clinically useful heterogeneity beyond
   depth 1.

## What HMP2 Cannot Validate Cleanly

1. IBS, reflux, allergies, migraine, obesity, kidney disease, CVD, and other
   non-IBD AGP phenotype bins.
2. Stool-specific prevalence claims from AGP.
3. A perfect ecological match to citizen-science stool sampling.

## Obvious Blockers

No hard blockers were found for a **first validation pass**.

The main caveats are:

1. **Body-site mismatch**
   - AGP is stool/feces
   - current HMP2 sample labels are mucosal gut sites

2. **Repeated measures**
   - the cohort appears longitudinal
   - validation should eventually become subject-aware

3. **Potential comparability drift**
   - if dominant taxa differ because mucosal and stool ecosystems differ,
     failure to replicate should not be overinterpreted as failure of DCSTs

## Recommended Execution Path

### Pass 1: pragmatic sample-level validation

Use all `176` runs and ask:

- do interpretable DCSTs emerge?
- do broad IBD / CD / UC signals appear?
- does depth-2 help?

This is good enough for a first validation readout.

### Pass 2: stronger subject-aware validation

If Pass 1 is promising:

- derive a more explicit subject identifier
- collapse to one representative sample per subject or use a subject-aware
  model
- re-test the key IBD targets

## Concrete Next Steps

1. Subset `PRJNA398089` from:
   - the SILVA abundance table
   - the project metadata table
2. Create a phenotype map with four labels:
   - Healthy
   - CD
   - UC
   - combined IBD
3. Run the same DCST pipeline as in AGP Phase 1.
4. Generate:
   - depth-1 size distribution
   - depth-2 size distribution
   - IBD / CD / UC association tables
   - one compact validation summary report
5. If the first pass is promising, do the subject-aware rerun.

## Exact File Paths Used

- `/Users/pgajer/current_projects/gut_microbiome/data/prime_gut_project_manifest_2026-03-24.csv`
- `/Users/pgajer/current_projects/gut_microbiome/data/prime_gut_project_sample_metadata_2026-03-24.csv.gz`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/prime_species/prime_gut_projects_silva_species_absolute_2026-03-24.csv.gz`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/prime_species/prime_gut_projects_silva_species_relative_2026-03-24.csv.gz`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/prime_species/download_prime_gut_species_2026-03-24.log`
