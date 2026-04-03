# Validation Cohort Comparison

Date: 2026-03-25

## Findings-First Recommendation

After HMP2, the best local backup or second validation cohort is:

- `PRJEB84421` (`OFGCD-FI-2025`)

This is the cleanest stool-only inflammatory cohort visible in the current
local PRIME assets, and it complements HMP2 well because HMP2 is IBD-rich but
not stool/feces matched to AGP.

## Important Constraints

1. `MetaHIT` is not currently visible as an actionable local PRIME cohort.
   It appears in literature references, but not as a ready local project in the
   present manifest/report workflow.
2. `COMBO` is in the same situation.
3. There is no single large, clean, stool-only IBD cohort in the current local
   assets that simultaneously matches:
   - HMP2's disease richness, and
   - AGP's stool ecology.

So the practical strategy is:

- use `HMP2` first for disease-focused IBD validation,
- then use a stool-based inflammatory cohort as an ecological complement.

## Ranked Candidates

### 1. `PRJEB84421` (`OFGCD-FI-2025`)

Best use:

- second validation cohort after HMP2
- stool-based inflammatory / Crohn-like directional replication

Why it is strong:

- stool/feces only
- inflammatory focus
- locally visible and extractable
- useful bridge between AGP stool ecology and HMP2 IBD disease focus

Known local characteristics:

- `73` total runs
- metadata and abundance both appear to contain all `73`
- phenotype structure:
  - `24` Crohn
  - `20` Healthy
  - `29` pediatric-onset orofacial granulomatosis

Caveats:

- smaller than HMP2
- pediatric
- the third phenotype is not a clean UC substitute

## 2. `PRJNA1074072` (`ILCI-US-2024`)

Best use:

- targeted stool-based UC vs Healthy directional replication

Why it is useful:

- stool/feces only
- simplest local UC/Healthy backup

Known local characteristics:

- `35` total runs
- `20` UC
- `15` Healthy

Caveats:

- likely underpowered for anything beyond targeted checks

## 3. `PRJNA1216621` (`HEAF-CN-2025`)

Best use:

- non-IBD second validation for the cardiovascular / heart-failure branch

Why it is useful:

- stool/feces only
- locally visible
- reasonable sample size

Known local characteristics:

- `90` total runs
- `60` heart failure
- `30` Healthy

Caveat:

- does not strengthen the main IBD validation story

## Honorable Mention

### `PRJNA1152702` (`VADY-CN-2024`)

Potential use:

- diabetes / metabolic branch

Why not ranked higher:

- mixed body sites
- likely paired stool + vaginal design
- less clean as a first backup cohort

## Not Currently Actionable

### `MetaHIT`

- useful in principle
- not currently surfaced as a ready local PRIME cohort

### `COMBO`

- useful in principle
- not currently surfaced as a ready local PRIME cohort

## What To Do With This

### If the paper stays primarily IBD-centered

Recommended sequence:

1. `HMP2`
2. `PRJEB84421`

### If the paper expands the cardiovascular branch

Recommended sequence:

1. `HMP2`
2. `PRJEB84421`
3. `PRJNA1216621`

## Practical Recommendation

For now:

1. keep `HMP2` as the first validation cohort,
2. treat `PRJEB84421` as the most valuable stool-based follow-up,
3. do not spend time chasing `MetaHIT` or `COMBO` until a clear local access
   path exists.

## Exact Local Sources Used

- `/Users/pgajer/current_projects/gut_microbiome/data/prime_gut_project_manifest_2026-03-24.csv`
- `/Users/pgajer/current_projects/gut_microbiome/data/prime_gut_project_sample_metadata_2026-03-24.csv.gz`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/reports/prime_gut_projects_by_sample_count_2026-03-24.html`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/prime_species/prime_gut_projects_silva_species_absolute_2026-03-24.csv.gz`
