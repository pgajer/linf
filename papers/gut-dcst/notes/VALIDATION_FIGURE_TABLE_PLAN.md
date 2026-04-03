# Validation Figure and Table Plan

Date: 2026-03-25

## Goal

Add one compact external-validation unit to the paper without letting the
validation section overwhelm the main AGP story.

The validation visuals should communicate:

1. why two external cohorts were needed,
2. what each cohort contributes,
3. why neither cohort is a perfect AGP surrogate,
4. what actually generalized across cohorts.

## Recommended Validation Figure

### Figure V1: Two-cohort external validation overview

Recommended panels:

- panel A: HMP2 depth-1 DCST size distribution
- panel B: HMP2 depth-2 DCST size distribution
- panel C: PRJEB84421 depth-1 DCST size distribution
- panel D: compact directional replication panel summarizing the strongest
  HMP2 and PRJEB84421 findings

Existing inputs:

- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/hmp2_depth1_sizes.png`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/hmp2_depth2_sizes.png`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/hmp2_depth2_results.csv`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/hmp2/hmp2_run_info.csv`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/prjeb84421/prjeb84421_depth1_sizes.png`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/prjeb84421/prjeb84421_depth2_results.csv`
- `/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_validation/prjeb84421/prjeb84421_run_info.csv`

Recommended message:

- HMP2 supplies clinically grounded IBD validation but mostly at depth 2
- PRJEB84421 supplies stool-based ecological validation with more modest,
  cohort-specific signals
- together they support structural / directional generalization

Preferred design for panel D:

- a small forest or dot plot using:
  - HMP2 `Bacteroides__Faecalibacterium` for `IBD_vs_Healthy`
  - HMP2 `Bacteroides__Faecalibacterium` for `Crohn_vs_Healthy`
  - HMP2 `Bacteroides__Faecalibacterium` for `UC_vs_Healthy`
  - PRJEB84421 `RARE_DOMINANT` for `Inflammatory_vs_Healthy`
  - PRJEB84421 `RARE_DOMINANT` for `OFG_vs_Healthy`

## Recommended Validation Table

### Table V1: External validation cohort summary

Columns:

- cohort
- samples analyzed
- phenotype split
- ecological caveat
- strongest depth-1 signal
- strongest depth-2 signal
- interpretation

Proposed rows:

- HMP2 / IBDMDB
- PRJEB84421 / OFGCD-FI-2025

Suggested values:

- HMP2 / IBDMDB:
  - samples analyzed: `166`
  - phenotype split: `45 Healthy / 81 Crohn / 40 UC`
  - ecological caveat: `mucosal gut sites; repeated measures`
  - strongest depth-1 signal: `no q < 0.05 depth-1 replication`
  - strongest depth-2 signal:
    `Bacteroides__Faecalibacterium depleted in IBD_vs_Healthy, Crohn_vs_Healthy,
    and UC_vs_Healthy`
  - interpretation:
    `directional / structural replication at depth 2`

- PRJEB84421 / OFGCD-FI-2025:
  - samples analyzed: `73`
  - phenotype split: `20 Healthy / 24 Crohn / 29 OFG`
  - ecological caveat: `stool-based but pediatric and phenotype-mixed`
  - strongest depth-1 signal:
    `RARE_DOMINANT enriched in OFG_vs_Healthy`
  - strongest depth-2 signal:
    `RARE_DOMINANT enriched in inflammatory_vs_healthy`
  - interpretation:
    `stool-based structural complement with modest cohort-specific signal`

## Optional Comparison Table

### Table V2: Discovery-to-validation target map

Purpose:

- show which Phase 1 claims were testable in HMP2, PRJEB84421, or both

Suggested columns:

- discovery target
- testable in HMP2? yes/no
- testable in PRJEB84421? yes/no
- reason
- validation outcome

This table is especially useful if we want to be explicit that HMP2 cannot
validate IBS, reflux, obesity, or other non-IBD AGP branches, and that
PRJEB84421 is not a UC or full-IBD substitute.

## Recommended Placement In Manuscript

- put Figure V1 in the Results subsection on external validation
- put Table V1 either in the main text or supplement depending on space
- use Table V2 in the supplement if the main text becomes too heavy

## Design Principles

1. Keep the validation section compact.
2. Make the HMP2 body-site mismatch visible in text or caption.
3. Make the PRJEB84421 pediatric / OFG caveat equally visible.
4. Emphasize directional replication over exact duplication.
5. Do not let the absence of strong depth-1 HMP2 replication read as method
   failure without context.
