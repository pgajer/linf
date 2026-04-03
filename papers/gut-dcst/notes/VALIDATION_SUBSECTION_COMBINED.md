# Draft Validation Subsection: Two-Cohort External Validation

Date: 2026-03-25

This file combines the HMP2 and PRJEB84421 validation passes into one compact,
paper-ready validation unit.

## Results Draft

### External Validation Across Two Complementary Cohorts

To assess whether the AGP findings generalized beyond a single citizen-science
cohort, we carried out two complementary external validation passes using the
same PRIME-derived SILVA abundance framework. The first used HMP2 / IBDMDB
(`PRJNA398089`), which provides clinically grounded IBD phenotypes but is, in
the current local metadata, dominated by rectal, ileal, and colonic sites and
contains repeated measures. After phenotype filtering, the local HMP2 set
contained 166 samples: 45 healthy controls, 81 Crohn's disease samples, and 40
ulcerative colitis samples. The second used `PRJEB84421`
(`OFGCD-FI-2025`), a stool-only inflammatory cohort that is a much better
ecological match to AGP but is smaller and pediatric, with 73 retained samples
split across 20 healthy controls, 24 Crohn's disease samples, and 29
pediatric-onset orofacial granulomatosis samples. We therefore treated HMP2 as
the disease-focused IBD validation cohort and PRJEB84421 as the stool-based
ecological complement.

The two validation cohorts contributed different kinds of evidence. In HMP2,
the depth-1 DCST landscape was sparse: with `n0 = 50`, only a `Bacteroides`
dominant type and `RARE_DOMINANT` cleared threshold, and the corresponding
depth-1 disease contrasts were not strongly significant after multiple-testing
correction. The more informative signal emerged at depth 2, where a
`Bacteroides__Faecalibacterium` subtype was significantly depleted in
`IBD_vs_Healthy` (OR = 0.391, q = 0.016), `UC_vs_Healthy` (OR = 0.318,
q = 0.032), and `Crohn_vs_Healthy` (OR = 0.429, q = 0.041). This supports
directional and structural replication of the inflammatory subtype story, even
though the AGP depth-1 landscape is not reproduced literally.

In PRJEB84421, we used cohort-appropriate thresholds of `n0 = 8` at depth 1
and `n0 = 4` at depth 2 to preserve interpretable structure in a smaller
dataset. The resulting stool-based landscape was led by `Faecalibacterium`,
`Bacteroides`, `Subdoligranulum`, and `RARE_DOMINANT`. The strongest signals
were more modest and cohort-specific than in HMP2: at depth 1,
`RARE_DOMINANT` was enriched in `OFG_vs_Healthy` (OR = 4.93, q = 0.062), and
at depth 2 it was enriched in the combined inflammatory-vs-healthy comparison
(OR = 4.94, q = 0.084). A `Faecalibacterium__Bacteroides` subtype showed
depletion in inflammatory contrasts, but not as strongly as the analogous HMP2
depth-2 signal.

Taken together, these external analyses support a more nuanced validation claim
than simple taxon-by-taxon replication. HMP2 shows that DCST subtype structure
captures clinically relevant inflammatory heterogeneity in a longitudinal IBD
cohort, while PRJEB84421 shows that the DCST framework remains interpretable in
a direct stool-based inflammatory cohort. We therefore view the validation as
supporting **structural and directional generalization** of the DCST framework
rather than exact reproduction of the AGP discovery landscape.

## Shorter Version

We next performed external validation in two complementary cohorts. HMP2 /
IBDMDB (`PRJNA398089`) provided clinically grounded IBD validation but was not
a stool-to-stool match to AGP, because in the current local metadata it is
dominated by rectal, ileal, and colonic sites and includes repeated measures.
After phenotype filtering, it contained 166 samples (45 healthy, 81 Crohn's
disease, 40 ulcerative colitis). In this cohort, the main validation signal
appeared at depth 2 rather than depth 1: a `Bacteroides__Faecalibacterium`
subtype was depleted in `IBD_vs_Healthy` (OR = 0.391, q = 0.016),
`UC_vs_Healthy` (OR = 0.318, q = 0.032), and `Crohn_vs_Healthy`
(OR = 0.429, q = 0.041).

As a stool-based ecological complement, we analyzed `PRJEB84421`, a 73-sample
cohort containing healthy controls, Crohn's disease, and pediatric-onset
orofacial granulomatosis. This second cohort yielded a clear DCST landscape
centered on `Faecalibacterium`, `Bacteroides`, `Subdoligranulum`, and
`RARE_DOMINANT`, but the disease signals were more modest and cohort-specific,
with the strongest contrasts involving `RARE_DOMINANT`. Together, these two
cohorts support directional and structural generalization of the DCST
framework, even though exact AGP-style depth-1 replication is limited.

## Discussion Paragraph Draft

The two validation cohorts are useful precisely because they stress different
parts of the AGP story. HMP2 replaces self-reported citizen-science phenotypes
with a clinically grounded IBD cohort, but at the cost of a different
ecological context and repeated-measures structure. PRJEB84421 restores the
stool ecology but introduces a smaller pediatric inflammatory cohort with a
nonstandard phenotype mix. In this paired setting, the most defensible
validation claim is not exact duplication of AGP depth-1 taxa, but rather that
DCSTs continue to produce interpretable dominant-community structure across
external cohorts, and that subtype-level inflammatory contrasts, especially
those involving `Faecalibacterium`-related structure and rare-dominant tails,
recur in informative ways.

## Usage Notes

- Use this file as the manuscript-facing validation subsection.
- Keep the cohort-specific files as supporting notes:
  - `VALIDATION_SUBSECTION_HMP2.md`
  - `VALIDATION_SUBSECTION_PRJEB84421.md`
- The key phrase to preserve is:
  `structural and directional generalization rather than exact replication`.
