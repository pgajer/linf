# Draft Validation Subsection: HMP2 / IBDMDB

Date: 2026-03-25

This file contains manuscript-ready prose for the first external validation
subsection.

## Results Draft

### External IBD-Focused Validation in HMP2 / IBDMDB

To assess whether the AGP findings generalized beyond a single citizen-science
cohort, we performed a first external validation pass in the HMP2 / IBDMDB
cohort (`PRJNA398089`) using the same PRIME-derived SILVA abundance framework.
After phenotype filtering, the local validation set contained 166 samples:
45 healthy controls, 81 Crohn's disease samples, and 40 ulcerative colitis
samples. Two caveats are important. First, in the current local metadata this
cohort is dominated by rectal, ileal, and colonic sites rather than stool,
making it an IBD-focused validation cohort rather than a direct ecological
match to AGP. Second, the cohort is longitudinal, with repeated participant
identifiers, so this first pass is sample-level rather than subject-aware.

Under the same quality-control settings, 217 taxa remained after filtering.
At depth 1, the HMP2 DCST landscape was sparse: with `n0 = 50`, only a
`Bacteroides` dominant type and `RARE_DOMINANT` cleared threshold. The
corresponding depth-1 disease contrasts were not strongly significant after
multiple-testing correction. This lack of strong depth-1 replication should be
interpreted cautiously, because the HMP2 samples differ from AGP both in body
site and repeated-measures structure.

The more informative signal emerged at depth 2. A
`Bacteroides__Faecalibacterium` subtype was significantly depleted in
`IBD_vs_Healthy` (OR = 0.391, q = 0.016), `UC_vs_Healthy` (OR = 0.318,
q = 0.032), and `Crohn_vs_Healthy` (OR = 0.429, q = 0.041). Conversely,
`RARE_DOMINANT` depth-2 structure was enriched in the inflammatory comparisons.
Taken together, these results do not provide a literal taxon-by-taxon
replication of the AGP discovery landscape, but they do support the broader
claim that DCST structure and subtype-level inflammatory contrasts generalize
beyond AGP, and that the most useful external signal may appear at depth 2
rather than depth 1.

## Shorter Version

We next performed an IBD-focused external validation pass in HMP2 / IBDMDB
(`PRJNA398089`). After phenotype filtering, the local validation cohort
contained 166 samples (45 healthy, 81 Crohn's disease, 40 ulcerative colitis).
This cohort provides clinically grounded IBD validation, but it is not a
stool-to-stool match to AGP: in the current local metadata it is dominated by
rectal, ileal, and colonic sites, and it also contains repeated measures.

At depth 1, the HMP2 DCST landscape was sparse and did not yield strong
multiple-testing-significant disease contrasts. At depth 2, however, a
`Bacteroides__Faecalibacterium` subtype was significantly depleted in
`IBD_vs_Healthy` (OR = 0.391, q = 0.016), `UC_vs_Healthy` (OR = 0.318,
q = 0.032), and `Crohn_vs_Healthy` (OR = 0.429, q = 0.041). We therefore view
HMP2 as supporting directional and structural replication of the inflammatory
subtype story, rather than exact depth-1 reproduction of the AGP landscape.

## Discussion Paragraph Draft

The HMP2 validation is informative precisely because it is not a trivial copy
of AGP. It replaces self-reported citizen-science phenotypes with a clinically
grounded IBD cohort, but introduces a different ecological context: mucosal gut
sites rather than stool and longitudinal repeated sampling rather than largely
independent individuals. In that setting, the main depth-1 AGP landscape does
not replicate literally. However, the depletion of a
`Bacteroides__Faecalibacterium` subtype across inflammatory contrasts is
consistent with the broader interpretation that DCST subtype structure captures
meaningful inflammatory heterogeneity. This suggests that external validation
of DCST-based claims may depend less on exact taxon labels at depth 1 than on
whether similar higher-resolution subtype contrasts recur across cohorts.

## Usage Notes

- Use the longer version in the Results section if the paper includes a full
  external-validation subsection.
- Use the shorter version if external validation is presented briefly in the
  first preprint.
- The discussion paragraph can be adapted directly into the Discussion or
  Limitations sections.
