# Draft Validation Subsection: PRJEB84421 / OFGCD-FI-2025

Date: 2026-03-25

This file contains manuscript-ready prose for the stool-based inflammatory
complement to HMP2.

## Results Draft

### Stool-Based Inflammatory Complement in PRJEB84421

As a second validation pass, we analyzed `PRJEB84421` (`OFGCD-FI-2025`), a
clean stool/feces cohort that complements the mucosal, longitudinal HMP2
validation. Unlike HMP2, this cohort is a direct stool study and therefore is a
better ecological match to AGP. The local metadata retained 73 samples after
phenotype filtering: 20 healthy controls, 24 Crohn's disease samples, and 29
pediatric-onset orofacial granulomatosis samples. Because the phenotype mix is
not a perfect match to AGP or HMP2, the main goal here is directional and
structural replication, not phenotype-by-phenotype equivalence.

Using a cohort-appropriate threshold of `n0 = 8` at depth 1 and `n0 = 4` at
depth 2, the dominant landscape was led by `Faecalibacterium`,
`Bacteroides`, `Subdoligranulum`, and `RARE_DOMINANT`. Depth-1 contrasts were
modest: the most notable signal was `RARE_DOMINANT` enrichment in
`OFG_vs_Healthy` (OR = 4.93, q = 0.062), along with weaker inflammatory
contrasts driven by `Faecalibacterium` depletion and `Subdoligranulum`
differences. At depth 2, the strongest signal again involved
`RARE_DOMINANT`, which was enriched in the combined inflammatory-vs-healthy
comparison (OR = 4.94, q = 0.084). The `Faecalibacterium__Bacteroides`
subtype showed depletion in inflammatory comparisons, but did not survive
multiple-testing correction as strongly as the HMP2 depth-2 signal.

Taken together, PRJEB84421 is best interpreted as a stool-based complement to
HMP2 that confirms the general DCST framing in a second cohort, but with a more
modest and less cleanly replicated disease-association profile. That is still
useful: it suggests the DCST landscape is not an artifact of a single cohort
type, even though the exact inflammatory subtype signals remain cohort-specific.

## Shorter Version

PRJEB84421 is the stool-based complement to HMP2. It contains 73 stool samples
split across healthy controls, Crohn's disease, and pediatric-onset orofacial
granulomatosis, so it is a better ecological match to AGP than HMP2. Using
`n0 = 8` at depth 1 and `n0 = 4` at depth 2, the cohort yielded a
`Faecalibacterium`/`Bacteroides`/`Subdoligranulum`-centered landscape with a
modest inflammatory signal. The clearest result was `RARE_DOMINANT`
enrichment in `OFG_vs_Healthy` at depth 1 and in combined inflammatory-vs-
healthy at depth 2, but the cohort does not show a strong taxon-by-taxon
replication of the AGP landscape. We therefore use it as an ecological
complement and a robustness check, not as a second strong discovery cohort.

## Discussion Paragraph Draft

PRJEB84421 adds an important complementary validation because it is stool-based
where HMP2 is mucosal. In that sense it helps answer a different question: not
whether AGP-style stool associations replicate in a clinically grounded IBD
cohort, but whether the DCST formalism still yields interpretable structure in a
smaller, stool-only inflammatory study with pediatric-onset disease. The answer
is yes, but only modestly. The resulting landscape is centered on familiar gut
commensals such as `Faecalibacterium`, `Bacteroides`, and `Subdoligranulum`,
and the dominant signal is `RARE_DOMINANT` rather than a single dramatic
pathobiont. That is still scientifically useful because it reinforces the idea
that DCSTs capture structured dominance patterns across cohorts, even when the
exact disease associations differ.

## Usage Notes

- Use this subsection alongside the HMP2 validation subsection.
- The cohort should be described as a stool-based inflammatory complement,
  not as an independent high-power replication of AGP.
- The OFG phenotype should be treated cautiously and not equated with UC.
