---
title: "Dominant Community State Types Reveal Disease-Associated Dominance Patterns in the Human Gut Microbiome"
author:
  - |
    \parbox{0.9\textwidth}{\centering Pawel Gajer$^{1}$\thanks{Research reported in this publication was supported by grant number INV-048956 from the Gates Foundation.}\thanks{Corresponding author: pgajer@gmail.com} and Jacques Ravel$^{1}$\\[0.35em] $^{1}$University of Maryland School of Medicine}
date: ""
documentclass: article
fontsize: 11pt
geometry: margin=1in
linestretch: 1.15
bibliography: "/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/references.bib"
link-citations: true
header-includes:
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{array}
  - \usepackage{tabularx}
  - \usepackage{float}
  - \usepackage{titling}
  - \usepackage{caption}
  - \captionsetup{font=small,labelfont=bf}
  - \newcolumntype{Y}{>{\raggedright\arraybackslash}X}
  - \setlength{\droptitle}{-2.6em}
  - \pretitle{\begin{center}\Large}
  - \posttitle{\par\end{center}\vspace{0.55em}}
  - \preauthor{\begin{center}\normalsize}
  - \postauthor{\par\end{center}\vspace{-0.5em}}
  - \predate{}
  - \postdate{}
---

# Abstract

Gut microbiome community typing is often built around clustering, making the
resulting labels sensitive to modeling choices and harder to interpret at the
level of a single sample. We study a deterministic alternative based on
dominant-community state types (DCSTs), in which each sample is assigned by
within-sample taxon ranking rather than by cluster membership.
Using the American Gut Project stool cohort, we analyzed 34,711 samples
available in the SILVA species-level table and retained 30,290 after quality
control, yielding 274 taxa for downstream analysis. At depth 1, the cohort
contained 40 common DCSTs plus a rare-dominant tail; at depth 2, subtype
structure revealed additional heterogeneity within broad dominant groups. Among
452 adjusted depth-1 association tests, 28 remained significant after
Benjamini-Hochberg correction, and 23 remained significant after a
contamination-aware sensitivity analysis. Depth-2 modeling contributed 23
additional adjusted-significant associations, while a comparison against GG2
showed 91.6\% exact dominant-genus agreement across overlapping samples.
Disease-wise analysis indicated that the clearest signals concentrate in
inflammatory and gastrointestinal phenotypes, especially IBD-related structure.
External validation in two complementary cohorts, HMP2 and PRJEB84421,
supported structural and directional generalization rather than exact
AGP-style replication. Together, these results position DCSTs as a practical,
interpretable, and hierarchical framework for large-scale gut microbiome
analysis.

# 1. Introduction

Microbiome community typing is useful when it compresses a high-dimensional
abundance profile into a smaller number of biologically interpretable states.
In practice, however, many typing schemes are built around clustering. That can
work well, but it also ties the final labels to choices about transformations,
distance metrics, initialization, model class, and the number of clusters.
Those choices are not merely technical. They affect whether a reader can
understand what a given sample label means and whether two analyses can be
compared without a hidden dependence on tuning decisions.

This paper studies a simpler route. A dominant-community state type (DCST) is
defined by the within-sample ordering of taxa. At depth 1, the dominant taxon
determines the top-level state. At depth 2, the second-ranked taxon refines
that state into an ordered subtype. Because the construction depends only on
rank order, it is unchanged by sample-wide rescaling and does not require a
particular $L_p$ normalization. The resulting construction is deterministic,
hierarchical, and local to the sample: the label is read directly from the
ordered dominant coordinates rather than inferred from a fitted global
partition.

The gut is a particularly good stress test for this idea. In the vaginal
microbiome, strong dominance can be common and even canonical. In the gut,
strong dominance is less universal and often appears in dysbiotic or otherwise
unusual settings. That makes the gut harder to type cleanly, but also more
interesting. Literature on IBS, IBD, autoimmune disease, and cardiometabolic
conditions suggests that shifts in dominant taxa or dominant taxon pairs can
mark altered ecological states [@Tap2017_IBS; @Pittayanon2019_IBS;
@LloydPrice2019_IBD; @Gevers2014_CD; @Scher2013_Prevotella_RA; @Vatanen2018_T1D;
@Wang2011_TMAO; @Tang2013_TMAO].

Our aim is not to argue that dominant taxa alone explain disease biology. The
claim is narrower: recurring dominance patterns can be identified
deterministically across a cohort of this size without fitting clusters; the
resulting DCST labels can then be analyzed phenotype by phenotype rather than
only through one pooled disease-versus-not-disease contrast; and several of the
strongest DCST-phenotype associations remain visible after covariate
adjustment, contaminant-aware sensitivity analysis, taxonomy change, and
external validation.

> **What is a DCST?**
>
> A DCST is defined by the ordered dominant taxa in a sample. At depth 1, only
> the largest component matters. At depth 2, the second-largest component
> refines the parent class. Because only the within-sample ordering matters, the
> same DCST is obtained from read counts, proportions, or any sample-wide
> $L_p$ normalization. How low-frequency dominant patterns are handled is a
> separate modeling choice. In the `rare` policy used here, samples that do not
> meet the prevalence threshold for a named dominant type are grouped into
> `RARE_DOMINANT`.

![Conceptual schematic of the DCST construction. A sample abundance profile induces a within-sample taxon ranking, which yields a depth-1 state from the dominant taxon and a depth-2 state from the ordered dominant pair. Under the `rare` low-frequency policy used in this paper, dominant patterns below the prevalence threshold are grouped into `RARE_DOMINANT`; the software also supports an `absorb` policy that merges such samples into larger named states.](/Users/pgajer/current_projects/linf/dev/FIGURE_1_dcst_conceptual_schematic.png){ width=92% }

# 2. Methods

## 2.1 Cohort and preprocessing

The discovery analysis used AGP gut samples from the PRIME-derived SILVA
species-level abundance table [@McDonald2018_AGP; @Quast2013_SILVA]. The
starting table contained 34,711 samples, of which 30,290 passed
quality-control filters and were retained for analysis. The retained matrix
contained 274 taxa. Age was available for 27,372 samples (90.4\%), sex for
25,598 (84.5\%), and BMI for 26,541 (87.6\%), producing 22,632 complete cases
for adjusted modeling.

## 2.2 DCST construction

A depth-1 DCST was defined by the dominant taxon. A depth-2 DCST was defined by
the ordered pair consisting of the dominant and second-ranked taxa. This
construction depends only on within-sample rank order, so the resulting label
is unchanged by moving between raw counts, relative abundances, or any
sample-wide $L_p$ normalization. To keep the resulting state space
interpretable, low-frequency dominant patterns were handled by a separate
policy choice. In the software this corresponds to `low.freq.policy`. Under
`low.freq.policy = "rare"` (used here), only dominant patterns meeting a
minimum prevalence threshold $n_0$ were promoted to named DCSTs, and the rest
were grouped into `RARE_DOMINANT`. This preserves the purity of the named DCSTs
because rare cells are not absorbed into larger ones. The alternative
`low.freq.policy = "absorb"` merges low-frequency cells into larger named
states, which can yield a coarser partition but disturbs the original dominant
structure of those named DCSTs.

## 2.3 Association analysis

We first performed univariate screening with Fisher's exact test. We then fit
adjusted logistic models for disease indicators using the DCST label together
with available covariates. For a mathematically oriented reader, the odds ratio
may be read as follows: an odds ratio above 1 means the DCST is enriched in the
given phenotype relative to the reference group, while an odds ratio below 1
means depletion. Multiple testing was handled by the Benjamini-Hochberg
procedure [@BenjaminiHochberg1995], and q-values below 0.05 were treated as
statistically notable rather than as proof of mechanism.

## 2.4 Sensitivity analysis and taxonomy comparison

Because some nominal gut signals can be driven by oral- or skin-associated
contaminants, we repeated the adjusted analysis after excluding samples
dominated by a prespecified contaminant set. We also ran a parallel
dominant-type analysis under GG2 taxonomy [@McDonald2024_GG2] and compared
dominant-genus assignments between SILVA and GG2 on the overlapping sample set.

## 2.5 External validation

We used two complementary external cohorts. HMP2 / IBDMDB
(`PRJNA398089`) provided clinically grounded inflammatory bowel disease
phenotypes but, in the currently available local metadata, was dominated by
mucosal gut sites and included repeated measures. PRJEB84421 provided a
stool-based inflammatory complement, but was smaller, pediatric, and included
orofacial granulomatosis alongside Crohn's disease and healthy controls. The
same DCST pipeline was applied to both cohorts with cohort-appropriate
prevalence thresholds.

# 3. Cohort Overview

## 3.1 Shared AGP overview

The Phase 1 discovery analysis used 30,290 quality-controlled AGP gut samples.
At depth 1 with $n_0 = 50$, the cohort yielded 40 named DCSTs plus
`RARE_DOMINANT` (1,162 samples). The largest types were `Bacteroides` (10,358
samples), `Escherichia-Shigella` (6,857), `Prevotella_9` (1,870),
`Enterobacteriaceae` (1,564), and `Faecalibacterium` (1,320). In other words,
the gut landscape is neither one universal dominant state nor an
uninterpretable spray of rare labels.

Across the full screen, 59 of 492 univariate Fisher tests were significant
after correction. The covariate-adjusted depth-1 pass yielded 28 significant
associations across 452 tests, and 23 of those remained significant in the
contamination-aware sensitivity analysis. Depth 2 added 23 adjusted-significant
subtype associations across 1,340 tests. This common overview is useful because
the disease-wise sections below should be read against this shared statistical
background rather than as isolated anecdotes. The twelve most common depth-1
states already cover 87.6\% of the retained cohort, so the dominant-state
picture is concentrated rather than diffusely fragmented (Table 1).

**Table 1. Top depth-1 DCSTs in the full AGP cohort.**

\begingroup\small

| Rank | DCST | n | Share (%) | Cumulative (%) |
| --- | --- | --- | --- | --- |
| 1 | Bacteroides | 10358 | 34.2 | 34.2 |
| 2 | Escherichia-Shigella | 6857 | 22.6 | 56.8 |
| 3 | Prevotella_9 | 1870 | 6.2 | 63.0 |
| 4 | Enterobacteriaceae | 1564 | 5.2 | 68.2 |
| 5 | Faecalibacterium | 1320 | 4.4 | 72.5 |
| 6 | Rare dominant | 1162 | 3.8 | 76.4 |
| 7 | Streptococcus | 827 | 2.7 | 79.1 |
| 8 | Pseudomonas | 664 | 2.2 | 81.3 |
| 9 | Staphylococcus | 583 | 1.9 | 83.2 |
| 10 | Prevotella_7 | 476 | 1.6 | 84.8 |
| 11 | Unassigned | 444 | 1.5 | 86.2 |
| 12 | Lachnospiraceae | 416 | 1.4 | 87.6 |

\endgroup

![Full-cohort AGP landscape. Panel A shows the leading depth-1 DCSTs after truncation at n0 = 50, while Panel B shows the empirical dominance strength, defined as the relative abundance of the dominant species in each sample.](/Users/pgajer/current_projects/linf/dev/FIGURE_2_full_cohort_landscape.png){ width=92% }

![Adjusted depth-1 association overview. Panel A shows all taxa that survive adjustment at least once across the disease screen, with colored cells annotated by adjusted odds ratio and blank cells indicating no surviving adjusted signal. Panel B shows the strongest adjusted associations with 95% confidence intervals.](/Users/pgajer/current_projects/linf/dev/FIGURE_3_adjusted_association_overview.png){ width=98% }

Table 2 collects the single strongest adjusted depth-1 signal for each
phenotype with at least one surviving adjusted association. This is not meant
to replace the full association screen; it is meant to keep the shared overview
legible before the disease-wise sections unpack the same results in context.

\begin{table}[H]
\footnotesize
\caption{Headline adjusted depth-1 associations by phenotype. Here $C$ indicates whether the same association remains significant with $q < 0.05$ in the contamination-aware clean analysis.}
\centering
\begin{tabularx}{\textwidth}{@{} l l Y c >{\raggedright\arraybackslash}p{3.2cm} c c @{}}
\toprule
Phenotype & Direction & DCST & $n$ & Adjusted OR (95\% CI) & adj $q$ & $C$ \\
\midrule
IBS & Enriched & Escherichia-Shigella & 6857 & 1.23 (1.12-1.35) & 6.8e-04 & Yes \\
IBD & Enriched & Morganella & 53 & 16.94 (8.88-32.31) & 4.0e-15 & Yes \\
Autoimmune & Enriched & Bacteroides & 10358 & 1.39 (1.28-1.52) & 2.1e-11 & Yes \\
Acid reflux & Enriched & Bacteroides & 10358 & 1.28 (1.17-1.40) & 6.9e-06 & Yes \\
Seasonal allergies & Depleted & Prevotella\_9 & 1870 & 0.67 (0.59-0.76) & 3.3e-08 & Yes \\
CDI & Enriched & Escherichia-Shigella & 6857 & 1.48 (1.19-1.84) & 0.0112 & Yes \\
Lung disease & Enriched & Staphylococcus & 583 & 1.55 (1.17-2.05) & 0.0386 & No \\
Migraine & Enriched & Bacteroides & 10358 & 1.16 (1.05-1.26) & 0.034 & No \\
Kidney disease & Enriched & Lactobacillus & 117 & 5.85 (2.26-15.14) & 0.00779 & Yes \\
\bottomrule
\end{tabularx}
\end{table}

Looking only at the counts, the sensitivity analysis reduces the number of
adjusted-significant depth-1 associations from 28 to 23. A more informative
view is to compare which signals are lost, which are retained, and which become
clearer once likely oral- or skin-associated dominants are excluded. Twenty
adjusted-significant associations are retained, eight are lost, and three
become significant only in the clean analysis. The lost signals are concentrated
in `Staphylococcus`- and `Streptococcus`-linked contrasts, while the strongest
IBD, allergy, and several `Bacteroides` / `Prevotella_9` associations remain in
place. Table 3 shows how this transition looks phenotype by phenotype rather
than only in aggregate.

![Contamination-aware robustness summary. The clean analysis removes several oral-associated contrasts but preserves the strongest inflammatory and allergy-linked signals.](/Users/pgajer/current_projects/linf/dev/FIGURE_4_contamination_robustness.png){ width=98% }

\begin{table}[H]
\footnotesize
\caption{Disease-wise comparison of univariate, adjusted, and clean results. Here $n_F$, $n_A$, and $n_C$ denote the numbers of depth-1 associations with Fisher $q < 0.05$, adjusted $q < 0.05$, and clean-analysis $q < 0.05$, respectively.}
\centering
\begin{tabularx}{\textwidth}{l c c c Y Y}
\toprule
Phenotype & $n_F$ & $n_A$ & $n_C$ & Representative retained & Representative Fisher-only \\
\midrule
IBS & 7 & 4 & 3 & Escherichia-Shigella (1.29 -> 1.23) & Corynebacterium (0.29 -> 0.41) \\
IBD & 9 & 7 & 8 & Morganella (18.05 -> 16.94) & Prevotella\_7 (0.10 -> 0.44) \\
Autoimmune & 9 & 6 & 4 & Bacteroides (1.34 -> 1.39) & Escherichia-Shigella (1.19 -> 0.99) \\
Acid reflux & 10 & 3 & 1 & Bacteroides (1.29 -> 1.28) & Prevotella\_7 (0.12 -> 0.55) \\
CVD & 6 & 0 & 0 & -- & Pasteurellaceae (8.77 -> 1.74) \\
Obesity & 3 & 0 & 0 & -- & Prevotella\_7 (2.75 -> not estimable) \\
Seasonal allergies & 6 & 3 & 3 & Prevotella\_9 (0.75 -> 0.67) & Prevotella\_7 (0.22 -> 1.11) \\
Diabetes & 2 & 0 & 0 & -- & Streptococcus (0.28 -> 0.44) \\
CDI & 2 & 1 & 2 & Escherichia-Shigella (1.45 -> 1.48) & Prevotella\_9 (0.47 -> 0.35) \\
Lung disease & 3 & 1 & 0 & Staphylococcus (1.26 -> 1.55) & Prevotella\_7 (0.23 -> 0.89) \\
Migraine & 2 & 1 & 0 & Bacteroides (1.19 -> 1.16) & Corynebacteriaceae (0.20 -> 0.30) \\
Kidney disease & 0 & 2 & 2 & Lactobacillus (3.50 -> 5.85) & -- \\
\bottomrule
\end{tabularx}
\end{table}

Depth 2 matters most when a broad parent state hides materially different
substructures. The depth-2 case-study figure makes that concrete for the depth-1 `Bacteroides`
state in IBD. The ten most common depth-2 children cover 85\% of all
`Bacteroides`-dominant samples, yet their IBD odds ratios range from strong
enrichment (`Lachnospiraceae`, `Enterobacteriaceae`) to attenuation or partial
reversal (`Faecalibacterium`, `Agathobacter`, `Parabacteroides`). The depth-1
`Bacteroides` signal should therefore be read as an average over heterogeneous
subtypes rather than as a single ecological state.

![Depth-2 refinement case study for the broad Bacteroides parent in IBD. Panel A shows the most common depth-2 children inside the parent state. Panel B shows that those common children have materially different adjusted odds ratios, with the dashed line marking the coarser depth-1 parent effect.](/Users/pgajer/current_projects/linf/dev/FIGURE_6_depth2_refinement_case_study.png){ width=95% }

## 3.2 Taxonomy concordance: SILVA versus GreenGenes2

An independent dominant-type pass under GG2 identified 59 named depth-1 DCSTs
plus `RARE_DOMINANT`. On the overlapping 30,290 samples, the exact
dominant-genus agreement between SILVA and GG2 was 27,739 / 30,290 (91.6\%).
This does not mean that the taxonomies are interchangeable. It does mean that
the main dominant-typing picture is fairly stable across reference databases,
even though some clinically relevant groups, especially
`Escherichia-Shigella`, expose genuine naming and resolution differences. The
new concordance figure makes an important distinction explicit: harmonized
genus-level agreement is high, but exact short-label equality is much lower
because GG2 splits some large SILVA states, especially `Bacteroides`, into
multiple named groups and leaves `Escherichia-Shigella` largely unresolved.
Table 4A summarizes the headline concordance metrics, while Table 4B gives a
few representative remappings that matter biologically and interpretively.

![Cross-taxonomy concordance summary. Harmonized agreement between SILVA and GG2 is high, while exact short-label equality is much lower because a few large dominant states are reclassified rather than contradicted.](/Users/pgajer/current_projects/linf/dev/FIGURE_6_taxonomy_concordance.png){ width=95% }

**Table 4A. Agreement metrics for SILVA versus GG2 depth-1 DCST assignments.**

\begingroup\small

| Metric | Value | Interpretation |
| --- | --- | --- |
| Overlapping samples | 30,290 | Common sample set available for direct SILVA-GG2 comparison. |
| Harmonized dominant-genus agreement | 27,739 / 30,290 (91.6%) | Broad depth-1 landscape is stable after genus-level harmonization. |
| Exact short-label equality | 4,566 / 30,290 (15.1%) | Most disagreements are relabelings or taxonomic splits rather than ecological contradictions. |
| Named depth-1 states | SILVA: 40 + rare; GG2: 59 + rare | GG2 resolves some large SILVA states into more named groups. |

\endgroup

**Table 4B. Representative SILVA-to-GG2 remappings.**

\begingroup\small

| SILVA state | Main GG2 destination(s) | Row share | Reading |
| --- | --- | --- | --- |
| Akkermansia | Akkermansia | 99.1% | Near one-to-one mapping. |
| Faecalibacterium | Faecalibacterium | 81.7% | Mostly preserved with modest spillover. |
| Prevotella_9 | Prevotella | 90.5% | High genus-level concordance despite label change. |
| Bacteroides | Phocaeicola_A; Bacteroides_H_857956; unresolved | 46.6%; 17.4%; 8.4% | Large SILVA state is split across multiple GG2 labels. |
| Escherichia-Shigella | Unresolved | 99.7% | GG2 leaves this dominant state largely unresolved. |
| Enterobacteriaceae | Enterobacteriaceae_A_725029; Enterobacterales_737866 | 92.6%; 4.3% | Mostly preserved at family/order level with minor split. |

\endgroup

\clearpage

# 4. Disease-Wise Results and Interpretation

The subsections in this section follow a common pattern. Each begins with a
brief literature frame, then summarizes the AGP signal across SILVA, GG2, depth
2, and sensitivity analysis, and finally interprets the result at a level that
matches the strength of the evidence. This repeated structure is intentional: it
lets the reader compare diseases without losing track of which claims are
strong, which are tentative, and which are mainly useful for context. Several
AGP phenotypes, especially autoimmune disease, seasonal allergies, and lung
disease, are best understood as discovery-screen labels rather than harmonized
clinical diagnoses, so their discussion leans more on ecological interpretation
than on disease-specific inference.

## 4.1 Irritable bowel syndrome (IBS)

### Introduction

IBS is a functional gastrointestinal disorder in which the gut microbiota has
been implicated through immune signaling, short-chain fatty acid production, and
the microbiota-gut-brain axis [@Pittayanon2019_IBS; @Dinan2015_IBS;
@Cryan2012_gutbrain]. Because many IBS signatures are subtle rather than
catastrophically dysbiotic, it is a useful test of whether dominant-state labels
retain signal in a noisy population cohort.

### Results

In SILVA, the strongest adjusted depth-1 IBS association was enrichment of
`Escherichia-Shigella` (adjusted OR 1.23, q = 6.8e-04), accompanied by
protective associations for `Faecalibacterium` (adjusted OR 0.64),
`Staphylococcus` (adjusted OR 0.46), and `Streptococcus`
(adjusted OR 0.61). GG2 showed a qualitatively similar picture, with an
enriched unresolved enterobacterial signal and depleted `Staphylococcus` and
`Streptococcus` states. At depth 2, additional IBS-related subtype signal
appeared in `Escherichia-Shigella`- and `Bacteroides`-centered subtypes.
Sensitivity analysis retained `Faecalibacterium`, `Escherichia-Shigella`, and
`Prevotella_9`, suggesting that the main IBS pattern is not driven entirely by
putative oral contaminants.

### Discussion

These IBS findings are consistent with literature linking Enterobacteriaceae
expansion and loss of anti-inflammatory commensals to symptom burden
[@Pittayanon2019_IBS; @Tap2017_IBS]. The `Faecalibacterium` signal is
particularly interpretable for a mechanistic reader, because it points toward
loss of butyrate-associated ecological stability rather than toward one
IBS-specific pathogen. The oral-associated taxa should still be treated
cautiously, and the overall effect sizes are modest compared with IBD.

## 4.2 Inflammatory bowel disease (IBD)

### Introduction

IBD provides the most clinically grounded inflammatory phenotype in the AGP
screen. Landmark studies report reduced diversity, expansion of Proteobacteria,
and depletion of butyrate-associated commensals, especially
*Faecalibacterium prausnitzii* [@Frank2007_IBD; @Gevers2014_CD;
@LloydPrice2019_IBD; @Sokol2008_Fprausnitzii; @Manichanh2006_CD].

### Results

IBD produced the strongest adjusted depth-1 signal in the entire Phase 1 scan:
`Morganella` was highly enriched (adjusted OR 16.94, q = 4.0e-15). Additional
adjusted SILVA signals included enrichment of `Bacteroides`,
`Lachnospiraceae`, `Bifidobacterium`, `Enterobacterales`, and `Proteus`, along
with depletion of `Prevotella_9`. GG2 showed concordant inflammatory structure,
again highlighting `Morganella` and related enterobacterial signals. Depth 2
was even more striking: a `Bacteroides` subtype with a `Lachnospiraceae`
second-ranked component reached adjusted OR 3.54 with q = 1.6e-37, and multiple
IBD-associated subtypes appeared within the broad `Bacteroides` parent class.
Eight IBD associations survived sensitivity analysis.

### Discussion

IBD is the clearest case here where depth 2 matters both scientifically and
mathematically. The established IBD literature emphasizes depletion of obligate
anaerobes and butyrate-associated commensals, especially
*Faecalibacterium prausnitzii*, together with expansion of inflammatory
facultative taxa [@Frank2007_IBD; @Gevers2014_CD; @LloydPrice2019_IBD;
@Sokol2008_Fprausnitzii; @Sartor2008_IBD]. The present DCST results fit that
picture in a compressed rank-based language: a rare `Morganella` state marks
severe enterobacterial dysbiosis, whereas the broad `Bacteroides` parent
resolves into subtypes with materially different odds ratios. That matters for
interpretation because the enriched depth-1 `Bacteroides` label is not the
opposite of the classic `Faecalibacterium`-depletion story; rather, it partly
reflects loss of protective second-rank partners inside Bacteroides-dominant
samples. Section 5 strengthens this reading by showing directional depletion of
a `Bacteroides__Faecalibacterium` subtype in the external IBD cohorts.

## 4.3 Autoimmune diseases

### Introduction

The AGP autoimmune phenotype is an umbrella label rather than a single disease,
so it should be interpreted more cautiously than IBD. Even so, autoimmune
disease is an important category because the literature repeatedly links loss of
immune-regulatory commensals and expansion of pathobionts to altered Treg/Th17
balance [@Scher2013_Prevotella_RA; @Atarashi2013_Clostridia_Treg;
@DeLuca2019_autoimmune; @Vatanen2018_T1D].

### Results

In SILVA, autoimmune disease showed adjusted enrichment of `Bacteroides`
(adjusted OR 1.39, q = 2.1e-11) and `Morganella`
(adjusted OR 6.88, q = 1.9e-07), together with depletion of `Prevotella_9`,
`Streptococcus`, `Faecalibacterium`, and `Staphylococcus`. GG2 recovered the
same broad pattern: enriched inflammatory or enterobacterial structure and
depleted `Prevotella` and `Streptococcus` states. Depth 2 again sharpened the
signal, especially within `Bacteroides`-centered subtypes, where a major
subtype reached adjusted OR 2.25 with q = 1.2e-23. Four autoimmune associations
survived sensitivity analysis.

### Discussion

Because this AGP label pools multiple autoimmune conditions, the safest reading
is an inflammatory umbrella rather than a disease-specific microbiome
signature. The observed `Bacteroides` / `Morganella` enrichment and
`Faecalibacterium` / `Prevotella_9` depletion are broadly compatible with
disease-specific literature in RA and T1D, but the component mechanisms differ:
`Prevotella` expansion in RA is not the same biological claim as early-life
dysbiosis in T1D, even though both can be placed within a broader Treg/Th17 or
immune-education framework [@Scher2013_Prevotella_RA; @Vatanen2018_T1D;
@Atarashi2013_Clostridia_Treg; @DeLuca2019_autoimmune]. The parallel with IBD
is therefore informative at the level of shared inflammatory ecology, not at
the level of one causal autoimmune taxon.

## 4.4 Acid reflux / GERD

### Introduction

GERD is increasingly discussed in microbiome terms rather than purely as an
acid-mediated disorder, with gram-negative-enriched states and loss of
gram-positive dominance implicated in esophageal disease
[@Deshpande2018_esophageal; @DeSouza2021_GERD; @Wang2024_GERD_MR].

### Results

In SILVA, acid reflux showed adjusted enrichment of `Bacteroides`
(adjusted OR 1.28, q = 6.9e-06) and depletion of `Streptococcus`
(adjusted OR 0.45, q = 4.2e-03) and `Staphylococcus`
(adjusted OR 0.41, q = 1.1e-02). GG2 replicated the depletion of
`Streptococcus` and `Staphylococcus` together with a depleted `Prevotella`
signal. Depth 2 produced additional `Parabacteroides`- and `Alistipes`-linked
signals. After contaminant-aware filtering, only the `Bacteroides` association
remained.

### Discussion

The reflux results are best interpreted as gut-side correlates of a literature
that is itself centered on esophageal and oral compartments
[@Deshpande2018_esophageal; @DeSouza2021_GERD]. The surviving `Bacteroides`
signal therefore suggests a modest reflux-spectrum association rather than a
direct esophageal mechanism, whereas the oral-associated taxa are more exposed
to compartment choice, PPI effects, and ecological confounding
[@Wang2024_GERD_MR]. The clean analysis is useful here because it narrows the
claim to one stool-side association rather than leaving a mixed oral-plus-gut
story on the table.

## 4.5 Cardiovascular disease

### Introduction

CVD is one of the most mechanistically developed microbiome literatures because
of the TMAO pathway, oral-gut cross-talk, and inflammation-thrombosis links
[@Wang2011_TMAO; @Tang2013_TMAO; @Jie2017_ASCVD; @Zhu2016_TMAO_thrombosis].

### Results

The univariate CVD signals were strong, especially `Pasteurellaceae`,
`Staphylococcus`, and `Neisseria`, but none survived covariate adjustment at
depth 1. GG2 showed analogous univariate enrichment of oral-associated or
aerotolerant taxa. Depth 2 nonetheless yielded a notable `Agathobacter`
association (adjusted OR 3.47, q = 2.6e-05). No CVD signal survived the
sensitivity analysis.

### Discussion

CVD illustrates why the covariate-adjusted pass matters. Without adjustment,
one could easily overstate the oral-associated signals. After adjustment and
contaminant-aware filtering, the evidence is much weaker, so this section is
best read as hypothesis-generating rather than as a primary biological claim.

## 4.6 Obesity

### Introduction

Obesity has a large and often contradictory microbiome literature, ranging from
energy-harvest hypotheses to gene-richness models and `Akkermansia`-linked
barrier effects [@Turnbaugh2006_energy; @Ley2006_obesity;
@LeChatelier2013_richness; @Everard2013_Akkermansia; @Walters2014_meta_obesity].

### Results

For obesity, the evidence remained at the univariate depth-1 level, led by
`Prevotella_7` enrichment (OR 2.75) and weaker shifts in
`Escherichia-Shigella` and `Lactobacillus`. No depth-1 association survived
covariate adjustment, GG2 mainly reproduced a broad `Prevotella`-side signal,
no depth-2 associations reached significance, and nothing survived the
sensitivity analysis.

### Discussion

These results suggest that dominant-state typing captures some obesity-related
structure, but that the signal is either confounded by demographic variables or
too diffuse to stabilize at the level of top-taxon dominance. That is itself
useful: it shows where the DCST framework should remain modest rather than
forcing every phenotype into a clean association story.

## 4.7 Seasonal allergies

### Introduction

Seasonal allergy is another umbrella phenotype, but one with unusually clear
microbiome priors through the hygiene hypothesis and early-life microbiota
maturation literature [@Strachan1989_hygiene; @Hua2016_AGP_allergy;
@Hoskinson2023_allergy].

### Results

In SILVA, seasonal allergies showed adjusted depletion of `Prevotella_9`
(adjusted OR 0.67, q = 3.3e-08) and adjusted enrichment of
`Escherichia-Shigella` (adjusted OR 1.13, q = 5.1e-03) and
`Enterobacteriaceae` (adjusted OR 1.21, q = 2.7e-02). GG2 showed the same broad
pattern, especially strong `Prevotella` depletion. Depth 2 added
`Bacteroides`-centered subtype signal, and all three adjusted depth-1
associations survived sensitivity analysis.

### Discussion

Among the non-IBD outcomes, seasonal allergy is one of the cleaner
screen-level signals in the AGP scan. The protective `Prevotella` pattern fits
well with the broader allergy literature, which links allergic phenotypes to
lower richness, delayed microbiota maturation, and weaker SCFA-associated
immune education [@Hua2016_AGP_allergy; @Hoskinson2023_allergy;
@Strachan1989_hygiene]. The Proteobacteria-side enrichment is likewise
compatible with a more dysbiotic, less tolerogenic gut ecosystem. Because the
signal survives sensitivity analysis and replicates qualitatively in GG2, it is
one of the stronger non-IBD findings. Even so, it is best read as an ecological
bridge between adult AGP allergy phenotypes and the early-life maturation
story, not as a definitive allergic-rhinitis-specific signature.

## 4.8 Diabetes

### Introduction

Diabetes is mechanistically rich but clinically heterogeneous. The literature
links diabetes to reduced butyrate producers, endotoxemia, medication effects,
and altered fermentation ecology [@Qin2012_T2D; @Karlsson2013_T2D_women;
@Cani2007_endotoxemia; @Forslund2015_metformin; @Takeuchi2023_carbohydrate].

### Results

At depth 1, diabetes yielded only univariate depletion of `Streptococcus` and
`Staphylococcus`, with no adjusted-significant depth-1 signal. GG2 again
recovered only modest, taxonomically heterogeneous univariate shifts. Depth 2,
however, produced a strong `Akkermansia`-associated subtype signal
(adjusted OR 4.88, q = 2.1e-04). No diabetes association survived the
sensitivity analysis.

### Discussion

The `Akkermansia` depth-2 result is biologically interesting precisely because
it is paradoxical if interpreted naively. Since metformin can enrich
`Akkermansia`, the cross-sectional signal may encode medication rather than
disease mechanism. The broader diabetes literature is dominated by T2D,
glucose-tolerance gradients, and treatment confounding, so a clean disease-only
signal is not expected in a cross-sectional screen of this kind
[@Qin2012_T2D; @Karlsson2013_T2D_women; @Forslund2015_metformin]. The present
pattern therefore reads sensibly: no stable adjusted depth-1 signal, but one
striking `Akkermansia`-associated depth-2 lead that likely mixes disease
status, treatment, and carbohydrate-metabolic ecology
[@Forslund2015_metformin; @Takeuchi2023_carbohydrate; @Cani2007_endotoxemia].

## 4.9 *Clostridioides difficile* infection (CDI)

### Introduction

CDI is a dysbiosis-driven disease in which colonization resistance is lost after
antibiotic disruption, making it a natural target for community-state analysis
[@Schubert2014_CDI; @Seekatz2014_FMT; @Sehgal2021_CDI].

### Results

In SILVA, CDI showed adjusted enrichment of `Escherichia-Shigella`
(adjusted OR 1.48, q = 1.1e-02). GG2 highlighted related inflammatory or
facultative-anaerobe structure, including `Proteus`. Depth 2 added a strong
`Sutterella`-linked subtype signal (adjusted OR 9.85, q = 1.2e-05). After
sensitivity analysis, `Escherichia-Shigella` remained enriched and
`Prevotella_9` showed a protective pattern.

### Discussion

CDI is one of the stronger non-umbrella disease sections because the disease
itself is defined by collapse and recovery of colonization resistance
[@Schubert2014_CDI; @Seekatz2014_FMT; @Sehgal2021_CDI]. The
enterobacterial enrichment is in line with recurrent-CDI dysbiosis, and the
protective `Prevotella_9` direction is at least ecologically compatible with
restoration of anaerobic fermentation and bile-acid-converting community
structure after FMT. The `Sutterella` depth-2 signal is interesting, but it is
best treated as exploratory until it can be tied more directly to recurrence,
bile acids, or treatment status.

## 4.10 Lung disease

### Introduction

The gut-lung axis links microbial fermentation, immune tone, and airway
phenotypes, but the AGP "lung disease" label is broad and probably pools several
different conditions [@Budden2017_gutlung; @Trompette2014_fiber_asthma;
@Horn2022_Prevotella_lung].

### Results

In SILVA, the main adjusted depth-1 lung-disease signal was `Staphylococcus`
enrichment (adjusted OR 1.55, q = 3.9e-02), while univariate depletion appeared
for `Prevotella_7`, `Streptococcus`, and `Pasteurellaceae`. GG2 reproduced the
depleted `Prevotella` / `Streptococcus` / `Pasteurellaceae` side of the signal.
Depth 2 added `Agathobacter` and `Staphylococcus` subtype associations. None of
the signals survived sensitivity analysis.

### Discussion

This phenotype is too broad to support a strong disease-specific claim, and the
comparison literature itself is split between gut-origin asthma studies, airway
COPD studies, and respiratory-infection models [@Budden2017_gutlung;
@Trompette2014_fiber_asthma; @Horn2022_Prevotella_lung]. The surviving
`Staphylococcus` signal is therefore difficult to interpret as a gut-lung-axis
result per se, especially because it disappears in the clean analysis. The most
defensible reading is that the AGP screen detects a weak, partly
compartment-mixed respiratory signal rather than a coherent lung-disease
microbiome signature.

## 4.11 Migraine

### Introduction

Migraine is increasingly discussed through gut-brain signaling, inflammation,
and microbial metabolite effects, but the current microbiome evidence base is
thinner than for IBD or allergy [@Cryan2012_gutbrain; @Gorenshtein2025_migraine;
@Kappeter2023_migraine_FMT].

### Results

Migraine produced one adjusted depth-1 signal in SILVA: `Bacteroides`
enrichment (adjusted OR 1.16, q = 3.4e-02). GG2 showed a qualitatively similar
`Bacteroides`-side enrichment together with weaker secondary signals. No
depth-2 association reached significance, and nothing survived the sensitivity
analysis.

### Discussion

This remains exploratory. The direct migraine microbiome literature is still
small and heterogeneous, with current reviews emphasizing dysbiosis, SCFA and
barrier pathways, and gut-brain-axis signaling more than any one reproducible
taxon [@Gorenshtein2025_migraine; @Kappeter2023_migraine_FMT;
@Cryan2012_gutbrain]. The present `Bacteroides` enrichment is therefore best
read as a weak cross-disease inflammatory echo rather than as
migraine-specific evidence.

## 4.12 Kidney disease

### Introduction

The gut-kidney axis is mechanistically attractive because microbial metabolites
such as TMAO, indoxyl sulfate, and p-cresyl sulfate can worsen renal decline
[@Tang2015_TMAO_CKD; @Vaziri2013_CKD].

### Results

Kidney disease showed no univariate depth-1 signals, but the adjusted pass
identified strong enrichment of `Lactobacillus` (adjusted OR 5.85, q = 7.8e-03)
and `Prevotella_7` (adjusted OR 4.57, q = 1.0e-02). GG2 recovered a compatible
pattern through enriched `Prevotella` and a reclassified `Bacteroides` group.
Depth 2 retained a `Lactobacillus` signal, and both adjusted depth-1
associations survived sensitivity analysis.

### Discussion

Kidney disease offers a useful example of a signal that is statistically
stronger than its immediate biological label. The CKD literature emphasizes
TMAO, uremic toxins, and depletion of saccharolytic commensals more than it
emphasizes `Lactobacillus` as a universal pathogenic signal
[@Tang2015_TMAO_CKD; @Vaziri2013_CKD]. The observed `Lactobacillus` and
`Prevotella_7` enrichments may therefore encode treatment behavior, probiotic
use, or compensatory ecology superimposed on a true gut-kidney-axis signal.
What keeps this section prominent is robustness: unlike several weaker
phenotypes, the adjusted depth-1 associations survive the clean analysis.

# 5. External Validation Across Two Complementary Cohorts

## 5.1 HMP2 / IBDMDB

HMP2 / IBDMDB (`PRJNA398089`) contributed the clinically grounded inflammatory
validation arm. After phenotype filtering, the local cohort contained 166
samples: 45 healthy controls, 81 Crohn's disease samples, and 40 ulcerative
colitis samples. The caveat is ecological: in the local metadata, the retained
set is dominated by rectal, ileal, and colonic sites and includes repeated
measures. At depth 1, the HMP2 landscape was too sparse to yield strong
multiple-testing-corrected disease signal. At depth 2, however, a
`Bacteroides__Faecalibacterium` subtype was depleted in `IBD_vs_Healthy`
(OR = 0.391, q = 0.016), `UC_vs_Healthy` (OR = 0.318, q = 0.032), and
`Crohn_vs_Healthy` (OR = 0.429, q = 0.041).

## 5.2 PRJEB84421

PRJEB84421 supplied the stool-based complement. It retained 73 samples split
across 20 healthy controls, 24 Crohn's disease samples, and 29 pediatric
orofacial granulomatosis samples. This cohort is a better ecological match to
AGP stool than HMP2, but it is smaller and phenotype-mixed. Its strongest
signals were therefore more modest and cohort-specific: `RARE_DOMINANT` was
enriched in `OFG_vs_Healthy` at depth 1 (OR = 4.93, q = 0.062) and in the
combined inflammatory-vs-healthy comparison at depth 2 (OR = 4.94, q = 0.084).

**Table 5. External validation cohort summary.**

\begingroup\small

| Cohort | n | Phenotype split | Main caveat | Main signal | Interpretation |
| --- | --- | --- | --- | --- | --- |
| HMP2 / IBDMDB | 166 | 45 Healthy / 81 Crohn / 40 UC | Mucosal gut sites; repeated measures | No depth-1 q<0.05; depth-2 depletion of a Bacteroides/Faecalibacterium subtype | Directional / structural replication |
| PRJEB84421 / OFGCD-FI-2025 | 73 | 20 Healthy / 24 Crohn / 29 OFG | Stool-based but pediatric and phenotype-mixed | Rare-dominant enrichment in OFG vs healthy and inflammatory vs healthy | Stool-based complement with cohort-specific signal |

\endgroup

## 5.3 Interpretation

Taken together, the two external cohorts support a narrower but more defensible
validation claim than exact taxon-by-taxon replication. HMP2 shows that the
DCST framework captures clinically relevant inflammatory subtype structure in an
IBD-focused longitudinal cohort, while PRJEB84421 shows that the same framework
remains interpretable in a smaller stool-based inflammatory cohort. The
validation therefore supports **structural and directional generalization**
rather than exact AGP-style replication. In practice, this strengthens the
IBD-centered part of the Phase 1 story much more than it strengthens the more
tentative umbrella outcomes.

![Two-cohort external validation summary. HMP2 provides clinically grounded IBD validation, while PRJEB84421 provides a stool-based inflammatory complement. Together they support structural and directional generalization rather than exact replication.](/Users/pgajer/current_projects/linf/dev/FIGURE_V1_two_cohort_external_validation.png){ width=92% }

# 6. General Discussion

The disease-wise organization is not merely editorial. It prevents the paper's
biological interpretation from being flattened into one pooled list of
odds ratios. Once the results are viewed disease by disease, a more coherent
picture emerges. IBD, seasonal allergies, IBS, CDI, and kidney disease carry
the clearest depth-1 or depth-2 signals. Autoimmune disease is also strong
statistically, but its AGP phenotype is broad and should be interpreted as an
umbrella inflammatory label. CVD, obesity, lung disease, diabetes, and migraine
are more mixed: they contain usable signals, but those signals are either
attenuated by adjustment, sensitive to contamination-aware filtering, or too
coarse for confident mechanistic interpretation.

Several cross-cutting patterns survive this disease-wise reading. First,
`Bacteroides`-dominant structure repeatedly aligns with inflammatory or
inflammation-adjacent phenotypes, including IBD, autoimmune disease, acid
reflux, and migraine. Second, `Prevotella`-side dominance is repeatedly
protective in allergy and inflammatory settings, especially seasonal allergies
and IBD. Third, rare but striking enterobacterial states, particularly
`Morganella`, behave as markers of severe dysbiosis rather than as subtle
population shifts. Fourth, depth 2 adds real information by turning a
top-dominant label into an ordered ecological pair. For a mathematically
oriented reader, this is important: the refinement is explicit and interpretable
rather than latent.

Those broad patterns also line up with the stronger disease-specific literatures
in recognizable ways. In IBD, the combination of enterobacterial rare states,
`Prevotella_9` depletion, and depth-2 `Bacteroides` heterogeneity matches the
treatment-naive and multi-omics literature that emphasizes loss of butyrate-
associated anaerobes together with inflammatory expansion of facultative taxa
[@Gevers2014_CD; @LloydPrice2019_IBD; @Sokol2008_Fprausnitzii]. IBS shows a
weaker but directionally similar ecology, with Enterobacteriaceae-side
enrichment and `Faecalibacterium` depletion rather than a dramatic dysbiosis
collapse [@Pittayanon2019_IBS; @Tap2017_IBS]. Seasonal allergy aligns with the
microbiota-maturation and immune-education literature, in which reduced
tolerogenic fermentation and delayed ecological maturation accompany allergic
phenotypes [@Hua2016_AGP_allergy; @Hoskinson2023_allergy; @Strachan1989_hygiene].
CDI fits the colonization-resistance framework of bile acids, donor-like
recovery, and Proteobacteria blooms [@Schubert2014_CDI; @Seekatz2014_FMT;
@Sehgal2021_CDI]. Kidney disease is the least mechanistically settled of these
stronger signals, but it still sits plausibly within the broader gut-kidney
axis defined by TMAO and uremic-toxin biology [@Tang2015_TMAO_CKD;
@Vaziri2013_CKD].

The paper should still be careful about overclaiming. AGP phenotypes are
self-reported. Some disease bins are umbrella categories. The analysis is
cross-sectional, so reverse causation, medication effects, and residual
confounding remain live possibilities. The sensitivity analysis shows that some
signals rely on taxa that are plausibly oral or skin associated. And the
validation cohorts, while valuable, are intentionally heterogeneous rather than
perfect one-to-one replicas of AGP stool. Those limitations argue for a
measured conclusion: DCSTs organize the data in a biologically useful way, but
they do not by themselves establish mechanism or causality. For the same
reason, the disease-specific literature cited in Section 4 should be read as
mechanistic context for the discovery-screen phenotypes rather than as proof
that the AGP labels are clinically identical to those literature-defined
entities.

# 7. Conclusion

Taken together, this Phase 1 manuscript is best understood as a disease-wise
atlas built on a common dominant-state framework. DCSTs recover a structured gut landscape,
remain fairly stable across taxonomic schemas, and identify a nontrivial set of
phenotype associations that survive adjustment and, in many cases,
contamination-aware filtering. The strongest part of the overall story concerns
inflammatory phenotypes, especially IBD and related subtype structure, and that
part is now supported by two complementary external validation cohorts. This is
already enough to justify DCSTs as a useful analytical language for future gut
microbiome studies and for the broader `linf` application paper.

# 8. Data and Code Availability

The manuscript source, manuscript figures, and manuscript table builders used
for this draft are maintained in the public `linf` repository:
<https://github.com/pgajer/linf>. The discovery analysis is based on the
American Gut Project stool cohort (`PRJEB11419`) through PRIME-derived
species-level abundance tables, and the external validation analyses use
`PRJNA398089` (HMP2 / IBDMDB) and `PRJEB84421`. The specific processed result
tables used to generate the figures and Tables 1-4 are those analyzed in the
companion gut-microbiome working tree during manuscript preparation.

# 9. Funding

Research reported in this publication was supported by grant number
`INV-048956` from the Gates Foundation.

# 10. Competing Interests

The authors declare no competing interests.

# 11. Acknowledgments

We thank the participants and investigators of the American Gut Project, HMP2 /
IBDMDB, and PRJEB84421 for making these cohort resources available, and we
thank the PRIME project context that made the large-scale gut screen feasible.

# 12. References
