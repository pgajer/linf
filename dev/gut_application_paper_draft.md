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
  - \usepackage{float}
  - \usepackage{titling}
  - \usepackage{caption}
  - \captionsetup{font=small,labelfont=bf}
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
$L_\infty$ normalization and dominant-community state types (DCSTs), in which
each sample is assigned by rank structure rather than by cluster membership.
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
defined from the rank structure of a sample after $L_\infty$ normalization. At
depth 1, the dominant taxon determines the top-level state. At depth 2, the
second-ranked taxon refines that state into an ordered subtype. The resulting
construction is deterministic, hierarchical, and local to the sample: the label
is read directly from the ordered dominant coordinates rather than inferred from
a fitted global partition.

The gut is a particularly good stress test for this idea. In the vaginal
microbiome, strong dominance can be common and even canonical. In the gut,
strong dominance is less universal and often appears in dysbiotic or otherwise
unusual settings. That makes the gut harder to type cleanly, but also more
interesting scientifically. Literature on IBS, IBD, autoimmune disease, and
cardiometabolic conditions suggests that shifts in dominant taxa or dominant
taxon pairs can mark altered ecological states [@Tap2017_IBS;
@Pittayanon2019_IBS; @LloydPrice2019_IBD; @Gevers2014_CD;
@Scher2013_Prevotella_RA; @Vatanen2018_T1D; @Wang2011_TMAO; @Tang2013_TMAO].

Our aim is not to argue that dominant taxa alone explain disease biology. The
claim is narrower: a clustering-free dominance framework can recover structured
gut community states at scale, support disease-wise interpretation without
collapsing every phenotype into one pooled signal, and remain informative under
covariate adjustment, sensitivity analysis, taxonomy change, and external
validation.

> **What is a DCST?**
>
> A DCST is the ordered dominant part of a sample after $L_\infty$
> normalization. At depth 1, only the largest component matters. At depth 2, the
> second-largest component refines the parent class. Samples that do not meet
> the prevalence threshold for a named dominant type are grouped into
> `RARE_DOMINANT`.

![Conceptual schematic of the DCST construction. A sample is first normalized by its largest component, then assigned a depth-1 state by the top-ranked taxon and a depth-2 state by the ordered dominant pair. Rare dominant patterns below the prevalence threshold are grouped into `RARE_DOMINANT`.](/Users/pgajer/current_projects/linf/dev/FIGURE_1_dcst_conceptual_schematic.png){ width=92% }

# 2. Methods

## 2.1 Cohort and preprocessing

The discovery analysis used AGP gut samples from the PRIME-derived SILVA
species-level abundance table. The starting table contained 34,711 samples, of
which 30,290 passed quality-control filters and were retained for analysis. The
retained matrix contained 274 taxa. Age was available for 27,372 samples
(90.4\%), sex for 25,598 (84.5\%), and BMI for 26,541 (87.6\%), producing
22,632 complete cases for adjusted modeling.

## 2.2 DCST construction

Each sample was normalized by its largest component, so the dominant taxon took
value 1 and the remaining coordinates were interpreted relative to that
reference. A depth-1 DCST was defined by the dominant taxon. A depth-2 DCST was
defined by the ordered pair consisting of the dominant and second-ranked taxa.
To keep the resulting state space interpretable, only dominant patterns meeting
a minimum prevalence threshold $n_0$ were promoted to named DCSTs; the rest were
grouped into `RARE_DOMINANT`.

## 2.3 Association analysis

We first performed univariate screening with Fisher's exact test. We then fit
adjusted logistic models for disease indicators using the DCST label together
with available covariates. For a mathematically oriented reader, the odds ratio
may be read as follows: an odds ratio above 1 means the DCST is enriched in the
given phenotype relative to the reference group, while an odds ratio below 1
means depletion. Multiple testing was handled by the Benjamini-Hochberg
procedure, and q-values below 0.05 were treated as statistically notable rather
than as proof of mechanism.

## 2.4 Sensitivity analysis and taxonomy comparison

Because some nominal gut signals can be driven by oral- or skin-associated
contaminants, we repeated the adjusted analysis after excluding samples
dominated by a prespecified contaminant set. We also ran a parallel dominant-type
analysis under GG2 taxonomy and compared dominant-genus assignments between
SILVA and GG2 on the overlapping sample set.

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
background rather than as isolated anecdotes.

![Full-cohort AGP landscape. Panel A shows the leading depth-1 DCSTs after truncation at n0 = 50, while Panel B shows the empirical dominance strength before $L_\infty$ normalization.](/Users/pgajer/current_projects/linf/dev/FIGURE_2_full_cohort_landscape.png){ width=92% }

![Adjusted depth-1 association overview. Panel A shows all taxa that survive adjustment at least once across the disease screen, with colored cells annotated by adjusted odds ratio and blank cells indicating no surviving adjusted signal. Panel B shows the strongest adjusted associations with 95% confidence intervals.](/Users/pgajer/current_projects/linf/dev/FIGURE_3_adjusted_association_overview.png){ width=98% }

Looking only at the counts, the sensitivity analysis reduces the number of
adjusted-significant depth-1 associations from 28 to 23. A more informative
view is to compare which signals are lost, which are retained, and which become
clearer once likely oral- or skin-associated dominants are excluded. Twenty
adjusted-significant associations are retained, eight are lost, and three
become significant only in the clean analysis. The lost signals are concentrated
in `Staphylococcus`- and `Streptococcus`-linked contrasts, while the strongest
IBD, allergy, and several `Bacteroides` / `Prevotella_9` associations remain in
place.

![Contamination-aware robustness summary. The clean analysis removes several oral-associated contrasts but preserves the strongest inflammatory and allergy-linked signals.](/Users/pgajer/current_projects/linf/dev/FIGURE_4_contamination_robustness.png){ width=98% }

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

![Cross-taxonomy concordance summary. Harmonized agreement between SILVA and GG2 is high, while exact short-label equality is much lower because a few large dominant states are reclassified rather than contradicted.](/Users/pgajer/current_projects/linf/dev/FIGURE_6_taxonomy_concordance.png){ width=95% }

\clearpage

# 4. Disease-Wise Results and Interpretation

The subsections in this section follow a common pattern. Each begins with a
brief literature frame, then summarizes the AGP signal across SILVA, GG2, depth
2, and sensitivity analysis, and finally interprets the result at a level that
matches the strength of the evidence. This repeated structure is intentional: it
lets the reader compare diseases without losing track of which claims are
strong, which are tentative, and which are mainly useful for context.

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
cautiously.

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

IBD is the clearest case in the current paper where depth 2 matters scientifically and
mathematically. At depth 1, `Bacteroides` is a broad parent label. At depth 2,
that parent is split into subtype structure with materially different disease
behavior. The rare but very strong `Morganella` signal is also consistent with
an inflammatory oxygenation story in which facultative Enterobacteriaceae gain
ecological advantage in the inflamed gut [@Sartor2008_IBD; @Frank2007_IBD].
Section 5 returns to IBD with external validation, because this is the disease
family for which the validation cohorts are most informative.

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
is ecological rather than disease-specific. Even so, the `Bacteroides` /
`Prevotella` / `Faecalibacterium` pattern fits an immune-regulation narrative in
which SCFA-associated commensals help maintain tolerance while inflammatory
states favor more dysbiotic community structures. The parallel with IBD is
notable and may reflect shared inflammatory mechanisms rather than identical
etiology.

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

The reflux results are interesting but should be handled conservatively. The
surviving `Bacteroides` signal suggests a real gut-side association, whereas the
protective oral-associated taxa may partly reflect shared sampling or ecological
confounding. This is exactly the sort of setting where the sensitivity analysis
helps distinguish potentially robust gut patterns from signals that may be
driven by translocated oral structure.

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
contaminant-aware filtering, the evidence is much weaker. In the manuscript,
this section should therefore read as hypothesis-generating rather than as a
primary biological claim.

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
useful: it tells us where the DCST framework is modest rather than pretending
every phenotype produces a clean result.

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

Among the non-IBD outcomes, seasonal allergy is one of the cleaner signals in
the AGP screen. The protective
`Prevotella` pattern fits well with the allergy literature, while the
Proteobacteria-side enrichment is consistent with a more dysbiotic, less
tolerogenic gut ecosystem. Because the signal survives sensitivity analysis and
replicates qualitatively in GG2, it deserves a prominent but still cautious
place in the paper.

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
disease mechanism. It is therefore a useful reminder that the disease-wise
narrative must separate association from etiology.

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

CDI is one of the stronger non-umbrella disease sections in the current paper. The
enterobacterial enrichment is in line with the CDI literature, and the
protective `Prevotella_9` direction is at least ecologically compatible with the
idea that diverse anaerobic states support colonization resistance. The
`Sutterella` depth-2 signal is novel and should be presented as exploratory.

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

This phenotype is too broad to support a strong disease-specific claim. The signal
is interesting as a gut-lung-axis hint, but because the phenotype is umbrella
level and the sensitivity analysis removes the depth-1 association, it belongs
in the paper as secondary evidence rather than as a central discovery.

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

This is a modest result, but still a usable one. The safest framing
is that migraine contributes to the cross-disease pattern in which
`Bacteroides`-dominant states tend to align with inflammatory or
neuroinflammatory phenotypes, while the mechanistic interpretation remains
speculative.

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

Kidney disease offers a useful example of a signal that is statistically stronger than
its immediate biological interpretation. `Lactobacillus` is not an obvious CKD
pathogen, so the enrichment could reflect probiotic use, compensatory ecology,
or disease-associated treatment behavior. The most honest presentation is that the
result is robust and interesting, but mechanistically unresolved.

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

The paper should still be careful about overclaiming. AGP phenotypes are
self-reported. Some disease bins are umbrella categories. The analysis is
cross-sectional, so reverse causation, medication effects, and residual
confounding remain live possibilities. The sensitivity analysis shows that some
signals rely on taxa that are plausibly oral or skin associated. And the
validation cohorts, while valuable, are intentionally heterogeneous rather than
perfect one-to-one replicas of AGP stool. Those limitations argue for a
measured conclusion: DCSTs organize the data in a biologically useful way, but
they do not by themselves establish mechanism or causality.

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

# References
