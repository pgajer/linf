# Manuscript Scaffold: Gut DCST Application Paper

Date: 2026-03-25

## Working Title Options

1. **Dominant Community State Types Reveal Disease-Associated Dominance
   Patterns in the Human Gut Microbiome**
2. **A Clustering-Free, Rank-Based Framework for Gut Microbiome Community
   Typing**
3. **L-infinity Dominant Community State Types for Gut Microbiome Disease
   Association Analysis**

Recommendation: use option 1 for the first preprint unless the methods angle
needs to dominate more strongly.

## One-Sentence Thesis

L-infinity DCSTs provide a deterministic, hierarchical, clustering-free way to
type gut microbiome samples, and in a full American Gut Project cohort they
identify disease-associated dominance patterns that remain informative after
covariate adjustment, contamination sensitivity analysis, and cross-taxonomy
comparison.

## Intended Audience

Primary audience:

- mathematically literate readers
- statistical microbiome readers
- computational microbiome readers

Secondary audience:

- biologically oriented microbiome readers who may not know the `linf` method

## Tone And Exposition Rules

Write the paper so that a mathematician can follow it without prior
microbiome-specific intuition.

That means:

- define CST/DCST/enterotype language early
- explain why clustering-free assignment matters
- explain odds ratios, BH correction, and sensitivity analysis in plain terms
- separate association from mechanism
- keep microbiome jargon interpretable rather than assumed

## Central Claim Hierarchy

### Claim 1: methodological

DCSTs are a deterministic, rank-based, clustering-free community-typing
framework built from L-infinity normalization and dominant-feature assignment.

### Claim 2: empirical

In the full AGP gut cohort, depth-1 and depth-2 DCSTs recover structured,
interpretable dominance patterns rather than an uninformative blur.

### Claim 3: statistical

A nontrivial set of DCST-disease associations persists after covariate
adjustment and contamination-aware sensitivity analysis.

### Claim 4: robustness

The broad association landscape is not an artifact of a single taxonomy:
SILVA and GG2 give high concordance at the dominant-genus level.

### Claim 5: biological interpretation

The strongest and most credible DCST signals align with known or plausible
microbiome-disease mechanisms, while other signals remain exploratory and are
presented as such.

## What This Paper Is

- a methods-application paper
- a gut use-case paper for `linf`
- a bridge between the `linf` methods paper and future validation studies

## What This Paper Is Not

- not a causal paper
- not a definitive disease-mechanism paper
- not a replacement for clinically curated cohorts
- not a full literature review across all diseases

## Recommended Manuscript Structure

### Abstract

Open with the problem:

- gut community typing is often clustering-dependent and unstable
- dominant taxa can be biologically meaningful, especially in dysbiotic tails

Then state:

- full AGP sample size
- depth-1 and depth-2 DCST counts
- number of adjusted-significant associations
- number robust to contaminant exclusion
- cross-taxonomy concordance

Close with:

- what DCSTs add methodologically
- what remains for future validation

### 1. Introduction

#### 1.1 Why another community-typing framework?

- explain instability and arbitrariness of clustering-based workflows
- position DCSTs as deterministic and sample-local

#### 1.2 Why gut is the right stress test

- vaginal microbiome: dominance often healthy and canonical
- gut microbiome: dominance rarer, often pathology-linked, and harder to type

#### 1.3 What this paper contributes

- first full-cohort gut DCST application
- disease association testing at scale
- depth-2 refinement as subdominant context
- cross-taxonomy robustness check

### 2. Methods

#### 2.1 Cohort and preprocessing

- AGP cohort
- inclusion criteria
- filtering parameters
- metadata availability and complete-case counts

#### 2.2 L-infinity normalization and DCST construction

- depth-1 DCSTs
- depth-2 refinement
- role of `n0`
- why this is rank-based rather than cluster-based

#### 2.3 Disease phenotypes and models

- phenotype parsing
- logistic regression specification
- univariate vs adjusted analyses
- BH correction

#### 2.4 Sensitivity analysis

- oral/skin dominant taxa exclusion
- rationale and limitations

#### 2.5 Cross-taxonomy comparison

- SILVA vs GG2
- what concordance means and what it does not mean

### 3. Results

#### 3.1 The gut DCST landscape

Section claim:

The full cohort contains a structured DCST landscape with a few large dominant
types, a meaningful rare tail, and substantial variation in dominance strength.

Needed outputs:

- size distribution
- dominance-strength distribution
- top depth-1 DCST table

#### 3.2 Disease associations at depth 1

Section claim:

Several DCST-disease associations survive covariate adjustment, with strongest
signals concentrated in gastrointestinal, autoimmune, and metabolic/cardiorenal
phenotypes.

Needed outputs:

- adjusted-results heatmap or forest summary
- top association table
- clear distinction between headline and exploratory signals

#### 3.3 Contaminant-aware robustness

Section claim:

A meaningful subset of the adjusted associations remains after excluding
putative oral/skin contaminant-dominated samples.

Needed outputs:

- robust-vs-lost association summary
- explicit discussion of which associations disappear

#### 3.4 Depth-2 refinement reveals hidden heterogeneity

Section claim:

Subdominant context splits broad depth-1 types into subtypes with distinct
disease profiles, especially within large parent DCSTs.

Needed outputs:

- depth-2 association examples
- one or two subtype case studies

#### 3.5 Cross-taxonomy robustness

Section claim:

The dominant-typing picture is largely stable across SILVA and GG2, despite
important naming and resolution differences.

Needed outputs:

- GG2 size distribution
- concordance heatmap
- a short table of per-condition association counts

#### 3.6 External validation across two complementary cohorts

Section claim:

A paired external-validation strategy strengthens the paper more than either
cohort alone: HMP2 contributes clinically grounded IBD validation and
PRJEB84421 contributes a stool-based inflammatory complement, together
supporting structural and directional generalization beyond AGP.

Needed outputs:

- concise two-cohort description:
  - HMP2 caveats: mucosal sites, repeated measures
  - PRJEB84421 caveats: pediatric, OFG phenotype mix
- HMP2 summary emphasizing depth-2 inflammatory signal
- PRJEB84421 summary emphasizing stool-based structural complement
- one combined validation figure
- one combined validation table
- explicit statement that the result is structural / directional validation, not
  exact AGP replication

### 4. Discussion

#### 4.1 What DCST adds

- deterministic assignment
- rank-based interpretation
- hierarchical refinement
- straightforward handling of rare dominance

#### 4.2 Biological interpretation of the most credible findings

Prioritize:

- IBD-related findings
- IBS / Enterobacteriaceae-related findings
- Faecalibacterium protective framing
- Bacteroides / Prevotella nuance
- metabolic/cardiorenal interpretation only where support is adequate

#### 4.3 Alternative explanations

- sample contamination
- phenotype ambiguity
- reverse causation
- unmeasured confounding
- 16S taxonomy artifacts

#### 4.4 Limitations

- self-reported AGP phenotypes
- cross-sectional design
- incomplete covariates
- prevalence filtering
- GG2 resolution issues for Escherichia-Shigella

#### 4.5 Future work

- external validation cohort
- medication and antibiotic metadata
- longitudinal follow-up
- shotgun/metagenomic refinement
- formal study of depth-2 subtype stability

### 5. Conclusion

Keep this tight:

- DCSTs are viable and interpretable in gut
- some associations are robust enough to be scientifically interesting
- the method is useful even where biology remains unresolved

## Must-Have Sidebars Or Short Explanatory Boxes

These can be literal boxes, short subsections, or visually distinct paragraphs.

1. **What is a DCST?**
2. **Why L-infinity rather than CLR/ILR here?**
3. **What does an odds ratio mean?**
4. **Why self-reported AGP phenotypes still have some value, but limited value**

## What To Leave Out Of The First Preprint

- a full equal-depth review of all disease buckets
- overconfident mechanistic claims for small DCSTs
- claims that depend entirely on contaminant-prone taxa
- any implication that the stratified package subset validates the full-cohort
  findings

## Immediate Writing Tasks

1. Draft the abstract from the report plus the claim inventory.
2. Draft the Introduction and Methods first.
3. Build Results directly from the existing full-cohort output files.
4. Only then pull in literature context from the disease-review library.

## Inputs To Use For The Draft

- `dev/Phase1_Full_Cohort_DCST_Analysis_Report.pdf`
- `../gut_microbiome/outputs/dcst_analysis/full_cohort_adjusted_results.csv`
- `../gut_microbiome/outputs/dcst_analysis/depth2_adjusted_results.csv`
- `../gut_microbiome/outputs/dcst_analysis/sensitivity_clean_adjusted_results.csv`
- `../gut_microbiome/outputs/dcst_analysis/silva_vs_gg2_crosstab.csv`
- disease-review assets under
  `../gut_microbiome/outputs/dcst_analysis/disease_reviews/`
