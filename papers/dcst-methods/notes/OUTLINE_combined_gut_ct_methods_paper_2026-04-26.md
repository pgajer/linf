# Proposed Outline: Combined Gut+CT dCST Methods Paper

## Recommendation

Create a new umbrella manuscript in:

`/Users/pgajer/current_projects/linf/papers/dcst-methods`

rather than trying to keep extending:

`/Users/pgajer/current_projects/linf/papers/gut-dcst`

This should be treated as a new paper, not as a direct revision of the gut application manuscript.

## Why this combination makes sense

The current gut manuscript demonstrates that dCSTs can provide an interpretable representation of heterogeneous gut dysbiosis and can be transferred across external cohorts. Its main weakness is that the strongest signal is sample-level and attenuates under one-sample-per-subject sensitivity.

The CT manuscript provides a complementary strength:

- a second microbiome niche
- a different data layer and analysis regime
- stronger biologically interpretable outcome-associated states
- richer follow-up through functional annotation

Together, the two applications support a broader methods claim:

> dCSTs provide a deterministic, interpretable representation of microbiome state structure that is useful across distinct biological settings, data types, and downstream inferential goals.

This is substantially stronger than a gut-only submission.

## What to borrow from the original Linf paper

Borrow only the minimum framework needed to motivate the method:

- deterministic sample-wise state assignment
- rank-based invariance to within-sample rescaling
- stability relative to adding or removing other samples
- avoidance of clustering-defined label instability

Do **not** import the full geometric program, extensive cube-embedding development, or broad L-infinity theoretical machinery unless it is directly needed to justify one method step in the manuscript.

The new paper should read as:

- a methods-forward microbiome paper with two strong applications

not as:

- a mathematical theory paper with two appended case studies

## Core manuscript thesis

### One-sentence thesis

dCSTs are a deterministic, rank-based alternative to clustering-defined CST/enterotype-style labels that yields interpretable sample-level microbiome states, supports cross-cohort transfer, and remains useful across distinct microbiome domains.

### Main claims

1. dCSTs define sample-level states without dependence on clustering model choice or global sample composition.
2. The framework is portable enough to support external transfer in gut IBD, while also revealing where portability weakens.
3. The framework is biologically informative in an independent vaginal metagenomic CT application, where it recovers function-linked outcome-associated states.
4. The combined evidence supports dCSTs as a reusable representation layer rather than as a niche gut-only construct.

## What each application should contribute

### Gut application should contribute

- motivation from the instability and interpretability limits of enterotypes
- deterministic taxonomic dCST construction in a large AGP-derived ecological background
- transfer to diagnosis-stronger external cohorts
- an honest view of strengths and limits of transfer
- one or two clinically interpretable within-disease follow-ups

The gut application should carry the paper's:

- portability claim
- taxonomy-universe harmonization claim
- sensitivity-analysis realism

### CT application should contribute

- a clearly independent validation domain
- a non-gut microbiome niche
- a stronger biologically structured association signal
- functional interpretation of signal-bearing states
- evidence that dCSTs are useful beyond taxonomic stool dysbiosis summaries

The CT application should carry the paper's:

- biological richness
- outcome-linked interpretability
- multi-layer follow-up value

## Proposed manuscript structure

## 1. Title options

### Option A

Deterministic Dominant-Community State Types Provide an Interpretable Microbiome Representation Across Gut and Vaginal Applications

### Option B

Dominant-Community State Types as a Portable and Interpretable Framework for Microbiome State Analysis

### Option C

A Deterministic Rank-Based Framework for Interpretable Microbiome State Typing Across Cohorts and Niches

## 2. Abstract

Recommended abstract logic:

1. Start with the general problem:
   clustering-defined microbiome state labels are useful but unstable and difficult to compare across analyses.
2. Introduce dCSTs as a deterministic rank-based alternative.
3. State the two applications:
   gut IBD external transfer and vaginal CT outcome association.
4. Summarize the different proof points:
   gut demonstrates transfer and its limits; CT demonstrates biologically interpretable and function-linked state associations.
5. Close with the representation-layer claim rather than a diagnostic-biomarker claim.

## 3. Introduction

### Section goals

- motivate why microbiome studies need sample-level state labels
- explain why clustering-defined labels are unstable
- introduce dCSTs in a concise, accessible way
- justify the need for more than one application domain

### Recommended subsection flow

1. CSTs and enterotypes as useful but clustering-defined constructs
2. Need for deterministic sample-level state labels
3. Short methods motivation from Linf ideas
4. Why two applications are used here
   - gut for transfer and cohort portability
   - CT for independent biological and functional generalization

## 4. Methods

### 4.1 dCST framework

Include:

- retained-feature ranking
- depthwise dominance-lineage construction
- absorb policy
- why assignment is deterministic
- invariance and interpretability properties

Keep this section clear and compact. This is where a few carefully chosen ideas from the original Linf paper belong.

### 4.2 Gut application design

Include:

- AGP local-QZA SILVA taxonomy harmonization
- hierarchy learning cohort
- external validation cohorts
- direct-transfer and rebuilt-cohort modes
- subject-level sensitivity

### 4.3 CT application design

Include:

- CT labeled cohort
- CT+reference background cohort
- VOG-cluster and full-VOG dCST layers
- outcome testing restricted to labeled subset
- functional follow-up

### 4.4 Statistical framework

Include:

- label-level association tests
- panel-level summaries
- multiple-testing strategy
- subject-level sensitivity where relevant
- Bayesian sparse-state follow-up where relevant

Do not force both applications into identical downstream statistics if that would weaken either one. The methods section should present a common dCST core plus application-specific inference layers.

## 5. Results

### 5.1 Framework overview and conceptual figure

One main conceptual methods figure should explain:

- dCST assignment
- absorb-policy retained states
- why deterministic rank-based labels differ from clustering labels

### 5.2 Gut application: external transfer and portability

Recommended scope:

- AGP discovery hierarchy summary
- Halfvarson as the main stool external test
- HMP2 and Gevers as boundary cases
- one compact within-disease extension

This section should be shorter and cleaner than the current gut paper. It should emphasize:

- what transfers
- what attenuates under subject-aware sensitivity
- what this teaches us about portability

### 5.3 CT application: outcome-associated states and functional follow-up

Recommended scope:

- VOG-cluster layer as lower-resolution view
- full-VOG dCST layer as primary high-resolution result
- clearance-oriented and persistence-oriented branches
- taxonomy, product, and KEGG follow-up

This section is where the paper gets biological weight.

### 5.4 Cross-application synthesis

This is the section that makes the combined paper coherent.

It should explicitly compare the two applications:

- gut emphasizes transfer across cohorts and clinical-label boundaries
- CT emphasizes internal biological and functional refinement
- in both settings dCSTs provide interpretable sample-level states without cluster instability

This section is essential. Without it, the paper risks reading like two unrelated case studies.

## 6. Discussion

### Key points to hit

1. What dCSTs do better than clustering-defined state labels
2. Why deterministic assignment matters
3. What the gut application teaches about portability and its limits
4. What the CT application teaches about biological interpretability
5. Why the framework should be viewed as a reusable representation layer

### Important limitation language

Be explicit that:

- not all signals survive subject-aware sensitivity in gut IBD
- transfer depends on a shared taxonomy or shared feature universe
- portability and interpretability are distinct strengths
- dCSTs are not by themselves causal or mechanistic evidence

## 7. Figures

### Figure 1

Conceptual dCST framework figure.

### Figure 2

Gut application overview:

- AGP discovery summary
- external transfer summary
- subject-aware sensitivity layer

### Figure 3

CT application overview:

- VOG-cluster and full-VOG outcome-associated states
- branch-level clearance versus persistence story

### Figure 4

CT biological follow-up:

- taxonomy/product/KEGG interpretation

### Figure 5

Cross-application synthesis or methods comparison figure:

- deterministic label properties
- portability versus biological-resolution comparison across the two domains

## 8. Supplement

Put in supplement:

- additional AGP exploratory phenotype screens
- secondary gut cohorts and sensitivity tables
- additional CT capped-feature robustness analyses
- more technical lineage tables
- extended Bayesian outputs
- any limited theory note imported from Linf that is helpful but not needed in the main text

## What to leave out

### Leave out from the gut manuscript

- extended side-branches that do not materially support the method claim
- too much reviewer-driven calprotectin detail
- overly long transfer audit narrative in the main text

### Leave out from the CT manuscript

- details that are important only for a CT-specific standalone biology paper
- excessive branch-by-branch enumeration that does not help the cross-application methods argument

### Leave out from the Linf paper

- most of the general geometric framework
- cube embedding details
- broad non-microbiome ambitions

## Suggested writing strategy

Build this manuscript as a new paper with:

- a fresh abstract
- a fresh introduction
- a new unifying discussion

Then import and rewrite:

- selected gut Results subsections
- selected CT Results subsections
- a short methods motivation distilled from Linf

Avoid trying to merge source text verbatim. It will produce tone drift and duplicated motivation.

## Suggested folder structure

Use:

`/Users/pgajer/current_projects/linf/papers/dcst-methods`

with this layout:

- `manuscript/`
- `notes/`
- `assets/`
- `scripts/`
- `build/`
- `archive/`

This is better than placing the combined paper under `gut-dcst`, because the combined paper is no longer a gut-only application paper.

## My recommendation

If the goal is to produce a stronger `mSystems` methods paper, this combined gut+CT manuscript is a good idea.

If we do it, we should commit to the following editorial choice:

- the paper is about the **dCST framework**
- gut and CT are the **two principal applications**
- Linf contributes only the minimum conceptual framing needed to explain why the framework is deterministic and interpretable

That structure has a much better chance of feeling like one coherent paper than either:

- submitting the current gut manuscript alone, or
- trying to paste the full CT and Linf papers into it.

## Immediate next steps

1. Decide whether CT contributes one main-results section or two.
2. Decide whether the main paper should keep only one gut within-disease follow-up, most likely Crohn location.
3. Draft a new title, abstract, and introduction before migrating any large results text.
4. Inventory which figures can be reused, which need redesign, and which should stay supplement-only.
5. Create the new manuscript source in `manuscript/` only after the cross-application thesis is finalized.
