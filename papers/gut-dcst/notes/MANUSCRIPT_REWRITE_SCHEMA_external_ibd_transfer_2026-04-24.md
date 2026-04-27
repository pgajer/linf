# Manuscript Rewrite Schema: AGP Structure Learning Plus External IBD Transfer

Date: 2026-04-24

## Why this rewrite is needed

The current manuscript still relies too heavily on AGP self-reported disease labels for the main IBD-facing story. That framing is vulnerable because AGP provides a very large and useful ecological discovery cohort, but not diagnosis-grade IBD phenotyping. The strongest available biological evidence in this project comes instead from external cohorts with clinically curated IBD labels, especially when CD and UC are separated where the metadata allow it.

The revised manuscript should therefore make a cleaner distinction between:

- **structure learning** in AGP, and
- **disease association testing** in external cohorts with stronger IBD diagnosis data.

This rewrite does not abandon AGP. It repositions AGP as the cohort used to learn a phenotype-agnostic dominance-state structure and then tests whether that structure remains biologically informative when transferred to clinically annotated IBD cohorts.

## Core paper claim after rewrite

The paper should argue the following:

1. A large AGP gut cohort can be used to learn an interpretable dominant-taxon hierarchy without relying on clinical disease labels.
2. That AGP-derived hierarchy can be transferred to external cohorts with stronger IBD diagnosis metadata.
3. In those external cohorts, AGP-derived dominant-taxon states capture reproducible IBD-associated, CD-associated, and UC-associated dominance patterns.
4. The contribution is therefore not discovery of new IBD biology, but an interpretable and portable representation of known heterogeneous IBD dysbiosis.

## Main narrative center

The manuscript should read as an **IBD-centered transfer paper** rather than an AGP phenotype-screen paper.

The paper should answer this question:

> Can a dominant-taxon hierarchy learned in a large general gut cohort be transferred to clinically characterized IBD cohorts, where it resolves reproducible CD- and UC-associated dominance patterns?

That question is stronger, cleaner, and more biologically defensible than asking readers to trust AGP self-reported IBD as the primary disease evidence.

## Data roles after rewrite

### AGP

Use AGP for:

- deriving the absorb dCST hierarchy;
- describing the ecological coverage of that hierarchy;
- optionally showing an exploratory AGP self-reported phenotype screen as motivation for prioritizing IBD.

Do **not** use AGP as the primary source of admissible IBD diagnosis data.

### External main-text IBD cohorts

Use external cohorts for the primary disease-association claims:

- **Halfvarson 2017**:
  best stool-based main validation cohort;
  strongest main-text support for `IBD_vs_Healthy`, `Crohn_vs_Healthy`, and `UC_vs_Healthy`.

- **HMP2 / IBDMDB**:
  strongest clinically grounded IBD cohort;
  useful for high-quality inflammatory validation even though the available PRIME slice is mucosal and longitudinal rather than stool-only.

- **Gevers 2014**:
  strong Crohn-focused cohort;
  useful as an external CD anchor and as a new-onset / treatment-naive complement.

- **PRJEB84421**:
  secondary inflammatory complement only;
  not a clean CD/UC cohort, so it should not carry the central subtype story.

### External secondary cohorts

Keep out of the main center unless later analyses materially improve:

- Jacobs IBS
- AGP autoimmune and other broad non-IBD phenotype branches

These may remain in supplement or notes, but they should not dilute the IBD transfer story.

## Main-text versus supplement boundary

### Main text

Keep in main text:

- dCST concept figure and brief hierarchy explanation;
- AGP-based hierarchy creation and ecological coverage;
- external transfer of AGP-derived labels into Halfvarson, HMP2, and Gevers;
- CD / UC split wherever the external cohort supports it;
- cautious mention of rebuilt-cohort validation only as a secondary comparator.

### Supplement

Move or keep in supplement:

- AGP self-reported phenotype-wide association screens;
- AGP self-reported IBD label-level association tables and sensitivity extensions;
- autoimmune and IBS AGP follow-up narratives;
- rebuilt-cohort external analyses if they are not the primary inferential mode;
- PRJEB84421 details unless needed in the main text as a narrow inflammatory complement;
- all broader non-IBD exploratory material.

### Working notes only

Keep out of the manuscript entirely:

- development-history narration;
- analysis portfolio language;
- future-cohort acquisition ideas;
- Vestergaard-style pooled-resource planning unless and until it becomes real analysis.

## Recommended section-by-section rewrite

### Title

The current title can remain:

- `Interpretable Dominant-Taxon States for Characterizing IBD and Related Gut Microbiome Dysbiosis`

It is still acceptable if the main text is clearly IBD-centered and “related” is interpreted as secondary inflammatory context rather than equal narrative weight.

### Abstract

The abstract should be rebuilt around this sequence:

1. Known problem:
   IBD dysbiosis is heterogeneous and hard to summarize as interpretable sample-level states.
2. Study design:
   derive a dominant-taxon hierarchy in AGP without depending on clinical disease labels.
3. Primary test:
   transfer AGP-derived labels into clinically annotated external IBD cohorts.
4. Main findings:
   strongest support in Halfvarson, selective support in HMP2, Crohn-focused support or informative boundary-setting in Gevers, subtype separation where available.
5. Interpretation:
   dCSTs provide an interpretable representation of established IBD dysbiosis rather than claiming new microbial biomarkers.

The abstract should no longer read as if AGP self-reported IBD is the primary discovery engine for the biological claim.

### Introduction

The Introduction should be tightened around three ideas:

1. IBD microbiome dysbiosis is established but heterogeneous.
2. The key open problem is how to represent that heterogeneity as interpretable sample-level states.
3. AGP provides a large gut ecology resource for learning a dominance hierarchy, which can then be tested in clinically annotated external IBD cohorts.

The current dCST concept figure should remain early because the construct is central and the reader needs it before interpreting transfer results.

The current “Figure 2 preview” language should be revised so the reader is prepared for an external-IBD-centered paper rather than an AGP phenotype-screen paper.

### Results

Recommended new Results structure:

1. **AGP-derived dominance-state hierarchy**
   Short section.
   Focus on hierarchy construction, ecological coverage, and why AGP is useful as a structure-learning cohort.
   If the AGP phenotype-wide screen is mentioned at all, frame it as exploratory prioritization that motivated the external IBD follow-up.

2. **Transfer of AGP-derived dominance states to external IBD cohorts**
   This becomes the heart of the paper.
   Emphasize the frozen AGP-derived hierarchy as the main inferential mode.

3. **Halfvarson 2017: stool-based CD and UC transfer**
   This should likely be the strongest main-text subsection because it is stool-based and supports both CD and UC contrasts.

4. **HMP2 / IBDMDB: clinically grounded inflammatory validation**
   Keep the caveat explicit that this cohort is mucosal and longitudinal in the current slice.
   Present as high-quality inflammatory support, not as a stool-equivalent replication.

5. **Gevers 2014: Crohn-focused external transfer**
   Present as the Crohn-specific external anchor, even if the signal is weaker or more directional than in Halfvarson.
   Weak or null transfer can still be informative if interpreted as a boundary on portability rather than as failure of the framework.

6. **Secondary inflammatory complement**
   If retained in main text, PRJEB84421 should be brief and explicitly secondary.
   Otherwise move it to supplement.

The current main-text subsections on IBS and autoimmune disease should be removed from Results and replaced by one sentence directing readers to the supplement.

### Materials and Methods

Methods should clearly separate:

- AGP hierarchy derivation,
- AGP exploratory phenotype-screen logic,
- external cohort transfer analysis.

The primary external method should be stated as:

- freezing the AGP-derived absorb dCST hierarchy;
- mapping external samples into that hierarchy;
- testing `IBD_vs_Healthy`, `Crohn_vs_Healthy`, `UC_vs_Healthy`, and related clinically admissible contrasts where the cohort permits.

Rebuilt-cohort validation should remain in Methods but be presented as a secondary comparison rather than the central inferential mode.

The AGP phenotype-coding paragraph should stay, but it should explicitly say that these are exploratory self-reported fields and not the main source of admissible IBD subtype inference.

### Discussion

The Discussion should emphasize:

- AGP was useful because it provided broad gut ecological structure, not diagnosis-grade IBD labels;
- the strongest disease-facing evidence comes from transfer into clinically annotated cohorts;
- Halfvarson provides the clearest stool-based CD/UC support;
- HMP2 provides clinically rich but ecologically non-identical support;
- Gevers helps define the Crohn-focused portability boundary.

The Discussion should no longer spend real estate on IBS or autoimmune disease beyond a brief note that broader applications remain exploratory.

### Conclusion

The conclusion should close on this narrower but stronger claim:

- AGP-derived dominant-taxon states provide an interpretable externalizable representation of known heterogeneous IBD dysbiosis, with strongest current support in clinically annotated external IBD cohorts.

## Figure and table guidance for this pass

The goal of this pass is a narrative rewrite, not a full figure rebuild. Therefore:

- keep the current concept figure;
- retain existing overview and external validation figures if they can support the rewritten narrative with revised captions and framing;
- avoid introducing new figure-building work unless a figure directly blocks the new story;
- treat rebuilt-cohort visuals as secondary if the main inferential story is AGP-derived transfer.

If later needed, a future pass can create a cleaner biology-facing results figure centered on:

- Halfvarson CD,
- Halfvarson UC,
- HMP2 pooled IBD,
- Gevers Crohn.

## Claims the revised manuscript should avoid

Avoid claiming:

- that AGP self-reported IBD provides diagnosis-grade evidence;
- that the paper discovers fundamentally new IBD microbial biology;
- that exact label portability is universal across cohorts;
- that IBS and autoimmune analyses are co-equal with the external IBD story.

## Claims the revised manuscript can defend

The revised manuscript can defend:

- AGP is effective for learning a phenotype-agnostic gut dominance hierarchy;
- AGP-derived states remain biologically informative in external IBD cohorts;
- transferred states capture meaningful CD-, UC-, and pooled-IBD-associated dominance structure where cohort design and ecology allow;
- deterministic dominance-lineage labels provide a readable representation of heterogeneous gut dysbiosis.

## Immediate implementation consequence

The live manuscript should be rewritten so that:

- AGP self-reported IBD results become contextual and exploratory;
- external IBD transfer results become the main Results centerpiece;
- CD and UC are split explicitly whenever the external cohort allows it;
- IBS and autoimmune analyses move out of the main Results narrative.
