# CT Results Subsection Map for the Combined dCST Methods Paper

## Purpose

This note defines the intended internal structure of the CT Results section in
the combined methods manuscript.

Primary upstream source:

- `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

## Role of the CT application in the new paper

The CT application should carry the paper's:

- independent-domain generalization
- stronger biological structure
- functional interpretability layer

It is the part of the combined paper that makes the framework feel broader than
a gut-only transfer method.

## Main-text subsection flow

## 3.3.1 CT cohort logic and outcome-association setup

### Keep

- the distinction between the labeled CT cohort and the broader CT+reference cohort
- explicit statement that the broader cohort is used for structure learning
- explicit statement that outcome tests remain restricted to the labeled CT subset

### Purpose

- this mirrors the gut structure-learning versus validation separation, but in a different biological setting

## 3.3.2 Lower-resolution VOG-cluster result

### Keep

- one concise description of the mixed persistence-associated neighborhood
- why the lower-resolution layer helps interpretation

### Keep brief

- this should orient the reader, not dominate the CT section

## 3.3.3 Full-VOG primary result

### Keep

- the high-resolution separation of clearance-oriented and persistence-oriented branches
- the strongest branch-level outcome associations
- the reason the full-VOG layer is the primary CT analysis

### Main biological anchors

- the \taxon{Lactobacillus crispatus}/\taxon{Lactobacillus jensenii} clearance-oriented branch
- persistence-oriented \taxon{L. iners} and \taxon{Ca. Lachnocurva vaginae}-linked structure

## 3.3.4 Functional follow-up

### Keep

- taxonomic interpretation where useful
- product annotation highlights
- KEGG-linked functional summaries

### Purpose

- this is the strongest “why the framework is biologically useful” layer in the combined paper

## 3.3.5 CT takeaway paragraph

This should close the CT Results section and set up the cross-application synthesis.

### It should state

- dCSTs are useful not only for taxonomic stool-cohort transfer
- the framework can also organize high-dimensional vaginal metagenomic outcome structure
- the method supports interpretable follow-up at a functional level, not only at the state-label level

## What to leave out of main text

- exhaustive branch-by-branch enumeration
- long CT-specific biological storytelling that belongs in a future standalone CT paper
- excessive cap-sensitivity or robustness-detail narrative

## Candidate paragraph import order

When migrating prose, prioritize this order:

1. CT cohort and analysis-layer setup
2. VOG-cluster orienting result
3. full-VOG primary outcome-associated result
4. functional follow-up paragraph set
5. a newly written transition paragraph into cross-application synthesis
