# Section Map: Combined Gut+CT dCST Methods Paper

## Purpose

This document defines the intended main-text architecture before prose migration begins.

It is meant to answer:

- what sections belong in the new paper
- what each section must accomplish
- which source manuscript contributes to each section

## Main text target length logic

The paper should read as:

- one methods-forward microbiome paper

with:

- two major applications

not as:

- a theory paper plus two appendices

or:

- two application papers stapled together

## Recommended main-text section map

## Title

Working function:

- state that this is a deterministic / interpretable microbiome-state framework paper
- signal that there are multiple applications or domains

Working examples:

- Dominant-Community State Types as a Deterministic and Interpretable Framework for Microbiome State Analysis
- Deterministic Dominant-Community State Types Support Interpretable Microbiome State Analysis Across Gut and Vaginal Applications

## Abstract

Write last.

Required components:

1. problem with clustering-defined microbiome state labels
2. short dCST solution
3. gut transfer application
4. CT functional/outcome application
5. overall representation-layer conclusion

## 1. Introduction

Write late, after Results and Discussion are stable.

### Goal

Motivate the need for deterministic sample-level microbiome state labels and justify why two applications are used.

### Subsections

#### 1.1 Why microbiome studies use state labels

Cover:

- CSTs
- enterotypes
- usefulness of coarse state labels

#### 1.2 Why clustering-defined labels are limited

Cover:

- dependence on distance metric
- dependence on clustering model
- dependence on sample set
- difficulty of comparing labels across studies

#### 1.3 dCST idea in one concise paragraph

Cover:

- rank-based deterministic assignment
- absorb-policy retained states
- interpretability

#### 1.4 Why two applications are needed

Cover:

- gut as portability / external-transfer application
- CT as independent biological / functional application

### Primary source inputs

- gut paper for motivation
- Linf paper for the concise deterministic/stability framing
- CT paper only for the second-application rationale

## 2. Methods

This is the backbone of the new paper.

## 2.1 dCST framework

### Must accomplish

- define retained-feature ranking
- define depthwise dominance-lineages
- define absorb policy
- explain what makes the labels deterministic

### Must not become

- a full Linf theory section

### Source inputs

- gut methods core
- selected Linf conceptual statements

## 2.2 Shared inferential framework

### Must accomplish

- explain association testing at label/panel level
- explain multiple-testing correction
- explain subject-aware sensitivity where applicable
- explain that downstream inference can differ by application

### Source inputs

- gut paper
- CT paper

## 2.3 Gut application design

### Must accomplish

- define AGP local-QZA SILVA structure-learning cohort
- define external diagnosis-grade validation cohorts
- define transfer and rebuilt modes
- define the role of subject-level sensitivity

### Source inputs

- gut paper

## 2.4 CT application design

### Must accomplish

- define labeled CT cohort
- define CT+reference cohort
- define VOG-cluster and full-VOG layers
- define the follow-up annotation layer

### Source inputs

- CT paper

## 3. Results

This should be the longest section, but it must feel structured and cumulative.

## 3.1 Framework overview

### Must accomplish

- introduce the conceptual figure
- explain what a dCST label is in concrete microbiome terms
- establish the framework before applications begin

### Figure support

- conceptual framework figure

### Source inputs

- gut Figure 1 logic
- minimal Linf-derived conceptual framing

## 3.2 Gut application: portability across cohorts

### Must accomplish

- show that the framework can learn interpretable taxonomic state structure in a large ecological background
- show that external transfer is possible
- show where it weakens under subject-aware scrutiny

### Recommended internal flow

1. AGP hierarchy and exploratory prioritization
2. external transfer to Halfvarson
3. HMP2 and Gevers as boundary cases
4. one compact within-disease extension

### Must avoid

- overloading the paper with every gut sensitivity branch

### Source inputs

- gut paper main results

## 3.3 CT application: outcome-associated and function-linked states

### Must accomplish

- show that the same framework is informative in an independent microbiome niche
- show that dCST states can be linked to biologically interpretable outcome structure
- show that a function-rich data layer can support meaningful follow-up

### Recommended internal flow

1. VOG-cluster lower-resolution result
2. full-VOG primary result
3. functional follow-up of clearance and persistence branches

### Source inputs

- CT paper main results

## 3.4 Cross-application synthesis

### Must accomplish

- make the gut and CT results explicitly talk to each other
- answer why both applications are in one paper

### Recommended content

- gut demonstrates portability and limits
- CT demonstrates biological and functional interpretability
- both demonstrate deterministic sample-level state labeling without clustering instability

### Source inputs

- newly written synthesis only

## 4. Discussion

## 4.1 What the framework contributes

Cover:

- deterministic assignment
- interpretable labels
- sample-set stability
- reusability across domains

## 4.2 What the gut application teaches

Cover:

- transfer can work
- portability has real limits
- subject-aware attenuation matters

## 4.3 What the CT application teaches

Cover:

- framework can recover outcome-associated structure in a different niche
- function-linked follow-up strengthens biological interpretability

## 4.4 Limits and scope

Cover:

- dCSTs are representation tools, not causal claims
- downstream signal quality depends on cohort design and feature universe
- not every application will support equally strong transfer

## 4.5 Future directions

Cover:

- multi-omics integration
- longitudinal/subject-aware extensions
- richer clinical validation

## 5. Data and code availability

Keep concise and current.

## Figure map

## Figure 1

Conceptual framework figure.

### Function

- define dCSTs visually

## Figure 2

Gut application figure.

### Function

- summarize AGP-derived taxonomic state structure and external transfer

## Figure 3

CT application figure.

### Function

- summarize CT-associated VOG/VOG-cluster state structure

## Figure 4

CT biological follow-up figure.

### Function

- show taxonomic / functional interpretation

## Figure 5

Cross-application synthesis figure.

### Function

- visually unify the manuscript

## Table map

## Table 1

High-level gut results summary.

## Table 2

High-level CT results summary.

## Optional Table 3

Framework property or application-comparison table if needed.

## Supplement map

## Supplement should contain

- additional gut exploratory phenotype screens
- additional gut sensitivity tables
- additional transfer audit details
- additional CT robustness outputs
- longer annotation tables
- any small Linf-derived technical note if still useful

## Supplement should not contain

- whole duplicated application narratives
- theory sections that the main paper does not need

## Section-by-source matrix

## Sections mostly from gut paper

- 1.1 to 1.2 motivation backbone
- 2.3 gut application design
- 3.2 gut results
- 4.2 gut-specific discussion

## Sections mostly from CT paper

- 2.4 CT application design
- 3.3 CT results
- 4.3 CT-specific discussion

## Sections lightly informed by Linf paper

- 1.3 concise conceptual framing
- 2.1 shared framework definition

## Sections that must be newly written

- 1.4 two-application rationale
- 3.4 cross-application synthesis
- 4.1 overall framework contribution
- 4.4 limits and scope
- 4.5 future directions

## Editorial guardrails

## Guardrail 1

Do not let the gut application dominate the new paper simply because its source text is already in LaTeX.

## Guardrail 2

Do not let the CT application dominate the new paper simply because it has richer biology.

## Guardrail 3

Do not let the Linf theory expand beyond what a methods-facing microbiome journal actually needs.

## Guardrail 4

Every major section should answer:

- what does this show about the framework?

not just:

- what happened in this one cohort?

## Recommended immediate implementation order

1. create a figure inventory note
2. create a source-to-section extraction note if needed
3. create the new manuscript skeleton `.tex`
4. write Methods first
5. migrate Results second
6. write Discussion third
7. write Introduction and Abstract last
