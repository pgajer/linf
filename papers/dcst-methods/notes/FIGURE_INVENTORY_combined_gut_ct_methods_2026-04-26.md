# Figure Inventory for the Combined Gut+CT dCST Methods Paper

## Purpose

This inventory lists likely main-text figures for the new combined manuscript
and records whether each should be:

- reused with light adaptation
- redesigned from existing source material
- created new

## Figure 1: Conceptual dCST framework figure

### Status

- reuse with adaptation

### Likely source

- conceptual logic from:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`

### Function in new paper

- explain retained-feature dominance
- explain depthwise dominance-lineages
- explain absorb-policy retained states
- visually distinguish dCST labels from clustering-defined CST/enterotype labels

### Needed changes

- make the examples less gut-specific
- make the caption explicitly cross-application

## Figure 2: Gut application overview

### Status

- redesign from existing source material

### Likely sources

- gut Figure 2 logic and supporting assets:
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/figures/FIGURE_2_ibd_results_overview.png`
  `/Users/pgajer/current_projects/linf/papers/gut-dcst/assets/figures/FIGURE_4_external_validation_modes.png`

### Function in new paper

- show AGP exploratory IBD prioritization
- show external transfer summary
- show subject-aware attenuation in a compact methods-paper way

### Needed changes

- simplify around the new umbrella thesis
- retain Halfvarson as the main positive external test
- keep HMP2 and Gevers as boundary cases without overloading the figure

## Figure 3: CT application overview

### Status

- redesign from CT source material

### Likely source

- main CT manuscript figures in:
  `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

### Function in new paper

- show the VOG-cluster orienting result
- show the full-VOG primary outcome-associated branches
- make the CT section visually parallel to the gut section without forcing identical graphics

### Needed changes

- standardize design language with the new manuscript
- emphasize framework output rather than CT-only local narrative

## Figure 4: CT biological follow-up figure

### Status

- redesign from CT source material

### Likely source

- taxonomy/product/KEGG follow-up figures or tables in:
  `/Users/pgajer/current_projects/CT_clearance/docs/ct_dcst_manuscript.pdf`

### Function in new paper

- demonstrate that dCST signal-bearing states can be interpreted biologically
- give the combined paper a strong functional layer

### Needed changes

- compress to the most framework-relevant biological contrasts
- avoid becoming a second CT-paper centerpiece

## Figure 5: Cross-application synthesis figure

### Status

- create new

### Function in new paper

- make the manuscript feel like one coherent methods paper
- contrast what each application contributes to the overall framework claim

### Possible designs

#### Option A

A three-column synthesis panel:

- left: deterministic dCST framework properties
- middle: gut portability / transfer proof point
- right: CT biological / functional proof point

#### Option B

A matrix-style comparison figure with rows such as:

- deterministic assignment
- portability
- subject-aware robustness
- functional interpretability
- application domain

and columns:

- gut
- CT

## Figures likely to remain supplement-only

- AGP exploratory non-IBD follow-up figures
- full gut transfer audit tables converted to figure-like summaries
- CT robustness or cap-sensitivity figures
- extended branch-level lineage figures

## Immediate next steps for figures

1. identify the exact source scripts for the current gut conceptual and overview figures
2. identify the exact CT figure source files, not only the PDF output
3. decide whether CT needs one main figure or two in the main text
4. sketch the new cross-application synthesis figure before prose migration goes too far
