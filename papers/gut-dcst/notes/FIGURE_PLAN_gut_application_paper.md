# Figure Plan: Gut DCST Application Paper

Date: 2026-03-25

## Figure Strategy

The main paper figures should do four jobs:

1. explain the DCST idea to a mathematician,
2. show the full-cohort gut landscape,
3. show which disease associations are robust,
4. show why depth-2 refinement and cross-taxonomy comparison matter.

Use existing analysis outputs wherever they are already strong, but prefer
redrawing panels in a consistent style for the manuscript.

## Main Figures

### Figure 1: Conceptual DCST schematic

Purpose:

- explain L-infinity normalization
- explain dominant-feature assignment
- explain depth-1 vs depth-2 refinement

Status:

- new figure required

Source:

- no current image file; should be drawn specifically for the paper

Why it matters:

- this is the most important figure for mathematician-friendly onboarding
- it reduces the need for microbiome-specific prior knowledge

Notes:

- keep it original
- do not reuse a literature figure here

### Figure 2: Full-cohort gut DCST landscape

Recommended panels:

- panel A: depth-1 DCST size distribution
- panel B: dominance-strength distribution

Existing assets:

- `~/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_dcst_sizes.png`
- `~/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_dominance_strength.png`
- older versions also exist in `linf/`

Status:

- existing panels available
- recommended action: redraw or re-export in manuscript style

Message:

- the gut is not dominated by one universal type
- a few large DCSTs coexist with a meaningful tail of rare dominance

### Figure 3: Adjusted association overview

Recommended panels:

- panel A: adjusted heatmap across depth-1 DCSTs and conditions
- panel B: forest plot of top adjusted associations

Existing assets:

- `~/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_heatmap.png`
- `~/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_forest_plot.png`

Status:

- existing panels available
- likely usable after style cleanup

Message:

- the main claim is not "many p-values"
- the main claim is that a nontrivial set of associations survives adjustment

### Figure 4: Contamination-aware robustness

Recommended design:

- panel A: robust vs lost associations after contaminant exclusion
- panel B: highlight examples that disappear vs examples that persist

Existing inputs:

- `~/current_projects/gut_microbiome/outputs/dcst_analysis/full_cohort_adjusted_results.csv`
- `~/current_projects/gut_microbiome/outputs/dcst_analysis/sensitivity_clean_adjusted_results.csv`

Status:

- new figure required

Message:

- the analysis is not blindly trusting all gut-dominant taxa
- contamination-aware filtering changes some signals but not all of them

### Figure 5: Depth-2 refinement case study

Recommended panels:

- panel A: depth-2 size distribution or subtype overview
- panel B: one worked case study, likely Bacteroides subtypes and IBD
- panel C: one additional example if space allows

Existing assets:

- `~/current_projects/gut_microbiome/outputs/dcst_analysis/dcst_depth2_distribution.png`
- `~/current_projects/gut_microbiome/outputs/dcst_analysis/depth2_adjusted_results.csv`

Status:

- one panel exists
- the key interpretation panel should be built new from the CSV outputs

Message:

- depth-2 is not decorative
- subdominant structure can change disease interpretation materially

### Figure 6: Cross-taxonomy robustness

Recommended panels:

- panel A: GG2 depth-1 size distribution
- panel B: SILVA vs GG2 concordance heatmap

Existing assets:

- `~/current_projects/gut_microbiome/outputs/dcst_analysis/gg2_dcst_sizes.png`
- `~/current_projects/gut_microbiome/outputs/dcst_analysis/silva_vs_gg2_concordance_heatmap.png`

Status:

- existing panels available

Message:

- the dominant-typing story is robust even when taxonomic naming changes
- the Escherichia-Shigella/GG2 gap should be shown, not hidden

## Main Tables

### Table 1: Top depth-1 DCSTs

Inputs:

- full-cohort assignment/size files

Purpose:

- orient the reader to the ecosystem

### Table 2: Headline adjusted associations

Inputs:

- `full_cohort_adjusted_results.csv`

Purpose:

- show the core findings cleanly
- include OR, CI, q-value, and whether robust to contaminant exclusion

### Table 3: Surviving univariate + adjusted overlaps

Inputs:

- `univariate_vs_adjusted_comparison.csv`

Purpose:

- show what adjustment preserves versus attenuates

### Table 4: Per-condition SILVA vs GG2 summary

Inputs:

- existing report summary / `gg2_univariate_results.csv`

Purpose:

- summarize the taxonomy comparison without overloading the reader

## Supplementary Material

Good supplementary tables/figures:

- full adjusted result table
- full depth-2 result table
- contaminant exclusion definitions
- expanded SILVA vs GG2 crosstab
- logistic model details and complete-case counts

## Figure Style Rules

1. Use a single consistent manuscript style across all figures.
2. Keep taxon labels human-readable in the final paper.
3. Avoid overwhelming heatmaps with long taxonomy strings.
4. Prefer a few carefully annotated figures over many dense ones.
5. Explain every figure so a mathematician can read it without microbiome
   field conventions.

## Literature-Derived Figures

The paper can cite literature figures for context, but should assume that most
paper figures will be original.

Use literature figures only when:

- the license clearly permits reuse, or
- permission is obtained.

Otherwise:

- redraw the idea
- cite the source paper
- adapt the visual language to the DCST narrative

## Immediate Build Order

1. Figure 1 conceptual schematic
2. Figure 2 full-cohort landscape
3. Figure 3 adjusted association overview
4. Figure 4 contamination-aware robustness
5. Figure 5 depth-2 case study
6. Figure 6 cross-taxonomy robustness

## What I Should Review Closely

- whether Figure 3 headlines the right associations
- whether Figure 4 makes the contamination story honest and not defensive
- whether Figure 5 shows genuine subtype insight rather than decorative
  complexity
- whether Figure 6 explains concordance without implying perfect taxonomic
  equivalence
