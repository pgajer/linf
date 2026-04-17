# Action Plan: Biology-First Figure Sequence

Date: 2026-04-17

## Goal

Refactor the gut-dCST manuscript so it reads as an IBD/gut-microbiome biology
paper that uses dCSTs as an analysis tool, rather than as a dCST methods paper
with an IBD example.

## Rationale

The previous first main-text figure was a dCST construction schematic. That
made the paper's visual hierarchy method-first. The paper should instead show
the reader the biological result first: IBD is the clearest disease-associated
dominance-state signal, it resolves into interpretable deeper lineages, and it
has partial external support.

## Figure Plan

1. Figure 1: Biology-first IBD overview.
   - Panel A: depthwise Cramer's V for IBD, autoimmune disease, and IBS.
   - Panel B: representative adjusted odds ratios for common IBD lineages and
     sparse enterobacterial states.
   - Panel C: external validation contrast between rebuilt-cohort validation
     and direct AGP-derived label transfer.

2. Figure 2: IBD absorb-lineage refinement.
   - Use the existing IBD lineage/refinement figure content, but make it the
     detailed second figure rather than the first IBD result encountered after a
     method schematic.

3. Figure 3: External validation and label portability.
   - Use the current validation comparison figure content, renumbered and
     captioned as the main validation result.

4. Figure 4: dCST construction schematic.
   - Keep the current dCST construction schematic in the main text, but move it
     out of the Introduction and into Materials and Methods.

## Manuscript Plan

1. Remove the dCST construction figure from the Introduction.
2. End the Introduction with the biological question and study scope.
3. Reorder the body as Introduction, Results, Materials and Methods,
   Discussion, Conclusion, and back matter.
4. Insert Figures 1-3 in the Results before any method-construction figure.
5. Insert Figure 4 in Materials and Methods after the dCST construction
   paragraph.
6. Keep claims calibrated: Figure 1 should show strongest support for IBD,
   weaker related IBS/autoimmune branches, and external validation that supports
   dominance-pattern reproducibility more strongly than universal label
   portability.

## Verification

1. Rebuild all active figure assets from source scripts.
2. Rebuild the manuscript PDF.
3. Check that figure captions and in-text figure references are consistent.
4. Check the LaTeX log for undefined references/citations.
5. Render the relevant PDF pages and visually inspect figure order, margins, and
   readability.
