# Action Plan: Reviewer A-M Revision for Gut Submission (2026-04-25)

## Goal
Sharpen the manuscript for submission to *Gut* by addressing all reviewer comments A through M, with particular emphasis on stronger clinical framing in the Abstract and Introduction, clearer interpretation of within-disease results, and reduction of figure redundancy.

## A. Crohn-location biological interpretation
- Add one sentence in the Halfvarson within-disease Results paragraph stating that the ileal-involved Crohn pattern is directionally consistent with the classical ileal Crohn / adherent-invasive *E. coli* literature.
- Add the same clinical interpretation in the Discussion, positioned as recovery of an established ileal-CD-associated signature rather than a new discovery.
- Add a supporting reference for the AIEC / ileal Crohn linkage.

## B. Rebuilt-vs-transfer asymmetry for Crohn location
- Revise Results to say explicitly that within-disease Crohn structure is present in Halfvarson, but is captured more clearly by rebuilt mode than by AGP-frozen transfer.
- Add one sentence in Discussion stating that this is a framework boundary: AGP-frozen labels appear better suited to disease-vs-control transfer than to subtype resolution when the discovery hierarchy does not anchor the relevant child state.
- Tie this interpretation to Figure 2C and the bootstrap results.

## C. Calprotectin cutoff at 250 µg/g
- Add a Methods sentence explaining that 250 µg/g was used as a pragmatic high-inflammatory threshold rather than as a remission cutpoint.
- Support that sentence with a clinical calprotectin guidance reference indicating that values above 250 µg/g are commonly interpreted as likely active inflammation.
- Adjust Results wording so the calprotectin follow-up is framed as high-versus-lower inflammatory burden, not as a strict active-vs-remission classification.

## D. Montreal definition of “ileal-involved”
- Add a Methods sentence stating that ileal-involved Crohn comprised Montreal L1, L3, L1+L4, and L3+L4, whereas colonic-only Crohn comprised L2.
- Mirror that definition briefly in the S13B caption or summary sentence so the supplement is self-explanatory.

## E. Effect sizes for within-disease findings
- Add key OR and confidence-interval values to the Halfvarson within-disease Results paragraph:
  - rebuilt Crohn-location depth-2 effect,
  - AGP-transfer Crohn-location depth-1 or depth-2 effect,
  - rebuilt calprotectin depth-3 effect,
  - AGP-transfer calprotectin depth-2 effect,
  - weak UC-extent effect as context.
- Keep the detailed per-depth table in the supplement, but make the main text clinically legible without forcing readers to hunt for S13A-S13D.

## F. HMP2 within-disease follow-up not performed
- Add an explicit sentence in Results or Methods stating that HMP2 within-disease follow-up was not emphasized because the current analyzed slice is mucosal-only, longitudinal, and small after CD/UC/site stratification.
- Add a sentence in Discussion flagging HMP2 within-disease analysis as an obvious next extension once a stool-equivalent or larger clinically resolved slice is available.

## G. Explain the mucosal-only HMP2 PRIME slice
- Add one concrete sentence in Results/Methods stating that the current analyzed PRIME-derived HMP2 slice contains rectal, ileal, and colonic mucosal samples rather than stool.
- State that this choice reflects the available local PRIME slice used in this manuscript.
- Tighten the HMP2 interpretation so readers understand that mucosal/stool differences are part of the reason this cohort is a complement rather than a stool-transfer replication.

## H. PSC-IBD overlap in Halfvarson
- Check whether PSC status is available in the harmonized Halfvarson metadata used for the paper.
- If PSC is identifiable:
  - report whether PSC-IBD subjects were retained, excluded, or absent from the primary branch,
  - add a brief limitation note if retained without stratification.
- If PSC is not resolved in the current harmonized slice:
  - say so explicitly in the manuscript,
  - note that unresolved PSC overlap remains one possible source of residual within-IBD heterogeneity.

## I. Bootstrap consensus interpretive scale
- Add a Methods sentence defining a working interpretive scale for bootstrap recurrence:
  - around or above 50% = working evidence of subject-level stability,
  - 25–50% = suggestive / borderline,
  - below 25% = not robust to repeated-sample selection.
- Recalibrate Results and Discussion language so 49.6% is described as near-threshold or strongest borderline support rather than simply “stable.”
- Recalibrate the UC transfer language similarly so 29.9% is described as the clearest transferred subject-aware signal but only suggestive by this scale.

## J. Abstract update for within-disease findings
- Add one sentence to the Abstract noting that rebuilt-mode Halfvarson follow-up also separated ileal-involved from colonic-only Crohn disease.
- Make clear that this subtype distinction was stronger in rebuilt mode than in AGP-frozen transfer.
- Keep the abstract concise and honest about the exploratory status of the within-disease analyses.

## K. Multiple-testing statement for within-disease analyses
- Add an explicit Methods sentence stating that the UC extent, Crohn location, and calprotectin analyses were exploratory and BH correction was applied within each comparison’s label panel rather than across the three follow-up analyses jointly.
- Add a short Discussion reminder that these follow-ups are hypothesis-generating.

## L. Figure 3 redundancy vs Figure 2B
- Remove the sample-level signal panel from the current Figure 3 and keep Figure 3 as the depth-coverage / mapping-coverage figure only.
- Update the figure-builder script rather than hand-editing the PNG.
- Revise the Figure 3 caption and the main-text callout so Figure 2B is the overview-with-bootstrap figure and Figure 3 is the dedicated mapping-coverage figure.

## M. Figure 2C reading guide
- Add a one-line in-panel reading guide to Figure 2C, for example:
  - `cell = depth / sample q -> first-subject q / boot %`
- Keep the caption explanation, but reduce the burden on first-pass readers by making the encoding self-contained in the figure.

## Additional clinical-framing upgrades for Gut submission
- Tighten the Abstract opening so the paper reads more like a clinically anchored IBD microbiome paper and less like a generic methods note.
- Sharpen the Introduction paragraph around the external cohorts so it explicitly says the paper tests both:
  - disease-vs-control transfer using CD and UC diagnosis-grade cohorts,
  - and exploratory within-disease structure where metadata permit.
- Ensure the Results hierarchy mirrors this framing:
  - main external IBD transfer,
  - then within-disease Halfvarson follow-up,
  - then HMP2 and Gevers boundary interpretation.

## Execution order
1. Audit manuscript, supplement, figure scripts, and metadata against A-M.
2. Add any needed bibliography entries for AIEC and calprotectin guidance.
3. Patch Abstract, Introduction, Results, Methods, and Discussion.
4. Patch supplement captions/wording where needed.
5. Refactor Figure 3 to remove redundancy and update Figure 2C guide text.
6. Rebuild both PDFs and verify the affected pages/sections.
7. Summarize how each reviewer point was addressed.
