# Action Plan: IBD Reviewer Strengthening Pass

Date: 2026-04-17

Scope: address the current IBD-expert reviewer recommendations for the gut-dCST manuscript while keeping the paper submission-facing and avoiding unsupported claim expansion.

Canonical manuscript sources:

- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex`
- `/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper_supplement.tex`

## Reviewer Concerns To Address

1. Repeated samples in external validation could inflate apparent evidence if sample-level tests are interpreted as subject-level replication.
2. Direct AGP-derived label transfer can be affected by case/control differences in mapping rates, especially in Halfvarson and at deeper depths.
3. AGP antibiotic metadata should be used where available rather than leaving antibiotic exposure only as an untested limitation.
4. Rare \emph{Morganella} and \emph{Proteus} dominance states should remain exploratory; the manuscript should not imply severity or disease activity without supporting metadata.
5. The inferential hierarchy should be explicit: BH-adjusted frequentist analyses support primary claims; Bayesian models stabilize sparse-label estimates but do not replace multiplicity-controlled discovery.
6. Absorb-policy robustness should be shown directly rather than assumed.
7. The phrase "species-level table" needs qualification because retained labels mix genus, family, order, species, and unresolved features.
8. External validation tables should make repeated-subject structure and transfer-mapping caveats easier for reviewers to see.

## Execution Plan

### Stage 1. Build reviewer-facing sensitivity assets

Create a reproducible builder that generates manuscript-facing tables for:

- AGP IBD antibiotic-status sensitivity using age/sex/BMI complete cases, plus coarse antibiotic-status adjustment, and exclusion of recent antibiotic users.
- External validation first-sample-per-subject sensitivity for cohorts with repeated subjects.
- Case/control mapping-balance summaries for direct AGP-derived label transfer.
- Absorb versus pure-policy robustness for selected IBD labels.
- Taxonomic rank transparency for selected IBD label components.

Done criterion:

- Tables are generated under `papers/gut-dcst/assets/tables/` and can be regenerated from one script.

### Stage 2. Recalibrate the main-text IBD claims

Revise the IBD Results and Discussion so that:

- Common *Bacteroides*-centered dCST dominance-lineages are presented as the main reproducible result.
- Sparse *Morganella* and *Proteus* states are described as exploratory facultative-anaerobe/Enterobacterales examples, not as markers of severity.
- Antibiotic-status sensitivity is added as supporting evidence for the main AGP IBD signal, while still acknowledging residual medication and disease-activity confounding.
- The distinction between sample-level validation and subject-level sensitivity is made explicit.
- Halfvarson is described as strong sample-level support with subject-level attenuation for some transferred comparisons; HMP2 is described as selective but subject-level-stable transfer support.

Done criterion:

- A reviewer can identify exactly which IBD claims are strong, which are selective, and which are exploratory.

### Stage 3. Strengthen Methods and Supplement

Revise Methods and Supplement so that:

- Antibiotic-status sensitivity, subject-level validation sensitivity, mapping balance, absorb-policy robustness, and taxonomic-rank transparency are documented.
- Bayesian models are explained as shrinkage/stabilization for sparse labels, with BH-adjusted frequentist results remaining the primary multiple-testing framework.
- External validation modes are explained without notebook-style process narration.
- New supplementary tables S9-S13 summarize the sensitivity assets.

Done criterion:

- The main paper is calibrated, and the supplement contains enough detail for a reviewer to audit the calibration without reading pipeline logs.

### Stage 4. Rebuild and verify

Regenerate manuscript PDFs and check for:

- Successful main and supplement LaTeX builds.
- No undefined references or build-stamp failures.
- No remaining "severe dysbiosis" or "frozen transfer" terminology.
- Clean diff and no accidental `.DS_Store` or temporary PDF-render artifacts.

Done criterion:

- Updated PDFs reflect the source edits and verification commands pass or any residual warnings are documented.

## Execution Status

- Stage 1: completed. Added `build_ibd_reviewer_sensitivity_assets.R` and generated Supplementary Tables S9-S13 under `papers/gut-dcst/assets/tables/`.
- Stage 2: completed. Recalibrated the main-text IBD claims around common *Bacteroides*-centered dCST dominance-lineages, exploratory sparse enterobacterial states, antibiotic-status robustness, and subject-level external-validation sensitivity.
- Stage 3: completed. Added supplement methods language and compact Supplementary Tables S9-S13 for antibiotic sensitivity, one-sample-per-subject validation, mapping balance, absorb-policy robustness, and retained-feature taxonomic ranks.
- Stage 4: completed. Rebuilt the main manuscript and supplement PDFs, checked LaTeX logs, searched the rendered PDFs for stale terminology, and visually inspected the updated validation and supplement pages.

## Key Outcomes From Execution

- The leading AGP IBD dCST dominance-lineages were stable after adding coarse antibiotic status and after excluding samples reporting antibiotic use in the previous week, month, or six months.
- Halfvarson remains the strongest sample-level external validation cohort, but one-sample-per-subject sensitivity attenuates the pooled IBD and Crohn direct-transfer comparisons; the UC direct-transfer signal remains significant.
- HMP2 provides selective but useful AGP-derived label-transfer support: the depth-2 *Bacteroides*/*Faecalibacterium* signal remains significant in the subject-level sensitivity analysis.
- Direct AGP-derived label transfer has phenotype-associated mapping imbalance in Halfvarson at deeper depths, so the manuscript now treats label portability as selected and calibrated rather than universal.
- Sparse *Morganella* and *Proteus* states remain exploratory; disease severity or activity claims were removed because those metadata are not available in AGP.
- The manuscript now states that retained dCST labels mix taxonomic ranks and should be interpreted as retained-feature dominance labels rather than species-only labels.
