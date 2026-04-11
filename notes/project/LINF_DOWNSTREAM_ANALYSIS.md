# LINF Downstream Analysis Pipeline

This note records the current downstream analysis pattern used with the
`linf` package so that future projects can reuse it without reconstructing the
gut DCST paper workflow from scattered scripts.

The intended split is:

- `linf`: core methods package and reusable analysis logic
- `gut_microbiome`: project workspace holding discovery inputs, cohort metadata,
  derived outputs, and runnable analysis scripts

## Core Idea

The `linf` package defines the DCST hierarchy itself:

- `filter.asv()` for library-size and prevalence filtering
- `normalize.linf()` for within-sample rank-preserving normalization
- `linf.csts()` for depth-1 DCST assignment
- `refine.linf.csts()` for recursive depth-2+ refinement

The downstream pipeline begins after DCST labels exist. In practice that means
the package gives us the hierarchy, while project scripts add:

- phenotype extraction
- omnibus screening
- label-level effect estimation
- covariate adjustment
- sensitivity analyses
- taxonomy-concordance checks
- external validation
- figure and table generation

## Reference Implementations

Current gut-project reference scripts:

- rare-policy AGP discovery: [`/Users/pgajer/current_projects/gut_microbiome/scripts/run_agp_dcst_analysis.py`](/Users/pgajer/current_projects/gut_microbiome/scripts/run_agp_dcst_analysis.py)
- absorb-policy AGP discovery depth scan: [`/Users/pgajer/current_projects/gut_microbiome/scripts/run_agp_absorb_depthscan.R`](/Users/pgajer/current_projects/gut_microbiome/scripts/run_agp_absorb_depthscan.R)
- absorb label-level frequentist + Bayesian follow-up: [`/Users/pgajer/current_projects/gut_microbiome/scripts/run_absorb_label_followup.R`](/Users/pgajer/current_projects/gut_microbiome/scripts/run_absorb_label_followup.R)

The important convention is that scripts live under `gut_microbiome/scripts/`;
`gut_microbiome/outputs/` is for products, not source.

## Pipeline Stages

### 1. Discovery Inputs

Typical required inputs are:

- abundance matrix
- sample metadata
- project/run inclusion rule
- phenotype definitions

In the gut project, the main discovery matrix is the PRIME-derived SILVA
species abundance table, with AGP restricted by `BioProject == PRJEB11419`.

### 2. Preprocessing

The current default preprocessing pattern is:

- drop samples with total library size below `1000`
- retain taxa with count `>= 2` in at least `5%` of retained samples
- then normalize with `normalize.linf()`

These values should be treated as project defaults, not immutable package
constants.

### 3. DCST Construction

Two low-frequency policies are currently important:

- `rare`: low-support dominant states collapse into a `RARE_DOMINANT` bucket
- `absorb`: low-support states are absorbed into the nearest retained parent

For deeper analyses, absorb-mode hierarchies must be built in absorb mode from
the start. A rare-policy hierarchy cannot later be “converted” into a true
absorb hierarchy.

Typical depth schedule for absorb scans:

- depth 1: `n0 = 50`
- depth 2: `25`
- depth 3: `12`
- depth 4: `10`
- exploratory depths 5 and 6: `8`, `5`

### 4. Omnibus Screening

Before drilling into individual labels, phenotype-by-depth panels are screened
with an omnibus contingency-table test.

Current absorb reference implementation:

- Monte Carlo Fisher’s exact test
- Cramer’s V as the effect-size summary
- optional adaptive simulation counts:
  - exploratory pass at `B = 10000`
  - confirmation rerun at `B = 100000` for panels whose exploratory q-value is
    below a trigger threshold

This stage answers:

- which phenotypes show any DCST-structure association at a given depth
- whether deeper levels add information or merely fragment the cohort

### 5. Label-Level Frequentist Follow-up

For phenotype-depth panels worth following up, the current frequentist layer is:

- unadjusted `2x2` Fisher tests for single labels against the rest
- Benjamini-Hochberg correction
- covariate-adjusted logistic regression for supported labels

The gut-paper style adjusted model is:

`outcome ~ label_present + Host_Age + Host_Sex + Host_BMI`

Practical gating is important. We do not fit adjusted models to every tiny deep
label; instead we require minimum label support and minimum case/control counts.

### 6. Bayesian Follow-up Layer

This is the major recent extension to the downstream pipeline.

The current reference Bayesian layer has two parts:

- unadjusted prevalence comparison with a Jeffreys-prior Beta-Binomial model
- adjusted logistic regression with `arm::bayesglm()`

This layer is not intended to replace the frequentist screen. Instead it is a
parallel estimation layer that is especially useful when:

- deeper labels are sparse
- separation makes ordinary logistic fits unstable
- we want graded evidence rather than only p-values

Current interpretation pattern:

- keep frequentist omnibus q-values as the main discovery gate
- keep frequentist adjusted models for comparability with the existing paper
- use Bayesian adjusted estimates to stabilize deeper or smaller labels and to
  identify plausible follow-up candidates that frequentist screening may treat
  conservatively

### 7. Structural and Robustness Outputs

Reusable outputs worth preserving for each discovery run:

- per-depth label counts
- reassignment summaries
- depth-transition crosstabs and heatmaps
- omnibus summary tables and heatmaps
- label-level result tables
- top-hit tables
- human-readable Markdown summaries

This makes each run self-describing and lets manuscripts pull from stable files
rather than hidden in-memory objects.

### 8. Validation and Cross-System Checks

Once discovery is stable, the current project pattern adds:

- taxonomy-system concordance checks, e.g. SILVA versus GG2
- contamination-aware sensitivity reruns
- external validation on compatible cohorts

The key rule is that validation should use the same DCST policy as discovery.
Rare-policy validation and absorb-policy validation are different branches, not
interchangeable labels.

## Reuse Checklist For New Projects

When porting this downstream pattern to another dataset or phenotype family:

1. Define the abundance table and metadata join key.
2. Define the phenotype columns or case-control contrasts.
3. Choose low-frequency policy: `rare` or `absorb`.
4. Choose the depth schedule and minimum support thresholds.
5. Run the omnibus scan first.
6. Limit label-level modeling to phenotype-depth panels that are supported and
   interpretable.
7. Add Bayesian follow-up when deeper or sparse labels matter.
8. Save every stage as explicit tables and plots.

## Suggested Long-Term Package Direction

The current Bayesian layer lives in project scripts, not in exported `linf`
functions. That is appropriate for now, but if the same pattern proves useful
across multiple projects, the natural next step would be a small package-level
downstream module or vignette covering:

- omnibus screening helpers
- label-level frequentist follow-up helpers
- Bayesian follow-up helpers
- standard output schemas

Until then, the gut-project scripts above should be treated as the canonical
reference implementation.
