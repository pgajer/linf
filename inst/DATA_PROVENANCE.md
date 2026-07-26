# Bundled data provenance and upstream terms

The package code is distributed under the package's MIT license. The bundled
microbiome data are derived from the upstream resources described below. The
package license does not replace or narrow applicable upstream attribution
terms.

## VALENCIA-derived objects

The `valencia2k`, `valencia_linf_hypercube_1k`,
`valencia13k_dcst_depth2_merged`, and `valencia13k_dcst_depth3_merged` objects
are derived from the VALENCIA training data:

- Source repository: <https://github.com/ravel-lab/VALENCIA>
- The primary article states that the training dataset is available in that
  repository and applies the Creative Commons CC0 public-domain dedication to
  data made available with the article unless otherwise stated.
- Upstream software repository license: MIT, copyright 2023 ravel-lab
- Primary reference: France et al. (2020), "VALENCIA: a nearest centroid
  classification method for vaginal microbial communities based on
  composition," <https://doi.org/10.1186/s40168-020-00934-6>

The source-preparation scripts in `data-raw/` record the deterministic
subsampling and dCST construction parameters. Original sample identifiers are
not retained in the merged 13k assignment objects.

## American Gut Project demonstration object

The `agp_gut` object is derived from public American Gut Project 16S records
and deidentified metadata:

- ENA/BioProject accession: PRJEB11419
  (<https://www.ebi.ac.uk/ena/browser/view/PRJEB11419>)
- Primary reference: McDonald et al. (2018), "American Gut: an Open Platform
  for Citizen Science Microbiome Research,"
  <https://doi.org/10.1128/mSystems.00031-18>
- The primary article is distributed under CC BY 4.0. It states that the
  project protocol permits public deposition of non-identifying data and that
  sequence data and deidentified participant responses are placed in the
  public domain without access restrictions.
- Processing source: PRIME SILVA species-level outputs dated 2026-03-24.

The full analysis first selected at most 5,000 samples with library size at
least 1,000 using random seed 42. `data-raw/create_agp_gut_subset.py` then
creates the 766-sample package object using an explicit deterministic rule:

1. include every sample assigned to `Prevotella_7`, `Pasteurellaceae`,
   `Akkermansia`, or `Staphylococcus`;
2. exclude depth-1 `Eukaryota` and `Unassigned` labels; and
3. fill the remaining slots by a seed-42 simple random sample from the
   eligible background.

Phenotype fields are joined only after membership is fixed and cannot
influence selection. Exact run-ID membership, selection reasons, dCST labels,
and derived phenotype fields are retained in
`inst/extdata/agp_gut_meta.csv`. The script also reconstructs the bundled
count matrix from the full abundance table and retained 5,000-sample dCST
assignments. `data-raw/build_agp_gut.R` then constructs the package object.

Because inclusion probabilities differ by dCST, `agp_gut` is suitable for
demonstrating data alignment, filtering, normalization, and dCST construction
only. Its phenotype frequencies, effect sizes, and p-values must not be
interpreted as population estimates.
