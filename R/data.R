#' Valencia 2k vaginal microbiome dataset
#'
#' A stratified subsample of 2,000 vaginal samples from the Valencia 13k
#' CST-classifier training set (France et al. 2020), bundled as an example
#' dataset for demonstrating L-infinity DCST construction.
#'
#' The subsample preserves proportional representation of all 13 Valencia
#' sub-CSTs and was drawn with \code{set.seed(42)}.
#'
#' @format A list with five components:
#' \describe{
#'   \item{rel}{Numeric matrix (2000 x 178). Compositional relative abundances;
#'     each row sums to 1. Rows are samples, columns are taxonomic features.}
#'   \item{cst}{Data frame (2000 x 3) with columns:
#'     \code{sample_id} (character),
#'     \code{Val_CST} (Valencia CST assignment: I, II, III, IV-A, IV-B, IV-C, V),
#'     \code{Val_subCST} (Valencia sub-CST assignment: I-A, I-B, II, III-A, III-B,
#'       IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V).}
#'   \item{reads}{Integer vector of length 2000. Per-sample read counts after
#'     taxonomic filtering. Use \code{sweep(valencia2k$rel, 1, valencia2k$reads, "*")}
#'     to reconstruct a count-like matrix.}
#'   \item{taxa}{Character vector of 178 taxon names (column names of \code{rel}).}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' The Valencia CST classifier (France et al. 2020) assigns vaginal microbiome
#' samples to community state types (CSTs) based on nearest-centroid
#' classification in relative-abundance space. The original training set
#' contains 12,881 samples and 178 taxonomic features after filtering.
#'
#' This 2,000-sample subsample is intended for vignette demonstrations. The
#' compositional matrix can be used directly with \code{\link{normalize.linf}}
#' and downstream DCST functions. For workflows that require count-like input,
#' reconstruct approximate counts using the \code{reads} vector.
#'
#' @references
#' France, M. T., Ma, B., Gajer, P., Brown, S., Humphrys, M. S., Holm, J. B.,
#' Waetjen, L. E., Brotman, R. M., & Ravel, J. (2020). VALENCIA: a nearest
#' centroid classification method for vaginal microbial communities based on
#' composition. \emph{Microbiome}, 8(1), 166.
#' \doi{10.1186/s40168-020-00934-6}
#'
#' @source
#' Subsampled from the Valencia CST-classifier training dataset.
#' See \code{data-raw/build_valencia2k.R} for the subsampling script.
#'
#' @examples
#' data(valencia2k)
#' dim(valencia2k$rel)           # 2000 x 178
#' table(valencia2k$cst$Val_CST) # CST distribution
#' head(valencia2k$taxa, 10)     # first 10 taxon names
#'
"valencia2k"


#' American Gut Project gut microbiome dataset
#'
#' A stratified subsample of 766 gut microbiome samples from the American Gut
#' Project (PRJEB11419, AGP-US-2015), bundled for demonstrating L-infinity
#' DCSTs in the gut ecosystem where single-species dominance is rare and often
#' associated with disease.
#'
#' The subsample is stratified to enrich rare disease-associated DCSTs
#' (Prevotella_7, Pasteurellaceae, Akkermansia, Staphylococcus) and
#' over-sample IBD, IBS, obesity, and cardiovascular disease cases.
#' Drawn with \code{set.seed(42)} from a 5,000-sample exploratory analysis.
#'
#' @format A list with four components:
#' \describe{
#'   \item{counts}{Integer matrix (766 x 314). Raw 16S V4 read counts.
#'     Rows are samples, columns are SILVA species-level taxa.}
#'   \item{meta}{Data frame (766 rows) with columns:
#'     \code{Run} (SRA run accession),
#'     \code{dcst_depth1}, \code{dcst_depth2} (pre-computed DCST labels),
#'     \code{IBS}, \code{IBD}, \code{Diabetes}, \code{Autoimmune},
#'     \code{Seasonal_allergies}, \code{Migraine}, \code{Acid_reflux},
#'     \code{Lung_disease}, \code{Cardiovascular_disease}, \code{Skin_condition},
#'     \code{Obesity} (binary disease indicators from self-reported AGP metadata),
#'     \code{BMI} (numeric, self-reported).}
#'   \item{taxa}{Character vector of 314 SILVA taxonomy strings.}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' The American Gut Project is a large citizen-science 16S rRNA survey of the
#' human microbiome. Health conditions are self-reported via questionnaire and
#' should be interpreted with appropriate caution.
#'
#' The count matrix can be used directly with \code{\link{filter.asv}},
#' \code{\link{normalize.linf}}, and downstream DCST functions.
#'
#' Note on Escherichia-Shigella: this genus is inflated in 16S V4 data due to
#' primer cross-reactivity and should be interpreted with caution.
#'
#' @references
#' McDonald, D., Hyde, E., Debelius, J. W., et al. (2018). American Gut: an
#' Open Platform for Citizen Science Microbiome Research. \emph{mSystems},
#' 3(3), e00031-18. \doi{10.1128/mSystems.00031-18}
#'
#' @source
#' Stratified subsample from AGP via PRIME pipeline (SILVA taxonomy).
#' See \code{data-raw/build_agp_gut.R} for the subsampling script.
#'
#' @examples
#' data(agp_gut)
#' dim(agp_gut$counts)                     # 766 x 314
#' table(agp_gut$meta$dcst_depth1)         # DCST distribution
#' sum(agp_gut$meta$IBS)                   # IBS cases
#'
"agp_gut"
