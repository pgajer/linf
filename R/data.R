#' Valencia 2k vaginal microbiome dataset
#'
#' A stratified subsample of 2,000 vaginal samples from the Valencia 13k
#' CST-classifier training set (France et al. 2020), bundled as an example
#' dataset for demonstrating dCST construction after L-infinity normalization.
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
#' and downstream dCST functions. For workflows that require count-like input,
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
#' Subsampled from the VALENCIA training data at
#' \url{https://github.com/ravel-lab/VALENCIA}. See
#' \code{data-raw/build_valencia2k.R} and the installed
#' \code{DATA_PROVENANCE.md} file for construction and licensing details.
#'
#' @examples
#' data(valencia2k)
#' dim(valencia2k$rel)           # 2000 x 178
#' table(valencia2k$cst$Val_CST) # CST distribution
#' head(valencia2k$taxa, 10)     # first 10 taxon names
#'
"valencia2k"


#' Valencia 1k four-component hypercube embedding example
#'
#' A stratified 1,000-sample subset of the Valencia 13k vaginal microbiome
#' training set, reduced to four selected phylotype coordinates and
#' L1-normalized over those coordinates. The object is bundled as a lightweight
#' example for the zero-aware homogeneous-coordinate hypercube embedding in
#' \code{\link{linf.hypercube.embedding}}.
#'
#' @format A list with four components:
#' \describe{
#'   \item{rel4}{Numeric matrix (1000 x 4). Rows are samples and columns are
#'     \code{Li}, \code{Lc}, \code{Gv}, and \code{Bv}. Each row sums to 1 after
#'     restricting the original Valencia profile to the four mapped taxa.}
#'   \item{meta}{Data frame (1000 rows) with columns:
#'     \code{sample_id} (anonymized bundled-data identifier),
#'     \code{source_row} (row number in the filtered Valencia 13k source),
#'     \code{Val_CST}, \code{Val_subCST}, \code{selected_mass} (original
#'     relative-abundance mass carried by the four selected taxa), and
#'     \code{dominant_component} (largest of \code{Li}, \code{Lc},
#'     \code{Gv}, \code{Bv} after renormalization).}
#'   \item{component_map}{Named character vector mapping \code{Li},
#'     \code{Lc}, \code{Gv}, and \code{Bv} to the original Valencia taxon names:
#'     \code{Lactobacillus_iners}, \code{Lactobacillus_crispatus},
#'     \code{Gardnerella_vaginalis}, and \code{BVAB1}.}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' Rows with zero total mass in the four selected taxa are removed before
#' renormalization. Sampling is stratified by the dominant selected component,
#' using \code{set.seed(20261604)}. The object is not intended to replace the
#' full Valencia matrix; it is a compact reproducible example for visualizing
#' compositional projective-space coordinate charts.
#'
#' @source
#' Generated from the VALENCIA training data at
#' \url{https://github.com/ravel-lab/VALENCIA}. See
#' \code{data-raw/build_valencia_linf_hypercube_1k.R} and the installed
#' \code{DATA_PROVENANCE.md} file.
#'
#' @examples
#' data(valencia_linf_hypercube_1k)
#' dim(valencia_linf_hypercube_1k$rel4)
#' table(valencia_linf_hypercube_1k$meta$dominant_component)
#' emb <- linf.hypercube.embedding(
#'   valencia_linf_hypercube_1k$rel4,
#'   reference = "Li"
#' )
#' dim(emb)
#'
"valencia_linf_hypercube_1k"


#' Valencia 13k merged depth-2 dCST assignments
#'
#' A lightweight bundled assignment asset containing merged depth-2
#' Dominant community state type (dCST) labels for the full Valencia
#' 13k vaginal microbiome training set.
#'
#' @format A list with five components:
#' \describe{
#'   \item{assignments}{Data frame with one row per source sample and columns
#'     \code{sample_id} (anonymized bundled-data identifier), \code{source_row}
#'     (row number in the filtered Valencia 13k source), \code{Val_CST},
#'     \code{Val_subCST}, \code{dcst_depth1}, and \code{dcst_depth2}.}
#'   \item{summaries}{Named list of per-depth summary tables. Each table
#'     contains \code{depth}, \code{dcst_label}, \code{n}, \code{prop}, and
#'     \code{path_length}.}
#'   \item{feature_labels}{Character vector of source taxon labels.}
#'   \item{params}{List recording the construction parameters.}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' The asset is computed from the Valencia 13k compositional matrix after
#' L-infinity normalization. dCSTs use \code{n0 = 50} and the merged
#' \code{low.freq.policy = "absorb"} view, so samples from low-support
#' provisional states are reassigned to retained states rather than stored as
#' explicit rare buckets.
#'
#' @source
#' Generated from the VALENCIA training data at
#' \url{https://github.com/ravel-lab/VALENCIA}. See
#' \code{data-raw/build_valencia13k_merged_dcst_depths.R} and the installed
#' \code{DATA_PROVENANCE.md} file.
#'
#' @examples
#' data(valencia13k_dcst_depth2_merged)
#' nrow(valencia13k_dcst_depth2_merged$assignments)
#' head(valencia13k_dcst_depth2_merged$summaries$depth2)
#'
"valencia13k_dcst_depth2_merged"


#' Valencia 13k merged depth-3 dCST assignments
#'
#' A lightweight bundled assignment asset containing merged depth-3
#' Dominant community state type (dCST) labels for the full Valencia
#' 13k vaginal microbiome training set.
#'
#' @format A list with five components:
#' \describe{
#'   \item{assignments}{Data frame with one row per source sample and columns
#'     \code{sample_id} (anonymized bundled-data identifier), \code{source_row}
#'     (row number in the filtered Valencia 13k source), \code{Val_CST},
#'     \code{Val_subCST}, \code{dcst_depth1}, \code{dcst_depth2}, and
#'     \code{dcst_depth3}.}
#'   \item{summaries}{Named list of per-depth summary tables. Each table
#'     contains \code{depth}, \code{dcst_label}, \code{n}, \code{prop}, and
#'     \code{path_length}.}
#'   \item{feature_labels}{Character vector of source taxon labels.}
#'   \item{params}{List recording the construction parameters.}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' The depth-3 asset extends \code{\link{valencia13k_dcst_depth2_merged}} by one
#' additional hierarchical dCST refinement. All levels use \code{n0 = 50} and
#' \code{low.freq.policy = "absorb"}. This object is intended as a reusable
#' source for selecting richer VALENCIA-derived component sets without
#' recomputing the full hierarchy from the 13k source matrix.
#'
#' @source
#' Generated from the VALENCIA training data at
#' \url{https://github.com/ravel-lab/VALENCIA}. See
#' \code{data-raw/build_valencia13k_merged_dcst_depths.R} and the installed
#' \code{DATA_PROVENANCE.md} file.
#'
#' @examples
#' data(valencia13k_dcst_depth3_merged)
#' nrow(valencia13k_dcst_depth3_merged$assignments)
#' head(valencia13k_dcst_depth3_merged$summaries$depth3)
#'
"valencia13k_dcst_depth3_merged"


#' American Gut Project gut microbiome dataset
#'
#' A stratified subsample of 766 gut microbiome samples from the American Gut
#' Project (PRJEB11419, AGP-US-2015), bundled for demonstrating L-infinity
#' dCST construction in a gut ecosystem.
#'
#' The subsample includes all samples assigned to four uncommon demonstration
#' dCSTs (Prevotella_7, Pasteurellaceae, Akkermansia, and Staphylococcus).
#' The remaining slots are a simple random sample, drawn with seed 42, from
#' the eligible background after excluding Eukaryota and Unassigned labels.
#' Phenotype fields are joined only after membership is fixed and do not
#' influence selection. Because inclusion probabilities differ by dCST, this
#' object is a computational demonstration dataset rather than a probability
#' sample of the underlying cohort.
#'
#' @format A list with four components:
#' \describe{
#'   \item{counts}{Integer matrix (766 x 314). Raw 16S V4 read counts.
#'     Rows are samples, columns are SILVA species-level taxa.}
#'   \item{meta}{Data frame (766 rows) with columns:
#'     \code{Run} (SRA run accession),
#'     \code{dcst_depth1}, \code{dcst_depth2} (pre-computed dCST labels),
#'     \code{IBS}, \code{IBD}, \code{Diabetes}, \code{Autoimmune},
#'     \code{Seasonal_allergies}, \code{Migraine}, \code{Acid_reflux},
#'     \code{Lung_disease}, \code{Cardiovascular_disease}, \code{Skin_condition},
#'     \code{Obesity} (binary disease indicators from self-reported AGP metadata),
#'     \code{BMI} (numeric, self-reported), and
#'     \code{selection_reason} (target-dCST inclusion or seeded background
#'     sampling).}
#'   \item{taxa}{Character vector of 314 SILVA taxonomy strings.}
#'   \item{source}{Character string documenting provenance.}
#' }
#'
#' @details
#' The American Gut Project is a large citizen-science 16S rRNA survey of the
#' human microbiome. Health conditions are self-reported via questionnaire and
#' should be interpreted with appropriate caution.
#'
#' The phenotype fields must not be used with this dCST-stratified subset for
#' population prevalence estimates, effect-size estimation, or association
#' testing. The exact selection is generated by
#' \code{data-raw/create_agp_gut_subset.py}; run-ID membership, selection
#' reasons, and derived annotations are retained in
#' \code{inst/extdata/agp_gut_meta.csv}.
#'
#' The count matrix can be used directly with \code{\link{filter.asv}},
#' \code{\link{normalize.linf}}, and downstream dCST functions.
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
#' Derived from the public American Gut Project records under ENA accession
#' \href{https://www.ebi.ac.uk/ena/browser/view/PRJEB11419}{PRJEB11419} via the
#' PRIME pipeline. See \code{data-raw/create_agp_gut_subset.py},
#' \code{data-raw/build_agp_gut.R}, and the installed
#' \code{DATA_PROVENANCE.md} file.
#'
#' @examples
#' data(agp_gut)
#' dim(agp_gut$counts)                     # 766 x 314
#' table(agp_gut$meta$dcst_depth1)         # dCST distribution
#' sum(agp_gut$meta$IBS)                   # IBS cases
#'
"agp_gut"
