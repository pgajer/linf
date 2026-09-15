# Run from the repository root. Optional --rebuilt=DIR and --prepared=DIR.
args <- commandArgs(TRUE)
stopifnot(all(grepl('^--(rebuilt|prepared)=', args)))
manifest <- read.csv('inst/DATA_MANIFEST.csv', stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(manifest$file), all(nchar(manifest$sha256) == 64L))
for (kind in c('asset', 'agp-prepared')) {
 rows <- manifest[manifest$kind == kind, ]
 prefix <- if (kind == 'asset') 'data' else 'inst/extdata'
 actual <- unname(tools::md5sum(file.path(prefix, rows$file)))
 if (!identical(actual, rows$md5)) {
   stop("Input/asset checksum mismatch: ", paste(rows$file[is.na(actual) | actual != rows$md5], collapse = ", "),
        "; use the original bytes (including the .gitattributes line-ending rules).")
 }
}
objects <- new.env()
for (f in list.files('data', '[.]rda$', full.names = TRUE)) load(f, objects)
v <- objects$valencia2k
stopifnot(nrow(v$cst)==2000L, !anyDuplicated(v$cst$source_row),
          identical(v$cst$sample_id, rownames(v$rel)))
for (name in c('valencia13k_dcst_depth2_merged','valencia13k_dcst_depth3_merged')) {
 asset <- objects[[name]]
 stopifnot(nzchar(asset$params$generator$version),
           identical(asset$params$generator$upstream_revision,
                     unique(manifest$upstream_revision[manifest$kind == 'upstream'])),
           !anyDuplicated(asset$assignments$source_row))
 j <- match(v$cst$source_row, asset$assignments$source_row)
 stopifnot(!anyNA(j), identical(v$cst$Val_CST,asset$assignments$Val_CST[j]),
           identical(v$cst$Val_subCST,asset$assignments$Val_subCST[j]))
 for (d in seq_along(asset$summaries)) {
  counts <- table(asset$assignments[[paste0('dcst_depth',d)]])
  summary <- asset$summaries[[d]]
  stopifnot(identical(unname(as.integer(counts[summary$dcst_label])),summary$n),
            all.equal(summary$prop, summary$n/nrow(asset$assignments)) == TRUE)
 }
}
for (arg in args) {
 path <- sub('^--[^=]+=', '', arg)
 if (startsWith(arg,'--rebuilt=')) {
  files <- list.files(path,'[.]rda$',full.names=TRUE)
  stopifnot(length(files)>0L)
  for (f in files) {
   e <- new.env(); name <- load(f,e); actual<-e[[name]]; expected<-objects[[name]]
   # The R patch version is an intentionally recorded environment difference.
   if (!is.null(actual$params$generator)) actual$params$generator$R <- expected$params$generator$R
   stopifnot(identical(actual,expected))
   cat('Rebuilt object matches:',name,'\n')
  }
 } else {
  tx <- readRDS(file.path(path,'tx.13k.rds'))
  cst <- readRDS(file.path(path,'cst.tx.13k.rds'))
  sf <- read.csv(file.path(path,'tx.13k_sample_filtering_summary.csv'))
  x <- tx[v$cst$source_row,,drop=FALSE]; rownames(x)<-rownames(v$rel)
  stopifnot(identical(x,v$rel),identical(cst$Val_subCST[v$cst$source_row],v$cst$Val_subCST),
            identical(as.integer(sf$reads_after_taxa_filter[match(rownames(tx)[v$cst$source_row],sf$sample_id)]),v$reads))
  cat('All 2,000 source rows, abundance values, sub-CST labels and read totals verified.\n')
 }
}
cat('Asset/input checksums, cross-object keys, generator records and saved summaries verified.\n')
