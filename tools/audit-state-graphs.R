options(rgl.useNULL=TRUE)
pkgload::load_all('.',quiet=TRUE)
ivue.source<-Sys.getenv('IVUE_SOURCE','../ivue')
if(file.exists(file.path(ivue.source,'DESCRIPTION'))) pkgload::load_all(ivue.source,quiet=TRUE)
out<-file.path(getwd(),'build','state-graphs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
e<-new.env(parent=globalenv())
rmarkdown::render('vignettes/articles/leading-taxa-state-graphs.Rmd',output_dir=out,envir=e,quiet=TRUE)
# Independently enumerate pure memberships and both sides of every admitted edge.
stopifnot(identical(as.integer(table(e$pair$pure.id)),c(5L,3L,1L,2L,1L)))
for(i in seq_len(nrow(e$g.all$edges))) {
  edge<-e$g.all$edges[i,]
  stopifnot(sum(e$pair$pure.id==edge$pair.a & e$triple$pure.id==edge$triple)==edge$i.a,
            sum(e$pair$pure.id==edge$pair.b & e$triple$pure.id==edge$triple)==edge$i.b)
}
# Subsetting the graph must not change full population marginals.
stopifnot(identical(e$g60$pair.frequencies,e$g.all$pair.frequencies),
          identical(e$g60$triple.frequencies,e$g.all$triple.frequencies))
write.csv(e$summary.table,file.path(out,'coverage-summary.csv'),row.names=FALSE)
write.csv(e$vertices90,file.path(out,'taxon-vertices-90.csv'),row.names=FALSE)
cat('PASS: pure counts, coverage tiers, absorption, vertex unions and triple evidence.\n')
