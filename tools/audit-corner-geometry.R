# Run on demand from the package root: make audit-corner-geometry
# Optional IVUE_SOURCE loads an authoritative local checkout instead of installed ivue.
options(rgl.useNULL = TRUE)
pkgload::load_all('.', quiet = TRUE)
ivue.source <- Sys.getenv('IVUE_SOURCE')
if (nzchar(ivue.source)) pkgload::load_all(ivue.source, quiet = TRUE)
out <- normalizePath('build/corner-geometry', mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
e <- new.env(parent = globalenv())
rmarkdown::render('vignettes/articles/simplex-corner-geometry.Rmd',
                  output_dir = out, envir = e, quiet = TRUE)
# Independent direct-vector checks for the equal-secondary derivation.
for (k in c(1L, 2L, 10L, 100L)) {
  x <- matrix(0, 2, 2*k+1)
  x[,1] <- c(0.8, 0.7)
  x[1,2:(k+1)] <- 0.2/k
  x[2,(k+2):(2*k+1)] <- 0.3/k
  q <- linf::normalize.linf(x)
  x[,1] <- x[,1]-1; q[,1] <- q[,1]-1
  u <- e$unit.rows(x); v <- e$unit.rows(q)
  stopifnot(abs(sum(u[1,]*u[2,]) - k/(k+1)) < 1e-12,
            abs(sum(v[1,]*v[2,])) < 1e-12)
}
# Repeating one secondary composition at different radial distances gives zero angles.
x <- rbind(c(-0.2,0.14,0.06),c(-0.4,0.28,0.12))
stopifnot(abs(sum(e$unit.rows(x)[1,]*e$unit.rows(x)[2,])-1) < 1e-12)
# The package dataset and the common projection do not change comparison membership.
stopifnot(nrow(e$P0)==2000L, ncol(e$P0)==178L,
          nrow(e$pairs)==84107L, length(unique(e$pair.group))==15L,
          all(e$partner[e$a] != e$partner[e$b]),
          !anyDuplicated(data.frame(e$pairs)),
          nrow(e$ZA)==nrow(e$ZB), ncol(e$ZA)==3L,
          max(abs(crossprod(e$basis)-diag(3)))<1e-12,
          all(rowSums(e$ZA^2)<=1+1e-12), all(rowSums(e$ZB^2)<=1+1e-12))
write.csv(e$pair.summary, file.path(out,'group-pair-angles.csv'))
write.csv(rbind(L1=e$projection.check(e$ZA,e$A),Max=e$projection.check(e$ZB,e$B)),
          file.path(out,'projection-check.csv'))
cat('PASS: article rendered; direct angles match identities; disjoint, overlapping,\n')
cat('unequal-mass examples and all 84,107 paired empirical comparisons checked.\n')
cat('Preview:',file.path(out,'simplex-corner-geometry.html'),'\n')
