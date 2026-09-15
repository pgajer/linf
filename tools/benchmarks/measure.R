# One fresh process per operation; invoked by run.py.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
source.root <- normalizePath(args[[1]])
case <- args[[2]]; operation <- args[[3]]; output <- args[[4]]
pkgload::load_all(source.root, quiet = TRUE)
set.seed(240915)
n <- switch(case, small = 1000L, sparse = 5000L, branching = 3000L)
p <- switch(case, small = 30L, sparse = 200L, branching = 300L)
X <- matrix(rexp(n * p), n, p)
if (case != "small") X[matrix(runif(n * p) > 0.03, n, p)] <- 0
# Supported parents with heterogeneous remaining features, including empty rows.
X[cbind(seq_len(n), rep(seq_len(p), length.out = n))] <- 10
X[seq(100L, n, 100L), ] <- 0
dimnames(X) <- list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
if (case == "sparse") X <- methods::as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
fit <- linf.csts(X, n0 = 2, low.freq.policy = "absorb")
for (i in 2:3) fit <- refine.linf.csts(X, fit, n0 = 2, refinement.factor = 2,
                                      low.freq.policy = "absorb", verbose = FALSE)
gc()
start <- proc.time()[["elapsed"]]
result <- switch(operation,
  fit = {f <- linf.csts(X, n0 = 2, low.freq.policy = "absorb");
         for (i in 2:3) f <- refine.linf.csts(X, f, n0 = 2, refinement.factor = 2,
                                          low.freq.policy = "absorb", verbose = FALSE); f},
  transfer = transfer.dcsts(X, fit),
  landmarks = linf.landmarks(X, fit, depth = 3),
  embedding = linf.hypercube.embedding(X, reference = 1))
elapsed <- proc.time()[["elapsed"]] - start
saveRDS(result, paste0(output, ".rds"))
write.csv(data.frame(case = case, operation = operation, samples = n, features = p,
                     nonzero.fraction = sum(X != 0) / (n * p), depth = fit$depth,
                     groups = length(unique(fit$lineage.id)), elapsed.seconds = elapsed,
                     input.bytes = as.numeric(object.size(X)),
                     output.bytes = as.numeric(object.size(result))),
          paste0(output, ".csv"), row.names = FALSE)
