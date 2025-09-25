## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(linf)
set.seed(1)

## -----------------------------------------------------------------------------
S.counts <- matrix(rpois(30, lambda = 5), nrow = 10, ncol = 3,
                   dimnames = list(paste0("s",1:10), paste0("t",1:3)))
Z <- normalize.linf(S.counts)
apply(Z, 1, max)  # nonzero rows => 1; zero rows => 0

## -----------------------------------------------------------------------------
# Give columns meaningful names
colnames(Z) <- c("Lactobacillus", "Gardnerella", "Anaerobes")

cells <- linf.cells(Z)
str(cells)

# Index vs label views
head(cells$index)
head(cells$label)

# Observed levels (subset of full levels, in column order)
cells$observed.levels

## -----------------------------------------------------------------------------
T <- rbind(c(0.7, 0.3, 0.0),
           c(0.0, 0.0, 0.0),
           c(0.5, 0.5, 0.0))  # tie -> first max
colnames(T) <- c("A","B","C")
linf.cells(T)

## -----------------------------------------------------------------------------
# Build a toy matrix with an obvious dominant cell
A <- matrix(c(5, 1, 0), nrow = 6, ncol = 3, byrow = TRUE)
B <- matrix(c(1, 5, 0), nrow = 2, ncol = 3, byrow = TRUE)
C <- matrix(c(1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
S <- rbind(A, B, C)
S <- sweep(S, 1, rowSums(S), "/")
colnames(S) <- c("Dom1", "Dom2", "Dom3")

res <- linf.csts(S, n0 = 5)
names(res)

# Kept cells (by index and label)
res$kept.cells.idx
res$kept.cells.lbl

# Final (truncated) assignments
table(res$cell.label, useNA = "ifany")

## -----------------------------------------------------------------------------
S2 <- diag(1, 3) # three singletons
colnames(S2) <- c("A","B","C")
res2 <- linf.csts(S2, n0 = 2)
res2$kept.cells.lbl
res2$cell.label

## ----eval=FALSE---------------------------------------------------------------
# filt <- filter.asv(S.counts, min.lib = 10, prev.prop = 0.1, min.count = 1)
# M <- normalize.linf(filt$counts)
# cst <- linf.csts(M, n0 = 50)

