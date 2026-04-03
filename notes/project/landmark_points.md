# Landmark Points for DCSTs in `linf`

## Overview

This note describes the idea of **landmark points** for dominant-cell / DCST
structures and explains how to compute them with the current `linf` package.

The key idea is simple: once a cell or a refined DCST leaf has been defined by
dominance of a feature, that same feature induces natural representative points
inside the cell. Endpoints are the clearest example, but the same logic also
gives minimum points and representative points near the cell mean or median.

In the current `linf` implementation:

- depth-1 dominant cells are computed by `linf.csts()`
- deeper DCST levels are computed by `refine.linf.csts()` or
  `refine.linf.csts.iter()`
- landmark points are computed by `linf.landmarks()`
- depth-1 landmarks can optionally be attached directly by
  `linf.csts(..., return.landmarks = TRUE)`

Throughout, let `X` be a nonnegative matrix with rows representing samples and
columns representing features.

## Concept

### Depth-1 cells

Let `C_i` denote the depth-1 DCST cell corresponding to feature `i`. In other
words, `C_i` is the set of rows whose dominant feature is `i`.

The natural landmark points of `C_i` are defined with respect to the same
feature `i`.

The **endpoint** of `C_i` is the sample in `C_i` where feature `i` is maximal:

```text
e_i = argmax_{x in C_i} X[x, i]
```

The corresponding **minimum endpoint** is the sample in `C_i` where feature `i`
is minimal:

```text
c_i = argmin_{x in C_i} X[x, i]
```

In the package API these are called:

- `endpoint.max`
- `endpoint.min`

I prefer these names to `co-endpoint`, because they are explicit and align well
with code.

### Depth-2 and deeper DCST cells

If a depth-2 cell is written as `C_{i,j}`, then it is obtained by first placing
rows in `C_i`, then refining that parent cell after excluding feature `i`, and
then assigning dominance of feature `j` among the remaining features.

The natural feature attached to `C_{i,j}` is therefore the **leaf feature**
`j`, not the whole path `(i, j)`. Its endpoint is

```text
e_{i,j} = argmax_{x in C_{i,j}} X[x, j]
```

and similarly the minimum endpoint is

```text
c_{i,j} = argmin_{x in C_{i,j}} X[x, j]
```

This extends recursively to arbitrary depth. If a leaf DCST cell has path

```text
C_{i_1, i_2, ..., i_k}
```

then its current natural landmark feature is the leaf feature `i_k`.

### Representative points

Besides extrema, a cell also has useful representative points near its central
tendency.

For a cell `C` with landmark feature `j`, define:

- the **mean representative** as the sample in `C` whose value `X[x, j]` is
  closest to `mean(X[C, j])`
- the **median representative** as the sample in `C` whose value `X[x, j]` is
  closest to `median(X[C, j])`

In the package API these are called:

- `mean.rep`
- `median.rep`

This first implementation is intentionally feature-specific. It uses the leaf
feature of the DCST path rather than the whole row vector `X[x, ]`. That keeps
the meaning tightly aligned with the way the cell itself is defined.

## Current Scope in `linf`

The current `linf.landmarks()` implementation supports:

- `endpoint.max`
- `endpoint.min`
- `mean.rep`
- `median.rep`

These landmarks are computed for:

- any requested DCST depth
- any requested CST/DCST view: `"active"`, `"rare"`, or `"absorb"`
- any `linf.csts` object that has been optionally refined to greater depth

The current implementation does **not** yet compute:

- local minima with basin-size constraints
- graph-based shortest paths between landmarks
- full-vector mean/median representatives based on `X[x, ]`

Those are natural future extensions, but they require additional structure that
is not yet part of the package.

## Why View Matters

`linf` supports two low-frequency views:

- `"rare"`: low-frequency cells remain explicit rare buckets
- `"absorb"`: low-frequency cells are reassigned into kept cells

Because of that, landmark computation must always be tied to a specific view.
The same row can belong to different effective leaves under `"rare"` and
`"absorb"`.

For this reason `linf.landmarks()` takes a `view` argument, and the CST object
stores parallel label and ID paths for the active, rare, and absorb views.

## How the Package Computes Landmark Points

### Step 1: Normalize the matrix if needed

The standard workflow begins with L-infinity normalization:

```r
Z <- normalize.linf(X)
```

This scales each row by its row maximum. Nonzero rows then have maximum `1`,
while zero rows remain zero.

### Step 2: Compute depth-1 DCSTs

Depth-1 cells are computed by `linf.csts()`:

```r
d1 <- linf.csts(
  Z,
  n0 = 50,
  low.freq.policy = "rare"
)
```

This returns a `linf.csts` object containing:

- active cell assignments
- the rare and absorb variants
- kept-cell summaries
- level-1 cell IDs and labels

### Step 3: Refine to depth 2 or beyond

Depth-2 DCSTs are obtained by refining the depth-1 object:

```r
d2 <- refine.linf.csts(
  Z,
  d1,
  n0 = 25,
  refinement.factor = 2,
  low.freq.policy = "rare",
  verbose = FALSE
)
```

Further depth can be added with either:

- repeated calls to `refine.linf.csts.iter()`
- repeated calls to `refine.linf.csts()` if that matches the workflow better

Internally, refinement drops the parent dominant feature(s) and recomputes
dominance on the remaining columns. This is exactly why a refined cell has a
well-defined leaf feature and therefore a natural landmark feature.

### Step 4: Compute landmark points

Use `linf.landmarks()` on the matrix and the DCST object:

```r
lm2 <- linf.landmarks(
  Z,
  d2,
  depth = 2,
  view = "rare",
  landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep")
)
```

The function identifies the leaf feature for each DCST cell from the cell ID
path, then computes the requested point(s) using the rows currently assigned to
that cell.

### Step 5: Use the small pipeline wrapper when you want all pieces together

For the common "normalize -> depth-1 DCST -> depth-2 refinement -> landmarks"
workflow, the package now provides a convenience wrapper:

```r
pipe <- linf.dcst.landmark.pipeline(
  X,
  feature.ids = paste0("asv_", seq_len(ncol(X))),
  feature.labels = format.linf.feature.labels(
    feature.ids = paste0("asv_", seq_len(ncol(X))),
    taxonomy = colnames(X)
  ),
  n0.depth1 = 50,
  n0.depth2 = 25,
  refinement.factor = 2,
  low.freq.policy = "rare",
  landmark.view = "absorb"
)
```

This returns a named list containing:

- `linf.rel`
- `dcst.depth1`
- `dcst.depth2`
- `landmarks.depth1`
- `landmarks.depth2`

The wrapper does not add any new DCST logic; it simply orchestrates the
existing `normalize.linf()`, `linf.csts()`, `refine.linf.csts()`, and
`linf.landmarks()` calls in a readable package-style interface.

## Meaning of the Returned Object

`linf.landmarks()` returns an object of class `linf.landmarks` with the
following main components.

### Metadata

- `depth`: the DCST depth used for computation
- `view`: the resolved view used for cell membership
- `sep`: the separator used in hierarchical IDs
- `rare.label`: the label used for explicit rare cells
- `feature.ids`
- `feature.labels`

### `cells`

`cells` is one row per DCST cell at the requested depth and view. It includes:

- `cell.id`
- `cell.label`
- `cell.size`
- `target.feature.id`
- `target.feature.label`
- `is.rare`
- `landmarks.computable`

This table is useful because not every cell necessarily has computable
landmarks. For example, explicit rare buckets are represented in `cells`, but
their landmarks are skipped because a rare bucket does not correspond to one
unique target feature.

### `landmarks`

`landmarks` is one row per computed landmark point. It includes:

- `cell.id`
- `cell.label`
- `landmark.type`
- `point.index`
- `point.name`
- `target.feature.id`
- `target.feature.label`
- `observed.value`
- `target.value`
- `abs.deviation`

Interpretation:

- for `endpoint.max` and `endpoint.min`, `target.value` equals
  `observed.value`
- for `mean.rep` and `median.rep`, `target.value` is the cell mean or median of
  the target feature and `abs.deviation` measures how close the chosen sample is
  to that target

## Example Workflow

```r
library(linf)

X <- matrix(
  c(
    10, 8, 2,
     9, 7, 3,
     8, 2, 7,
     7, 1, 8,
     2, 9, 1,
     3, 8, 1
  ),
  nrow = 6,
  byrow = TRUE,
  dimnames = list(
    paste0("s", 1:6),
    c("A", "B", "C")
  )
)

Z <- normalize.linf(X)

d1 <- linf.csts(
  Z,
  n0 = 2,
  low.freq.policy = "rare"
)

d2 <- refine.linf.csts(
  Z,
  d1,
  n0 = 2,
  refinement.factor = 2,
  low.freq.policy = "rare",
  verbose = FALSE
)

lm1 <- linf.landmarks(
  Z,
  d1,
  depth = 1,
  view = "rare",
  landmark.types = c("endpoint.max", "mean.rep")
)

lm2 <- linf.landmarks(
  Z,
  d2,
  depth = 2,
  view = "rare",
  landmark.types = c("endpoint.max", "endpoint.min", "median.rep")
)

lm1$cells
lm1$landmarks

lm2$cells
lm2$landmarks
```

## Convenience: Attaching Landmarks at Depth 1

For depth-1 CSTs, landmarks can be computed automatically during the initial
call:

```r
d1 <- linf.csts(
  Z,
  n0 = 50,
  low.freq.policy = "rare",
  return.landmarks = TRUE,
  landmark.types = c("endpoint.max", "endpoint.min"),
  landmark.view = "rare"
)

d1$landmarks
```

This is convenient when the main use case is depth-1 analysis. For depth-2 and
deeper structures, the recommended approach is still to call `linf.landmarks()`
explicitly after refinement.

## Practical Interpretation

These landmark points are useful because they give concrete observed samples
that summarize different parts of the geometry of a DCST cell:

- `endpoint.max` identifies the sample that most strongly expresses the leaf
  feature
- `endpoint.min` identifies the weakest expression of the leaf feature within
  the cell
- `mean.rep` gives a sample near the average level of the leaf feature
- `median.rep` gives a robust central representative

Together, these points provide a compact way to inspect how a DCST cell varies
internally without collapsing the cell to a single synthetic average vector.

## Recommended Terminology

For prose:

- endpoint
- minimum endpoint
- mean representative
- median representative

For code:

- `endpoint.max`
- `endpoint.min`
- `mean.rep`
- `median.rep`

I would avoid using `co-endpoint` as the primary API name. It is suggestive, but
less precise and less self-explanatory than `endpoint.min`.

## Future Extensions

The current implementation provides a good first layer, but several natural
extensions remain:

1. Basin-aware local minima, for example minima whose attraction basin contains
   at least `M` points.
2. Full-vector representatives based on distance to `mean(X[C, ])` or
   `median(X[C, ])`.
3. Paths between landmarks inside a nearest-neighbor or graph structure induced
   on each DCST cell.
4. Automatic attachment of deeper landmark summaries during refinement.

These are all compatible with the current design. The main architectural choice
that makes future work easier is that landmarks are implemented as a separate,
well-defined layer on top of the DCST object rather than being fused into the
core assignment logic.
