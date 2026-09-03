# Changelog

## linf 0.3.1

- Frozen transfer now correctly matches stable feature IDs that differ
  from display labels, including reordered columns and refined
  hierarchies.
- Transfer returns `NA` when no realized child has positive feature
  abundance; it no longer invents an unretained depth-1 state. The
  existing `carry.forward.terminal.depths` argument cannot enable an
  out-of-tree fallback.
- Dense and sparse absorption use the same original-column ordering for
  positive ties, including seeded random ties. With no positive retained
  value, the fallback remains the state with greatest reference support.
- Invalid feature-metadata lengths and non-finite normalization
  tolerances are rejected explicitly. Normalization documentation now
  distinguishes unscaled small rows, exact-zero rows, and downstream
  absorption.
- Refinement and transfer recover feature identities from fitted
  transitions, so taxonomy IDs or display labels containing `__` are not
  split apart.
- Sparse matrix preparation avoids deprecated direct
  triangular-to-general coercion. The DESCRIPTION citation now includes
  authors and year.

## linf 0.3.0

This release simplifies the public API and requires the explicit dCST
object structure introduced in `linf` 0.2.0.

- [`refine.linf.csts()`](https://pgajer.github.io/linf/reference/refine.linf.csts.md)
  now supports both automatic and explicit dominance-lineage selection
  through `lineages.to.refine`; it can be called repeatedly to add
  successive hierarchy depths.
- `refine.linf.csts.iter()` has been removed in favor of repeated calls
  to
  [`refine.linf.csts()`](https://pgajer.github.io/linf/reference/refine.linf.csts.md).
- The narrow convenience wrapper `asv.to.linf.csts()` has been removed.
  Use
  [`filter.asv()`](https://pgajer.github.io/linf/reference/filter.asv.md),
  [`normalize.linf()`](https://pgajer.github.io/linf/reference/normalize.linf.md),
  and
  [`linf.csts()`](https://pgajer.github.io/linf/reference/linf.csts.md)
  explicitly so that each stage and its dCST parameters remain visible.
- The standalone formatter `latex.linf.csts()` has been removed; dCST
  summary tables can be formatted with general reporting tools such as
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html).
- Pre-0.2 flat dCST objects are no longer accepted. Functions consuming
  fitted dCSTs now require the explicit active, pure-policy, and
  absorb-policy hierarchies produced by
  [`linf.csts()`](https://pgajer.github.io/linf/reference/linf.csts.md).

## linf 0.2.0

CRAN release: 2026-08-21

This release adopts the dominance-sample-set and dominance-lineage
terminology used in the methods paper. It is intentionally a breaking
release: obsolete cell-based names and compatibility aliases are not
retained.

- `linf.cells()` is replaced by
  [`linf.dominant.features()`](https://pgajer.github.io/linf/reference/linf.dominant.features.md).
- dCST assignments now use `lineage.id`, `lineage.label`, `lineage.ids`,
  and `lineage.labels`. Depth-1 feature assignments use
  `depth1.feature.index`.
- Retained and provisional depth-1 features now use the
  `retained.feature.*` and `provisional.feature.*` fields.
- Pure- and absorb-policy fields now end in `.pure` and `.absorb`,
  respectively.
- `cells.to.refine` is replaced by `lineages.to.refine`.
- Landmark summaries are returned in `lineages`, with `lineage.id`,
  `lineage.label`, and `lineage.size` columns.
- `collapse.rare()` and `expand.rare()` are replaced by
  `dcst.view(csts, view = "absorb")` and
  `dcst.view(csts, view = "pure")`.
- The deprecated `"rare"` low-frequency policy and view alias has been
  removed; use `"pure"`.
