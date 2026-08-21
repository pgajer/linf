# linf 0.3.0

This release simplifies the public API and requires the explicit dCST object
structure introduced in `linf` 0.2.0.

- `refine.linf.csts()` now supports both automatic and explicit
  dominance-lineage selection through `lineages.to.refine`; it can be called
  repeatedly to add successive hierarchy depths.
- `refine.linf.csts.iter()` has been removed in favor of repeated calls to
  `refine.linf.csts()`.
- The narrow convenience wrapper `asv.to.linf.csts()` has been removed. Use
  `filter.asv()`, `normalize.linf()`, and `linf.csts()` explicitly so that each
  stage and its dCST parameters remain visible.
- The standalone formatter `latex.linf.csts()` has been removed; dCST summary
  tables can be formatted with general reporting tools such as `knitr::kable()`.
- Pre-0.2 flat dCST objects are no longer accepted. Functions consuming fitted
  dCSTs now require the explicit active, pure-policy, and absorb-policy
  hierarchies produced by `linf.csts()`.

# linf 0.2.0

This release adopts the dominance-sample-set and dominance-lineage terminology
used in the methods paper. It is intentionally a breaking release: obsolete
cell-based names and compatibility aliases are not retained.

- `linf.cells()` is replaced by `linf.dominant.features()`.
- dCST assignments now use `lineage.id`, `lineage.label`, `lineage.ids`, and
  `lineage.labels`. Depth-1 feature assignments use `depth1.feature.index`.
- Retained and provisional depth-1 features now use the
  `retained.feature.*` and `provisional.feature.*` fields.
- Pure- and absorb-policy fields now end in `.pure` and `.absorb`, respectively.
- `cells.to.refine` is replaced by `lineages.to.refine`.
- Landmark summaries are returned in `lineages`, with `lineage.id`,
  `lineage.label`, and `lineage.size` columns.
- `collapse.rare()` and `expand.rare()` are replaced by
  `dcst.view(csts, view = "absorb")` and `dcst.view(csts, view = "pure")`.
- The deprecated `"rare"` low-frequency policy and view alias has been removed;
  use `"pure"`.
