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
