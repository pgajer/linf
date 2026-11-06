# linf 0.4.0

## New functionality

* Add `linf.udcsts()` for unordered top-feature sets. Opposite dominance orders
  pool before support filtering. Pure, symmetric set-based absorption, and
  filtering policies preserve original memberships; dense and sparse inputs
  agree. Frequency tables report cumulative coverage including frequency ties.
* Add `linf.udcst.graph()` for pair-state graphs supported by observed leading
  triples, with full-reference incidence counts, similarities, witnessed faces
  and sample membership. Edges have unit length; disconnected components and
  isolated states are preserved.
* Frozen transfer now returns lineage-ID matrices, a feature-match table and
  per-query stopping and tie diagnostics alongside its existing label matrices.
* Fitted objects record per-depth settings and refinement selections. Earlier
  settings that were not saved are marked unavailable when older fits are used.

## Correctness and compatibility

* Preserve stable feature identity independently of displayed lineage paths
  during fitting, refinement, transfer and landmark lookup. Feature IDs and
  labels containing `__`, suffix collisions, and real features sharing a rare
  category's label remain distinct. Only colliding readable lineage names gain
  a node suffix. Unambiguous older fits upgrade on use; ambiguous legacy fits
  must be rebuilt from their input matrix. Synthetic rare parents remain
  terminal in their stored policy view.
* Validate matrix dimensions and, for new fits, original input dimnames
  independently of custom feature IDs. Named axis reorderings fail explicitly;
  unnamed axes and legacy fits retain positional matching.
* Correct frozen transfer when stable feature IDs differ from display labels,
  including reordered query columns and refined hierarchies. Unnamed queries
  use the same synthetic IDs as fitting. Reject duplicate, missing or empty
  matching keys while allowing repeated display labels when matching by ID.
* Return `NA` when no realized child has positive query abundance. Transfer
  no longer invents an unretained depth-1 state; `carry.forward.terminal.depths`
  does not permit assignments outside the fitted hierarchy.
* Make dense and sparse absorption agree on positive ties, including seeded
  random ties. With no positive retained value, the fallback remains the state
  with greatest reference support.
* Stabilize hypercube calculations for extreme finite values and small positive
  radii. The optional trailing `log.lambda` argument and returned attribute
  allow reuse of scales outside the ordinary numeric range; `lambda` is `NA`
  for those scales. Inputs with fewer than two features fail explicitly.
* Validate count-filtering data and scalar thresholds, handle zero-total rows
  under relative prevalence, and preserve fields and indices in empty results.
  Fractional nonnegative count-like values remain supported.
* Report input, assigned, unassigned and rare sample counts in summaries.
  Undefined size statistics are `NA` without warnings; `total.samples` retains
  its historical assigned-sample meaning. Printing bounds the group listing
  and explains empty and all-rare fits.
* Reject invalid feature-metadata lengths and non-finite normalization
  tolerances. Sparse preparation avoids deprecated triangular-matrix coercion.

## Documentation, data and efficiency

* Add installed function and example-dataset guides with executable workflows,
  package-help navigation and keyboard-accessible tables. Clarify fitting
  versus frozen transfer, observed landmarks versus averaged profiles,
  normalization exceptions, and the limits of source-training-data comparisons.
* Expand the tutorials and companion articles on simplex-corner geometry,
  unordered states and graph interpretation. Installed guides link to available
  article sources; website builds link to rendered articles.
* Document the gut example's stratified sampling and limits on population
  inference. Compare support thresholds under both policies on the same input
  samples and cite AGP sample-handling evidence without attributing causes that
  cannot be established from the bundled subset.
* Pin the VALENCIA preparation recipe and install data provenance and checksums.
  Add `valencia2k$cst$source_row`; existing IDs, ordering, abundances and read
  totals are unchanged. Saved merged assignments record generator versions and
  correct tied choices affecting two depth-2 and ten depth-3 assignments across
  ten source rows. Their summaries are regenerated from those assignments.
* Reduce repeated work in transfer and landmark construction. All 36 complete
  benchmark results match the preceding implementation; measured timings and
  whole-process memory use are documented with their workload limits.
* The embedding article includes native coordinates, a static preview and
  connectivity information. Its optional nearest-neighbor graph includes cutoff
  ties, and its opt-in optimizer example uses the current `grip::edge.kk()` API
  while retaining optimizer diagnostics.

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
