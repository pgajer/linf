# linf 0.3.1 (unreleased)

- Fitted hierarchies now store explicit node paths independently of displayed
  lineage strings. Literal separator-containing feature IDs and refined paths
  remain distinct, including real features sharing the rare-category name.
  Only colliding readable IDs/labels gain a node suffix. Unambiguous older fits
  upgrade on use; ambiguous legacy fits require rebuilding from their matrix.
  Synthetic parents remain terminal within their own stored policy view.
- Count filtering validates data and scalar thresholds, handles zero-total rows
  under relative prevalence, and preserves fields and indices in empty results.
  Fractional nonnegative count-like values remain supported.
- Hierarchy summaries include input, assigned, unassigned and rare sample counts.
  Undefined size statistics are NA without warnings; total.samples retains its
  historical assigned-sample meaning. Printing bounds the group listing and
  explains empty and all-rare fits.
- Package help is visible in the installed topic index, and the guide audit can
  verify a built installation with --library=<path>. Remove the broken website
  PDF link and distinguish released from development documentation.
- The README computes and explains a first result near the opening. The function
  guide presents its small workflow earlier; the gut tutorial uses readable
  taxon labels and unclipped count annotations, retaining its pure policy.

- Landmark lookup now preserves complete feature IDs containing the lineage
  separator, including suffix collisions, refined and terminal lineages, and
  custom rare labels. It reads the target from explicit node metadata.
- Add installed task-oriented function and example-dataset guides, with
  executable workflows, catalog coverage checks and package-help navigation.
  Clarify normalization exceptions, fitting versus transfer, landmark scale,
  and the limits of comparisons using bundled source-training data.
- Frozen transfer now correctly matches stable feature IDs that differ from
  display labels, including reordered columns and refined hierarchies.
  Unnamed query matrices use the same synthetic IDs as hierarchy fitting.
- Transfer returns `NA` when no realized child has positive feature abundance;
  it no longer invents an unretained depth-1 state. The existing
  `carry.forward.terminal.depths` argument cannot enable an out-of-tree fallback.
- Dense and sparse absorption use the same original-column ordering for
  positive ties, including seeded random ties. With no positive retained value,
  the fallback remains the state with greatest reference support.
- Invalid feature-metadata lengths and non-finite normalization tolerances are
  rejected explicitly. Normalization documentation now distinguishes unscaled
  small rows, exact-zero rows, and downstream absorption.
- Refinement and transfer use explicit fitted feature paths,
  so taxonomy IDs or display labels containing `__` are not split apart.
- Sparse matrix preparation avoids deprecated direct triangular-to-general
  coercion. The DESCRIPTION citation now includes authors and year.

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
