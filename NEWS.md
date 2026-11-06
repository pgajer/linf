# linf 0.4.0

## New features

* Add `linf.udcsts()` to group samples by a chosen number of their most
  abundant features, ignoring their order. For example, samples led by A then B
  and samples led by B then A count towards the same `{A, B}` group. The minimum
  group size is applied after these samples are combined. Small groups can be
  marked as rare, excluded, or reassigned; the original memberships are also
  kept. The frequency table reports how many samples each group contains and
  what percentage of samples are covered by groups at least that large.
* Add `linf.udcst.graph()` to build a network of groups defined by their two
  leading features. By default, two pairs are linked only when both occur
  among samples with the same three leading features. Results include sample
  counts, measures of overlap, and each sample's group. Counts use all input
  samples, even when only some pairs are kept in the network. Triangles are
  reported when all three pairs are supported by the same three-feature group.
  Each link has length one. Unconnected pairs and separate parts of the network
  are kept.
* `transfer.dcsts()` assigns new samples to previously defined groups without
  changing their definitions. It now returns group IDs as well as readable group
  names, with one row per sample and one column per grouping level. These IDs
  can be used to match outputs from the same grouping result. Do not use them
  to match groups from independently created results.
* Assignment results also show which features were matched between the old and
  new data, why assignment stopped for each sample, and how often equal feature
  values required a choice between groups.
* Grouping results now record the settings used at each level and which groups
  were selected for further subdivision. Settings that were not recorded by
  older versions are marked as unavailable.

## Fixes and changes to existing behavior

* Keep feature IDs separate from displayed group names when forming groups,
  subdividing them, assigning new samples, and selecting representative samples.
  Names containing `__`, repeated names, and features named like the rare-sample
  category no longer cause groups to be confused. A suffix is added only when
  needed to distinguish otherwise identical group names. Older results are
  updated automatically when their groups can be identified unambiguously;
  otherwise, recreate the grouping from the original data. Categories that
  collect rare samples are not subdivided further.
* When a function needs the data used to form groups, check the numbers of rows
  and columns. For groupings created with this version, also check their
  original names and order. Reordering named rows or columns now gives an error.
  Unnamed rows and columns, and older results without saved names, continue to
  be matched by position.
* Correct assignment of new samples when feature IDs differ from displayed
  names, including when the new data's columns have been reordered or the
  original groups have been subdivided. Unnamed columns use the same generated
  IDs as when the groups were created. Names or IDs used to match features must
  not be missing, empty or repeated; displayed names may repeat when matching
  uses IDs.
* New samples can only be assigned to groups that already exist. If there is
  no eligible group, or none of its defining features has positive abundance,
  assignment stops and returns `NA` at that level and all later levels. This
  also applies when `carry.forward.terminal.depths = TRUE`.
* Make reassignment of samples from small groups agree between ordinary
  matrices and sparse matrices, which avoid storing most zero entries. This
  includes choices between equal values made using the same random seed. When
  all eligible feature values are zero, the sample still goes to the group
  with the most samples in the data used to define the groups.
* Make `linf.hypercube.embedding()` handle extremely large finite input values
  and very small nonzero output values more reliably. The new optional
  `log.lambda` argument and corresponding result allow scaling values to be
  reused even when they are too extreme to store directly; `lambda` is `NA`
  in those cases. Inputs with fewer than two features give an error.
* Check the data supplied to `filter.asv()` and require each cutoff to be a
  single valid value. Correctly handle samples with zero total counts when
  deciding whether a feature is present from its proportion of the sample.
  Empty results still contain the documented information and original row and
  column numbers. Nonnegative fractional counts remain accepted.
* Summaries now report the numbers of input, assigned, unassigned and rare
  samples. Statistics that cannot be calculated are returned as `NA` without
  warnings. The existing `total.samples` field continues to count assigned
  samples. Printed results limit the number of groups shown and explain results
  with no samples or only rare samples.
* Check that feature IDs and labels have one entry per column, and reject
  missing or infinite values for the rescaling cutoff. Preparing sparse
  matrices no longer uses an outdated conversion.

## Guides, example data and speed

* Add guides to the functions and example datasets, with runnable examples,
  links to package help, and tables that can be navigated with a keyboard.
  Explain how defining groups differs from assigning new samples to existing
  groups, how representative samples are selected from observed data, when
  scaling leaves rows unchanged, and why comparisons with the data used to
  define groups do not measure performance on new samples.
* Expand the tutorials and companion articles on the geometry of proportions,
  grouping by leading features regardless of order, and interpreting group
  networks. Installed guides link to available article source files; the
  website links to the formatted articles.
* Explain how the gut example samples were selected and why they cannot be
  used to draw conclusions about the broader population. Compare minimum group
  sizes with and without reassignment of rare samples, using the same inputs.
  Cite evidence about how American Gut Project samples were handled without
  claiming causes that cannot be established from this example subset.
* Record the exact version of the VALENCIA data-preparation script. Include
  data sources and file checksums for detecting changes. Add
  `valencia2k$cst$source_row` to identify the source row for each sample;
  existing sample IDs, order, abundances and read totals are unchanged.
  Saved results record the package versions that produced them. Corrections
  to choices between tied values change two assignments at level 2 and ten
  at level 3, across ten source rows. The saved summaries are updated to match.
* Reduce repeated calculations when assigning new samples and selecting
  representative samples. All 36 performance-test results match the preceding
  implementation. Report execution times and memory used by the whole R
  process, along with the examples on which those measurements were made.
* The visualization article includes the original coordinates, a static preview
  and information about which groups are connected. When selecting nearest
  neighbors, it includes all points tied at the final selected distance. Its
  optional example of adjusting point positions uses the current
  `grip::edge.kk()` function and records the algorithm's progress and stopping
  information.

# linf 0.3.0

This release simplifies the available functions and requires grouping results
created by `linf` 0.2.0 or later.

* `refine.linf.csts()` can choose groups to subdivide automatically, or you can
  choose them with `lineages.to.refine`. Call it repeatedly to add further
  levels of subgroups.
* `refine.linf.csts.iter()` has been removed. Call `refine.linf.csts()`
  repeatedly instead.
* `asv.to.linf.csts()` has been removed. Use `filter.asv()`, `normalize.linf()`,
  and `linf.csts()` in sequence to filter the data, rescale it, and form groups,
  choosing the settings for each step explicitly.
* `latex.linf.csts()` has been removed. Format group summary tables with
  general reporting tools such as `knitr::kable()`.
* Grouping results from versions before 0.2 are no longer accepted. Functions
  now require the results created by `linf.csts()`, which include the selected
  grouping and both treatments of small-group samples: kept in a rare category
  (`pure`) and reassigned to larger groups (`absorb`).

# linf 0.2.0

This release changes function and result names to match the methods paper.
It uses "dominance sample set" for samples sharing a leading feature, and
"dominance lineage" for a sequence of group subdivisions. Code using the old
names must be updated; the old names are no longer supported.

* `linf.cells()` is replaced by `linf.dominant.features()`.
* Group assignments now use `lineage.id`, `lineage.label`, `lineage.ids`, and
  `lineage.labels` for group IDs and readable names. At the first grouping
  level, `depth1.feature.index` records the column number of the leading feature.
* Information about features whose groups meet the minimum size uses names
  beginning with `retained.feature.`. Information before that size check uses
  names beginning with `provisional.feature.`.
* Result names end in `.pure` when small-group samples are kept in a rare
  category, and `.absorb` when they are reassigned to larger groups.
* `cells.to.refine` is replaced by `lineages.to.refine`, which names the groups
  to subdivide.
* Summaries of representative samples are returned in `lineages`, with
  `lineage.id`, `lineage.label`, and `lineage.size` columns identifying each
  group, its readable name, and its size.
* `collapse.rare()` and `expand.rare()` are replaced by
  `dcst.view(csts, view = "absorb")` and `dcst.view(csts, view = "pure")`,
  respectively, to select results with or without reassignment of rare samples.
* The `"rare"` option for handling small groups has been removed; use `"pure"`
  to keep rare samples separate from the retained groups.
