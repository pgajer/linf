# Working on linf

Use the task-oriented guide in `vignettes/function-guide.Rmd` to understand the
public workflow. Keep stable feature identity separate from display labels and
preserve equivalence of dense and sparse paths, including tie behavior.

## Source map

| Responsibility | Maintained sources |
|---|---|
| Matrix preparation, dense/sparse validation | `R/backend_helpers.R` |
| Filtering and normalization | `R/filter_asv.R`, `R/normalize.R` |
| Dominant feature assignment and display metadata | `R/dominant_features.R`, `R/label_format.R` |
| Depth-1 fits, refinement and validation | `R/fit.R`, `R/refine.R`, `R/validate_fit.R`, `R/fit_contracts.R` |
| Node identity and policy views | `R/lineage_nodes.R`, `R/policy_helpers.R`, `R/policy_views.R` |
| Printed summaries | `R/csts_methods.R` |
| Landmarks, frozen transfer, embedding | `R/landmarks.R`, `R/transfer_dcsts.R`, `R/hypercube_embedding.R` |
| Shared small helpers | `R/utils.R` |
| Dataset descriptions and builders | `R/data.R`, `data-raw/README.md`, `data-raw/*.R`, `data-raw/*.py` |
| Installed guides and tutorials | `vignettes/*.Rmd` |
| Website-only embedding article | `vignettes/articles/valencia-hypercube-embedding.Rmd` |
| README and website navigation | `README.Rmd`, `_pkgdown.yml` |

`man/*.Rd` and `NAMESPACE` come from roxygen comments; edit the R sources.
`README.md` comes from `README.Rmd`. `docs/` and installed vignette HTML are
build outputs. Bundled `data/*.rda` are versioned outputs: regenerate in a
separate directory and compare with `tools/audit-data.R` before replacing them.
Never add unpublished study/manuscript files to this public repository.

## Build and verify

```sh
make document
make readme
make audit-guides
make audit-data
make audit-embedding-article
make check-clean
```

`make check-clean` rebuilds a source archive and runs the full `--as-cran`
check with user compiler overrides disabled. Supply `R_TIDYCMD` if HTML Tidy
is not discoverable. To check installed navigation, run
`Rscript tools/audit-guides.R --library=linf.Rcheck` after the check.
Use `testthat::test_local(filter = "...")` for a focused development check;
keep regression files named after the behavior they exercise.

Build the site with `pkgdown::build_site(preview = FALSE)` and inspect the
changed pages in a browser, including standalone installed HTML and narrow
viewports. Ordinary article builds do not run the optional `grip` optimizer.
Its opt-in requirements and diagnostics are documented in the article.

The GitHub workflow checks five platforms. Pushes build the site but do not
publish it; publication is an explicit `workflow_dispatch` release action.
CRAN submission is a separate maintainer decision.

## Performance investigations

The fixed workloads and runner live in `tools/benchmarks/`. Run a baseline
from an isolated checkout before changing an algorithm:

```sh
python3 tools/benchmarks/run.py --source /path/to/baseline --output build/before
python3 tools/benchmarks/run.py --source . --output build/after
```

Keep workload generation unchanged between runs. Compare complete RDS results,
not just runtime. `inst/PERFORMANCE.md` describes the measured cases, memory
measurement boundary and limits of the current evidence. Benchmarks are not
ordinary unit tests or promised service-level performance.
