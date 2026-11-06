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
| Unordered dominance sets and support filtering | `R/udcsts.R` |
| Node identity and policy views | `R/lineage_nodes.R`, `R/policy_helpers.R`, `R/policy_views.R` |
| Printed summaries | `R/csts_methods.R` |
| Landmarks, frozen transfer, embedding | `R/landmarks.R`, `R/transfer_dcsts.R`, `R/hypercube_embedding.R` |
| Shared small helpers | `R/utils.R` |
| Dataset descriptions and builders | `R/data.R`, `data-raw/README.md`, `data-raw/*.R`, `data-raw/*.py` |
| Installed guides and tutorials | `vignettes/*.Rmd` |
| Website-only embedding article | `vignettes/articles/valencia-hypercube-embedding.Rmd` |
| README and website navigation | `README.Rmd`, `_pkgdown.yml` |

`man/*.Rd` and `NAMESPACE` come from roxygen comments; edit the R sources.
`README.md` comes from `README.Rmd`. Website and installed vignette HTML are
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
make audit-installed-guides
```

Developer outputs live outside the source tree. Set a checkout's output root
with `git config --local linf.qaDir /absolute/path/to/linf-qa`. Otherwise tools
use R's user cache directory for linf. `LINF_QA_DIR` overrides that setting;
Make also accepts `QA_DIR=/absolute/path/to/linf-qa`. Run `make qa-path` to
print the resolved location. Builds, checks and logs use its `build/`,
`checks/` and `logs/` subdirectories. Keep archived release evidence separate
from this working directory: `make clean` removes its current archive, check
directory and logs. It never removes archived releases outside that directory.

`make check-clean` rebuilds a source archive and runs the full `--as-cran`
check with user compiler overrides disabled. Supply `R_TIDYCMD` if HTML Tidy
is not discoverable. `make audit-installed-guides` checks navigation in the
installed package from that check.
Use `testthat::test_local(filter = "...")` for a focused development check;
keep regression files named after the behavior they exercise.

Build the local site with `make site`, validate it with `make audit-site`, and
inspect the changed pages in a browser, including standalone installed HTML and narrow
viewports. Ordinary article builds do not run the optional `grip` optimizer.
Its opt-in requirements and diagnostics are documented in the article.
The local site goes to `site/` under the QA root; `make vignette-previews`
uses `vignettes/`, and article audits use `articles/`. `make readme` updates
the tracked Markdown without leaving an HTML preview in the repository.
The GitHub website workflow retains its disposable `docs/` build directory.

The GitHub workflow checks five platforms. Pushes build the site but do not
publish it; publication is an explicit `workflow_dispatch` release action.
CRAN submission is a separate maintainer decision.

## Performance investigations

The fixed workloads and runner live in `tools/benchmarks/`. Run a baseline
from an isolated checkout before changing an algorithm:

```sh
python3 tools/benchmarks/run.py --source /path/to/baseline --output "$(make -s qa-path)/benchmarks/before"
python3 tools/benchmarks/run.py --source . --output "$(make -s qa-path)/benchmarks/after"
```

Keep workload generation unchanged between runs. Compare complete RDS results,
not just runtime. `inst/PERFORMANCE.md` describes the measured cases, memory
measurement boundary and limits of the current evidence. Benchmarks are not
ordinary unit tests or promised service-level performance.
