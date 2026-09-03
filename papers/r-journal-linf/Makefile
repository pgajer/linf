SHELL := /bin/sh

PACKAGE_TARBALL := package/linf_0.3.0.tar.gz
ARTICLE := linf.Rmd
LIBRARY ?= $(abspath library)
R_ENV := env R_MAKEVARS_USER=/dev/null R_LIBS="$(LIBRARY)" R_LIBS_USER="$(LIBRARY)"
TIDY_DIR := $(if $(wildcard /opt/homebrew/opt/tidy-html5/bin/tidy),/opt/homebrew/opt/tidy-html5/bin,)
CHECK_PROFILE := $(abspath scripts/check-profile.R)
PDF_OUT := output/pdf/linf-r-journal.pdf
HTML_OUT := output/html/linf-r-journal.html

.PHONY: dependencies package-source package-install render citation-check rj-check audit version-check package-check submission-audit draft-bundle clean

dependencies:
	mkdir -p "$(LIBRARY)"
	$(R_ENV) Rscript scripts/install-dependencies.R

package-source:
	python3 scripts/fetch-package.py

package-install: package-source
	mkdir -p "$(LIBRARY)"
	env R_MAKEVARS_USER=/dev/null R CMD INSTALL --library="$(LIBRARY)" $(PACKAGE_TARBALL)

render: package-install
	mkdir -p build output/pdf output/html
	$(R_ENV) Rscript scripts/render-paper.R
	cp linf.pdf $(PDF_OUT)
	cp linf.html $(HTML_OUT)
	python3 scripts/artifact-manifest.py record

citation-check: render
	python3 scripts/check-citation-verification.py \
		--source $(ARTICLE) \
		--bib RJreferences.bib \
		--html citation_verification.html

rj-check: render
	$(R_ENV) Rscript scripts/run-rjtools-checks.R

audit: citation-check rj-check
	python3 scripts/readiness-scan.py

version-check: package-install
	$(R_ENV) Rscript scripts/check-public-version.R

package-check: package-source
	mkdir -p build/package-check
	env R_MAKEVARS_USER=/dev/null R_PROFILE_USER="$(CHECK_PROFILE)" \
		LINF_TIDY_DIR="$(TIDY_DIR)" _R_CHECK_CRAN_INCOMING_=false \
		R CMD check --as-cran --output=build/package-check $(PACKAGE_TARBALL)
	Rscript -e 'stopifnot("Status: OK" %in% readLines("build/package-check/linf.Rcheck/00check.log"))'

submission-audit: audit version-check package-check

draft-bundle: audit
	python3 scripts/create-draft-bundle.py

clean:
	$(R_ENV) Rscript scripts/clean-generated.R
