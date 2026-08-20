SHELL := /bin/sh

PACKAGE_ROOT := ../..
ARTICLE := linf.Rmd
LIBRARY := $(abspath library)
R_ENV := env R_MAKEVARS_USER=/dev/null R_LIBS=$(LIBRARY)
TIDY_DIR := $(if $(wildcard /opt/homebrew/opt/tidy-html5/bin/tidy),/opt/homebrew/opt/tidy-html5/bin,)
CHECK_PROFILE := $(abspath scripts/check-profile.R)
PDF_OUT := output/pdf/linf-r-journal.pdf
HTML_OUT := output/html/linf-r-journal.html

.PHONY: package-install render citation-check rj-check audit version-check package-check submission-audit draft-bundle clean

package-install:
	mkdir -p $(LIBRARY)
	env R_MAKEVARS_USER=/dev/null R CMD INSTALL --library=$(LIBRARY) $(PACKAGE_ROOT)

render: package-install
	mkdir -p build output/pdf output/html
	$(R_ENV) Rscript scripts/render-paper.R
	cp linf.pdf $(PDF_OUT)
	cp linf.html $(HTML_OUT)

citation-check: render
	python3 scripts/check-citation-verification.py \
		--source $(ARTICLE) \
		--bib RJreferences.bib \
		--html citation_verification.html

rj-check: render
	$(R_ENV) Rscript scripts/run-rjtools-checks.R

audit: citation-check rj-check
	python3 scripts/readiness-scan.py

version-check:
	$(R_ENV) Rscript scripts/check-public-version.R

package-check:
	cd $(PACKAGE_ROOT) && env R_MAKEVARS_USER=/dev/null \
		R_PROFILE_USER=$(CHECK_PROFILE) LINF_TIDY_DIR=$(TIDY_DIR) \
		$(MAKE) check-clean

submission-audit: audit version-check package-check

draft-bundle: audit
	$(R_ENV) python3 scripts/create-draft-bundle.py

clean:
	$(R_ENV) Rscript scripts/clean-generated.R
