.PHONY: clean build build-verbose build-log check check-clean check-fast check-examples install document attrs winbuilder-release winbuilder-devel winbuilder-oldrelease

VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := $(shell grep "^Package:" DESCRIPTION | sed 's/Package: //')
QA_DIR ?= $(shell Rscript tools/qa-paths.R --print)
QA_DIR := $(QA_DIR)
ifeq ($(strip $(QA_DIR)),)
$(error Could not resolve an external QA directory; set LINF_QA_DIR or QA_DIR)
endif
export LINF_QA_DIR := $(QA_DIR)
BUILDDIR := $(QA_DIR)/build
CHECKDIR := $(QA_DIR)/checks
TARBALL := $(BUILDDIR)/$(PKGNAME)_$(VERSION).tar.gz
LOGDIR := $(QA_DIR)/logs

.PHONY: qa-path check-installed-guides site check-site
clean:
	rm -rf "$(CHECKDIR)/$(PKGNAME).Rcheck"
	rm -f "$(TARBALL)"
	rm -f "$(LOGDIR)"/*.log

qa-path:
	@printf '%s\n' "$(QA_DIR)"

attrs:
	@mkdir -p "$(LOGDIR)"
	@echo "No Rcpp usage detected; skipping compileAttributes()." > "$(LOGDIR)/$(PKGNAME)_rcppattrs.log"

document: attrs
	@mkdir -p "$(LOGDIR)"
	@echo "Running devtools::document()..."
	@R -q -e "devtools::document()" > "$(LOGDIR)/$(PKGNAME)_document.log" 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_document.log)"

build: repo-hygiene clean document
	@mkdir -p "$(LOGDIR)" "$(BUILDDIR)"
	@Rscript tools/check-guides.R
	@echo "Building package..."
	@cd "$(BUILDDIR)" && R CMD build "$(CURDIR)" > "$(LOGDIR)/$(PKGNAME)_build.log" 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: repo-hygiene clean document
	@mkdir -p "$(BUILDDIR)"
	cd "$(BUILDDIR)" && R CMD build "$(CURDIR)"

build-log: repo-hygiene clean document
	@mkdir -p "$(LOGDIR)" "$(BUILDDIR)"
	cd "$(BUILDDIR)" && R CMD build "$(CURDIR)" > "$(LOGDIR)/$(PKGNAME)_build.log" 2>&1
	@echo "Build output saved to $(LOGDIR)/$(PKGNAME)_build.log"

check: build
	@mkdir -p "$(CHECKDIR)"
	R CMD check "$(TARBALL)" --as-cran --output="$(CHECKDIR)"

check-clean: build
	@mkdir -p "$(CHECKDIR)"
	env R_MAKEVARS_USER=/dev/null R CMD check "$(TARBALL)" --as-cran --output="$(CHECKDIR)"

check-fast: build
	@mkdir -p "$(CHECKDIR)"
	R CMD check "$(TARBALL)" --as-cran --no-examples --no-tests --no-manual --output="$(CHECKDIR)"

check-examples: build
	@mkdir -p "$(CHECKDIR)"
	R CMD check "$(TARBALL)" --as-cran --no-tests --no-vignettes --no-manual --output="$(CHECKDIR)"

install: build
	R CMD INSTALL "$(TARBALL)"

winbuilder-release: build
	Rscript tools/check-win-builder.R release

winbuilder-devel: build
	Rscript tools/check-win-builder.R devel

winbuilder-oldrelease: build
	Rscript tools/check-win-builder.R oldrelease

.PHONY: repo-hygiene
repo-hygiene:
	python3 tools/check-no-manuscripts.py

.PHONY: readme check-guides vignette-previews
readme:
	Rscript tools/render-readme.R

check-guides:
	Rscript tools/check-guides.R

check-installed-guides:
	Rscript tools/check-guides.R --library="$(CHECKDIR)/$(PKGNAME).Rcheck"

site:
	Rscript tools/build-site.R

check-site:
	python3 tools/check-guide-html.py "$(QA_DIR)/site/articles"

vignette-previews: check-guides
	Rscript tools/render-vignette-previews.R

.PHONY: check-embedding-article
check-embedding-article:
	Rscript tools/check-embedding-article.R

.PHONY: check-data
check-data:
	Rscript tools/check-data.R

.PHONY: check-corner-geometry
check-corner-geometry:
	Rscript tools/check-corner-geometry.R

.PHONY: state-graphs-vignette
state-graphs-vignette:
	Rscript tools/check-state-graphs.R
