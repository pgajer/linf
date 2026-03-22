.PHONY: clean build build-verbose build-log check check-clean check-fast check-examples install document attrs winbuilder-release winbuilder-devel winbuilder-oldrelease

VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := $(shell grep "^Package:" DESCRIPTION | sed 's/Package: //')
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
LOGDIR := .claude
HOMEBREW_BIN := /opt/homebrew/bin

clean:
	rm -rf $(PKGNAME).Rcheck ..Rcheck *.Rcheck
	rm -f $(TARBALL) ../$(TARBALL)
	rm -f $(LOGDIR)/*.log

attrs:
	@mkdir -p $(LOGDIR)
	@echo "No Rcpp usage detected; skipping compileAttributes()." | tee $(LOGDIR)/$(PKGNAME)_rcppattrs.log

document: attrs
	@mkdir -p $(LOGDIR)
	@echo "Running devtools::document()..."
	@R -q -e "devtools::document()" > $(LOGDIR)/$(PKGNAME)_document.log 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_document.log)"

build: clean document
	@mkdir -p $(LOGDIR)
	@echo "Building package..."
	@R CMD build . > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: clean document
	R CMD build .

build-log: clean document
	@mkdir -p $(LOGDIR)
	R CMD build . > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Build output saved to $(LOGDIR)/$(PKGNAME)_build.log"

check: build
	PATH="$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran

check-clean: build
	env R_MAKEVARS_USER=/dev/null PATH="$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran

check-fast: build
	PATH="$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	PATH="$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran --examples

install: build
	R CMD INSTALL $(TARBALL)

winbuilder-release: build
	Rscript tools/check-win-builder.R release

winbuilder-devel: build
	Rscript tools/check-win-builder.R devel

winbuilder-oldrelease: build
	Rscript tools/check-win-builder.R oldrelease
