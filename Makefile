include Makefile.defs

PACKAGE = biostatUZH
VERSION = 2.2.2
TAR = $(PACKAGE)_$(VERSION).tar.gz



.PHONY: test-package update-src lib check-cran check covr manual \
	winbuild winbuild-devel webpage clean


update-src:
	sed -i -r -- 's/^Version:.*/Version: '$(VERSION)'/g' DESCRIPTION ;      
	sed -i -r -- 's/^Date:.*/Date: '`date +'%F'`'/g' DESCRIPTION ;
	$(RSCRIPT) -e "roxygen2::roxygenize(\".\")"

lib: update-src
	mkdir -p lib
	$(R) CMD INSTALL -l lib . --no-lock

test-package:
	$(RSCRIPT) -e "devtools::test('.')"

$(TAR): update-src
	$(RSCRIPT) -e "devtools::build(path = '.', args = '--compact-vignettes=both')"

check-cran: $(TAR)
	$(RSCRIPT) -e "devtools::check_built(path = './$(TAR)', cran = TRUE)"

check: $(TAR)
	$(RSCRIPT) -e "devtools::check_built(path = './$(TAR)', cran = FALSE)"

covr: 
	$(RSCRIPT) -e "covr::package_coverage()"

manual: update-src
	$(R) -e 'devtools::build_manual(pkg = ".", path = ".")'

winbuild: update-src
	$(RSCRIPT) -e "devtools::check_win_release()"

winbuild-devel: update-src
	$(RSCRIPT) -e "devtools::check_win_devel()"

webpage: update-src
	$(RSCRIPT) -e 'pkgdown::build_site()'

clean:
	rm -rf lib $(PACKAGE).Rcheck *.tar.gz

