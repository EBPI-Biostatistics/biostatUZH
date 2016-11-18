################################################################################
## Author: Sebastian Meyer [sebastian *.* meyer *a*t* ifspm *.* uzh *.* ch]
## Time-stamp: <[Makefile] 2016-11-18 12:20 (CET) by SM>
## Build/Check/Install a Package
################################################################################

R := R
PKG := biostatUZH
VERSION := $(strip $(shell grep "^Version:" ${PKG}/DESCRIPTION | cut -f 2 -d ":"))

build:
	$R CMD build ${PKG}

check: build
	_R_CHECK_FORCE_SUGGESTS_=FALSE $R CMD check ${PKG}_${VERSION}.tar.gz
	@cd ${PKG}.Rcheck; nwarn=`grep -c "^Warning" ${PKG}-Ex.Rout`; \
	if [ $$nwarn -gt 0 ]; then echo "\n\tWARNING: $$nwarn" \
        "warning(s) thrown when running examples,\n" \
	"\t         see file ${PKG}.Rcheck/${PKG}-Ex.Rout\n"; fi

install: build
	$R CMD INSTALL ${PKG}_${VERSION}.tar.gz

manual:
	$R CMD Rd2pdf --batch --force --output=${PKG}.pdf ${PKG}

## all targets are "phony"
.PHONY: build check install manual
