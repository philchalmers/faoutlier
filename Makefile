PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: install

build:
	cd ..;\
	R CMD build $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGSRC)

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran
	make clean

test:
	Rscript -e 'library("testthat",quietly=TRUE);library("faoutlier",quietly=TRUE);require("methods",quietly=TRUE);test_dir("tests/testthat")' 

paralleltest:
	Rscript -e 'library("testthat",quietly=TRUE);library("faoutlier",quietly=TRUE);require("methods",quietly=TRUE);setCluster();test_dir("tests/testthat")' 


clean:
	$(RM) src/*.o
	$(RM) src/*.so
	$(RM) ../$(PKGNAME)_$(PKGVERS).tar.gz
	$(RM) -r ../$(PKGNAME).Rcheck/


