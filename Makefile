######################
# 
# Created: 2016.03.26
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
######################

############### FLAGS ###############

VMAJOR 						 = 0
VMINOR 						 = 1
VPATCH  					 = 0
VDEV 							 = 
VERSION 					 = $(VMAJOR).$(VMINOR).$(VPATCH)$(VDEV)
TODAY 						:= $(shell date +%Y-%m-%d)

DOCKER 						?= $(shell which docker)

PKG_NAME 					:= fromo
PKG_LCNAME 				:= $(shell echo $(PKG_NAME) | tr 'A-Z' 'a-z')
PKG_VERSION				:= $(VERSION)
PKG_SRC 					:= $(shell basename $(PWD))

PKG_TGZ 					 = $(PKG_NAME)_$(PKG_VERSION).tar.gz
PKG_INSTALLED 		 = .$(basename $(basename $(PKG_TGZ))).installed
PKG_CRANCHECK 		 = $(basename $(basename $(PKG_TGZ))).crancheck

ALL_CPP 					 = $(wildcard src/*.cpp)
SRC_CPP 					 = $(filter-out src/RcppExports%,$(ALL_CPP))
EXPORTS_CPP				 = $(filter src/RcppExports%,$(ALL_CPP))

ALL_RD  					 = $(wildcard man/*.Rd)
PKG_DEPS 					 = $(ALL_CPP)
PKG_DEPS 					+= $(ALL_RD)

############## DEFAULT ##############

.DEFAULT_GOAL 	:= help

############## MARKERS ##############

.PHONY   : help
.PHONY   : build attributes
.SUFFIXES: 
.PRECIOUS: %.cpp .docker_img

############ BUILD RULES ############

help:  ## generate this help message
	@grep -E '^([a-zA-Z_-]+\s*)+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

build : $(PKG_TGZ) ## build the package .tar.gz file

attributes : $(EXPORTS_CPP) ## build the file src/RcppExports.cpp

$(EXPORTS_CPP) : $(SRC_CPP)
	r -l Rcpp -e 'compileAttributes(".")'

$(PKG_TGZ) : $(PKG_DEPS)
	r -l devtools -e 'build(".",path=".");'

$(ALL_RD) : 
	r -l devtools -e 'document(".");'

$(PKG_INSTALLED) : .%.installed : %.tar.gz
	r -e 'install.packages("$<");'

installed : $(PKG_INSTALLED) ## install the package

README.md : nodist/README.Rmd
	r -l Rcpp -l knitr -l devtools -e 'build();install();setwd("$(<D)");knit("$(<F)");'
	mv nodist/README.md $@

.docker_img : docker/Dockerfile
	$(DOCKER) build --rm -t $(USER)/$(PKG_LCNAME)-crancheck docker
	touch $@

%.crancheck : %.tar.gz .docker_img
	$(DOCKER) run -it --rm --volume $(PWD):/srv:ro $(USER)/$(PKG_LCNAME)-crancheck $< > $@

check: $(PKG_CRANCHECK) ## check the package as CRAN.

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
