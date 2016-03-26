######################
# 
# Created: 2016.03.26
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
######################

############### FLAGS ###############

ALL_CPP 					= $(wildcard src/*.cpp)
SRC_CPP 					= $(filter-out src/RcppExports%,$(ALL_CPP))
EXPORTS_CPP				= $(filter src/RcppExports%,$(ALL_CPP))

############## DEFAULT ##############

default : all

############## MARKERS ##############

.PHONY   : 
.SUFFIXES: 
.PRECIOUS: %.cpp

############ BUILD RULES ############

$(EXPORTS_CPP) : $(SRC_CPP)
	r -l Rcpp -e 'compileAttributes(".")'

README.md : nodist/README.Rmd
	r -l Rcpp -l knitr -l devtools -e 'build();install();setwd("$(<D)");knit("$(<F)");'
	mv nodist/README.md $@

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
