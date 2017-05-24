######################
# 
# Created: 2016.03.26
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
######################

############### FLAGS ###############

VMAJOR 						 = 0
VMINOR 						 = 1
VPATCH  					 = 3
VDEV 							 = .3001
PKG_NAME 					:= fromo

RPKG_USES_RCPP 		:= 1

include ./rpkg_make/Makefile

nodist/%.csv nodist/%.md : nodist/%.Rmd $(PKG_INSTALLED) | tools/figure
	r -l Rcpp -l knitr -l devtools -e 'setwd("$(<D)");if (require(knitr)) { knit("$(<F)") }'

nodist/timings_$(PKG_VERSION).csv : nodist/timings.csv
	cp $< $@

.PHONY : timings

timings : nodist/timings_$(PKG_VERSION).csv ## save timings for performance regression checking.

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
