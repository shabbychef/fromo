######################
# 
# Created: 2016.03.26
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
######################

############### FLAGS ###############

VMAJOR 						 = 0
VMINOR 						 = 2
VPATCH  					 = 0
#VDEV 							 = .6797
VDEV 							 =
PKG_NAME 					:= fromo

RPKG_USES_RCPP 		:= 1

include ./rpkg_make/Makefile

	# r -l Rcpp -l knitr -l devtools -e 'setwd("$(<D)");if (require(knitr)) { knit("$(<F)") }'
nodist/%.csv nodist/%.md : nodist/%.Rmd $(PKG_INSTALLED) 
	$(DOCKER) run -it --rm \
		--volume $(PWD)/nodist:/srv:rw \
		--volume $$(readlink -f $(RLIB_D)):/opt/R/lib:rw \
		$(DOCKER_ENV) \
		--entrypoint="r" $(USER)/$(PKG_LCNAME)-crancheck \
		"-l" "knitr" "-l" "$(PKG_NAME)" \
		"-e" 'setwd(".");if (require(knitr)) { knit("$(<F)") }'

nodist/timings_$(PKG_VERSION).csv : nodist/timings.csv
	cp $< $@

.PHONY : timings

timings : nodist/timings_$(PKG_VERSION).csv ## save timings for performance regression checking.

.PHONY : ref_timings

ref_timings : nodist/ref_timings.md  ## timing against reference implementations of Welford sd

# experimenting with building README.md in docker. 
# not working yet b/c I do not have the requisite packages in my docker image. sigh.
reame : $(PKG_INSTALLED) $(DOCKER_IMG)
	$(DOCKER) run -it --rm \
		--volume $(PWD):/srv:rw \
		--volume $$(readlink -f $(RLIB_D)):/opt/R/lib:rw \
		$(DOCKER_ENV) \
		--entrypoint="r" $(USER)/$(PKG_LCNAME)-crancheck \
		"-l" "knitr" "-l" "$(PKG_NAME)" \
		"-e" 'setwd(".");if (require(knitr)) { knit("README.Rmd") }'

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
