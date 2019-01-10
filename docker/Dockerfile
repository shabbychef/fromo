#
# dockerfile to CRAN-check with r-dev
#
# docker build --rm -t shabbychef/fromo-crancheck .
#
# docker run -it --rm --volume $(pwd):/srv:rw fromo-crancheck
#
# Created: 2016.01.10
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

#####################################################
# preamble# FOLDUP
FROM rocker/r-base:latest
MAINTAINER Steven E. Pav, shabbychef@gmail.com
# UNFOLD

ENV DOCKERFILE_REFRESHED_AT 2018-12-28
# see http://crosbymichael.com/dockerfile-best-practices.html
#RUN echo "deb http://archive.ubuntu.com/ubuntu precise main universe" > /etc/apt/sources.list

RUN (apt-get clean -y ; \
 apt-get update -y -qq; \
 apt-get update --fix-missing );

#RUN (apt-get dist-upgrade -y ; \
#apt-get update -qq ; \
RUN (DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true apt-get install -q -y --no-install-recommends \ 
  libgs9 texlive-base texlive-binaries libcupsimage2 libcups2 file curl wget qpdf pandoc ghostscript \
  texlive-latex-extra texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texinfo \
  libmariadbclient-dev libpq-dev liblapack-dev libblas-dev libssl-dev libcurl4-openssl-dev libgit2-dev libxml2-dev libcairo2-dev; \
  sync ; \
  mkdir -p /usr/local/lib/R/site-library ; \
  chmod -R 777 /usr/local/lib/R/site-library ; \
  sync )

# rinstall somethings...
RUN (install2.r --error --skipinstalled --deps TRUE \
  devtools roxygen2 testthat curl \
  knitr formatR codetools \
  Rcpp \
  readr dplyr tidyr ggplot2 \
  moments microbenchmark PDQutils roll RcppRoll multcomp; \
  r -l 'devtools' -e 'options(unzip="internal");install_github("andrewuhl/rollingwindow");' )

# get latest version of roll to test against. as of 2019-01-09, this is 1.1.2
# from github
RUN ( r -l 'devtools' -e 'options(unzip="internal");install_github("jjf234/roll");' )

# already in there, I think. weird.
# RUN groupadd -g 1000 spav && useradd -g spav -u 1000 spav;
# USER spav
# docker has id 1000 in rocker/r-base; so use that. ugh.
USER docker

WORKDIR /srv

#####################################################
# entry and cmd# FOLDUP
# these are the default, but remind you that you might want to use /usr/bin/R instead?
# always use array syntax:
ENTRYPOINT ["R","CMD","check","--as-cran","--output=/tmp"]

# ENTRYPOINT and CMD are better together:
CMD ["/srv/*.tar.gz"]
# UNFOLD

#for vim modeline: (do not edit)
# vim:nu:et:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=Dockerfile:ft=Dockerfile:fo=croql
