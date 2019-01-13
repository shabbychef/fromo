dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Type: Package
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Fast Robust Moments
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Fast, numerically robust computation of weighted moments via 'Rcpp'. 
   Supports computation on vectors and matrices, and Monoidal append of moments. 
   Moments and cumulants over running fixed length windows can be computed, 
   as well as over time-based windows.
   Moment computations are via a generalization of Welford's method, as described
   by Bennett et. (2009) <doi:10.1109/CLUSTR.2009.5289161>.
Imports:
    Rcpp (>= 0.12.3),
    methods
LinkingTo: Rcpp
Suggests:
    knitr,
    testthat,
    moments,
    PDQutils,
    e1071,
    microbenchmark
RoxygenNote: 5.0.1
URL: https://github.com/shabbychef/PKG_NAME()
VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:et:syn=m4:ft=m4
