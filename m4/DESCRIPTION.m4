dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Type: Package
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", role=c("aut","cre"), 
    email="shabbychef@gmail.com"))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Fast Robust Moments
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Fast computation of moments via Rcpp. Supports computation on
   vectors and matrices, and Monoidal append of moments (NYI).
License: LGPL-3
Imports:
    Rcpp (>= 0.12.3)
LinkingTo: Rcpp
Suggests:
    testthat,
    moments,
    microbenchmark
RoxygenNote: 5.0.1
URL: https://github.com/shabbychef/PKG_NAME()
dnl VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4
