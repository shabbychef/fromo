/*

  This file is part of fromo.
  
  fromo is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  fromo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with fromo.  If not, see <http://www.gnu.org/licenses/>.

  running operations.
 
  Created: 2024.11.20
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_TWO_RUNCODE__
#define __DEF_TWO_RUNCODE__

#include "common.h"
#include "two_welford.h"
#include "two_running.h"

#endif /* __DEF_TWO_RUNCODE__ */

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// running correlation, regression and so on.

//' @title
//' Compute covariance, correlation, regression over a sliding window
//' @description
//' Computes 2nd moments and comoments, as well as the means, over
//' an infinite or finite sliding window, returning a matrix with the 
//' correlation, covariance, regression coefficient, and so on.
//' 
//' @param x a vector
//' @param y a vector
//' @param window the window size. if given as finite integer or double, passed through.
//' If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, equivalent
//' to an infinite window size. If negative, an error will be thrown.
//' @param restart_period the recompute period. because subtraction of elements can cause
//' loss of precision, the computation of moments is restarted periodically based on 
//' this parameter. Larger values mean fewer restarts and faster, though less accurate
//' results. 
//' @param na_rm whether to remove NA, false by default.
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent moments from being computed on too few observations.
//' Defaults to zero, meaning no restriction.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the standard errors computation. These are subtracted from the number of
//' observations. 
//'
//' @details
//'
//' Computes the correlation or covariance, or OLS regression coefficients and 
//' standard errors. 
//' These are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//'
//' @return Typically a matrix, usually only one row of the output value. More specifically:
//' \describe{
//' \item{running_covariance}{Returns a single column of the covariance of \code{x} and \code{y}.}
//' \item{running_correlation}{Returns a single column of the correlation of \code{x} and \code{y}.}
//' \item{running_covariance_3}{Returns three columns: the variance of \code{x}, the covariance of \code{x} and \code{y}, and the
//' variance of \code{y}, in that order.}
//' \item{running_regression_slope}{Returns a single column of the slope of the OLS regression.}
//' \item{running_regression_intercept}{Returns a single column of the intercept of the OLS regression.}
//' \item{running_regression_fit}{Returns two columns: the regression intercept and the regression slope of the OLS regression.}
//' \item{running_regression_diagnostics}{Returns five columns: the regression intercept, the regression slope, the regression standard error, 
//' the standard error of the intercept, the standard error of the slope of the OLS regression.}
//' }
//'
//' @examples
//' x <- rnorm(1e5)
//' y <- rnorm(1e5) + x
//' rho <- running_correlation(x, y, window=100L)
//'
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @template param-heywood
//' @template note-wts
//' @template note-heywood
//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_correlation(SEXP x, SEXP y, 
                                  SEXP window = R_NilValue, 
                                  Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                  bool na_rm=false, int min_df=0, int restart_period=100,
                                  bool check_wts=false, 
                                  bool check_negative_moments=true) {
    int wins=get_wins(window);
    const double used_df = 0.0;
    const bool normalize_wts = false;
    NumericMatrix preval = two_runQMCurryTwo<ret_correlation>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                              na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//2FIX: maybe the used_df should be 1 here? to match R's builtin cov?
//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_covariance(SEXP x, SEXP y, 
                                 SEXP window = R_NilValue, 
                                 Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                 bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                                 bool check_wts=false, bool normalize_wts=true,
                                 bool check_negative_moments=true) {
    int wins=get_wins(window);
    NumericMatrix preval = two_runQMCurryTwo<ret_covariance>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                             na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_covariance_3(SEXP x, SEXP y, 
                                   SEXP window = R_NilValue, 
                                   Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                   bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                                   bool check_wts=false, bool normalize_wts=true,
                                   bool check_negative_moments=true) {
    int wins=get_wins(window);
    NumericMatrix preval = two_runQMCurryTwo<ret_covariance_matrix>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                    na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}


//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_regression_slope(SEXP x, SEXP y, 
                                       SEXP window = R_NilValue, 
                                       Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                       bool na_rm=false, int min_df=0, int restart_period=100,
                                       bool check_wts=false, 
                                       bool check_negative_moments=true) {
    int wins=get_wins(window);
    const double used_df = 0.0;
    const bool normalize_wts = false;
    NumericMatrix preval = two_runQMCurryTwo<ret_regression_slope>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                   na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_regression_intercept(SEXP x, SEXP y, 
                                           SEXP window = R_NilValue, 
                                           Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                           bool na_rm=false, int min_df=0, int restart_period=100,
                                           bool check_wts=false, 
                                           bool check_negative_moments=true) {
    int wins=get_wins(window);
    const double used_df = 0.0;
    const bool normalize_wts = false;
    NumericMatrix preval = two_runQMCurryTwo<ret_regression_intercept>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                       na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_regression_fit(SEXP x, SEXP y, 
                                     SEXP window = R_NilValue, 
                                     Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                     bool na_rm=false, int min_df=0, int restart_period=100,
                                     bool check_wts=false,
                                     bool check_negative_moments=true) {
    int wins=get_wins(window);
    const double used_df = 0.0;
    const bool normalize_wts = false;
    NumericMatrix preval = two_runQMCurryTwo<ret_regression_fit>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                 na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//' @rdname two_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_regression_diagnostics(SEXP x, SEXP y, 
                                             SEXP window = R_NilValue, 
                                             Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                             bool na_rm=false, int min_df=0, double used_df=2.0, int restart_period=100,
                                             bool check_wts=false, bool normalize_wts=true,
                                             bool check_negative_moments=true) {
    int wins=get_wins(window);
    NumericMatrix preval = two_runQMCurryTwo<ret_regression_diagnostics>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                         na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
