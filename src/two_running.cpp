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
//' of the centered moments computation. These are subtracted from the number of
//' observations. Not sure this is meaningful here yet.
//'
//' @details
//'
//' 2FIX
//' Computes the number of elements, the mean, and the 2nd through kth
//' centered (and typically standardized) moments, for \eqn{k=2,3,4}{k=2,3,4}. These
//' are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//'
//' Given the length \eqn{n} vector \eqn{x}, we output matrix \eqn{M} where
//' \eqn{M_{i,j}}{M_i,j} is the \eqn{order - j + 1} moment (\emph{i.e.}
//' excess kurtosis, skewness, standard deviation, mean or number of elements)
//' of \eqn{x_{i-window+1},x_{i-window+2},...,x_{i}}{x_(i-window+1),x_(i-window+2),...,x_i}.
//' Barring \code{NA} or \code{NaN}, this is over a window of size \code{window}.
//' During the 'burn-in' phase, we take fewer elements.
//'
//' @return 2FIX Typically a matrix, where the first columns are the kth, k-1th through 2nd standardized, 
//' centered moments, then a column of the mean, then a column of the number of (non-nan) elements in the input,
//' with the following exceptions:
//' \describe{
//' \item{running_cent_moments}{Computes arbitrary order centered moments. When \code{max_order_only} is set,
//' only a column of the maximum order centered moment is returned.}
//' \item{running_std_moments}{Computes arbitrary order standardized moments, then the standard deviation, the mean,
//' and the count. There is not yet an option for \code{max_order_only}, but probably should be.}
//' \item{running_cumulants}{Computes arbitrary order cumulants, and returns the kth, k-1th, through the second 
//' (which is the variance) cumulant, then the mean, and the count.}
//' }
//'
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
                                  bool check_wts=false, bool normalize_wts=true,
                                  bool check_negative_moments=true) {
    int wins=get_wins(window);
    int used_df=0;
    NumericMatrix preval = two_runQMCurryTwo<ret_correlation>(x, y, wts, wins, restart_period, min_df, used_df, 
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
                                  bool check_wts=false, bool normalize_wts=true,
                                  bool check_negative_moments=true) {
    int wins=get_wins(window);
    int used_df=0;
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
                                  bool check_wts=false, bool normalize_wts=true,
                                  bool check_negative_moments=true) {
    int wins=get_wins(window);
    int used_df=0;
    NumericMatrix preval = two_runQMCurryTwo<ret_regression_intercept>(x, y, wts, wins, restart_period, min_df, used_df, 
                                                                       na_rm, check_wts, normalize_wts, check_negative_moments);
    return preval;
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
