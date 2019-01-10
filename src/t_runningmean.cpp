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

  running sum and mean functions.

  Created: 2019.01.05
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#include "t_runningmean.h"

#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' Compute sums or means over a sliding time window.
//'
//' @description
//' Compute the mean or sum over 
//' an infinite or finite sliding time window, returning a vector the same size as the lookback
//' times.
//' 
//' @param v a vector.
//' @inheritParams t_running_sd
//'
//' @param restart_period the recompute period. because subtraction of elements can cause
//' loss of precision, the computation of moments is restarted periodically based on 
//' this parameter. Larger values mean fewer restarts and faster, though potentially less 
//' accurate results. Unlike in the computation of even order moments, loss of precision
//' is unlikely to be disastrous, so the default value is rather large.
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned,
//' only for the means computation.
//' This can be used to prevent moments from being computed on too few observations.
//' Defaults to zero, meaning no restriction.
//'
//' @details
//'
//' Computes the mean or sum of the elements, using a Kahan's Compensated Summation Algorithm,
//' a numerically robust one-pass method.
//'
//' Given the length \eqn{n} vector \eqn{x}, we output matrix \eqn{M} where
//' \eqn{M_{i,1}}{M_i,1} is the sum or mean 
//' of some elements \eqn{x_i} defined by the sliding time window.
//' Barring \code{NA} or \code{NaN}, this is over a window of time width \code{window}.
//'
//' @template sec-t-win
//'
//' @return A vector the same size as the lookback times.
//' @examples
//' x <- rnorm(1e5)
//' xs <- t_running_sum(x,time=seq_along(x),window=10)
//' xm <- t_running_mean(x,time=cumsum(runif(length(x))),window=7.3)
//'
//' @template etc
//' @template ref-romo
//' @template ref-kahan
//' @inheritParams sd3
//' @template note-wts
//' @rdname t_runningmean 
//' @export
// [[Rcpp::export]]
SEXP t_running_sum(SEXP v, 
                   Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                   Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                   SEXP window = R_NilValue, 
                   SEXP wts = R_NilValue,
                   Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                   bool na_rm=false, int min_df=0, int restart_period=10000,
                   bool variable_win=false, bool wts_as_delta=true, bool check_wts=false) {
    double wins = get_double_wins(window);
    return t_runningSumishCurryFour<ret_sum>(v,time,time_deltas,wins,wts,lb_time,
                                             na_rm,0,restart_period,
                                             variable_win,wts_as_delta,check_wts);
}

//' @rdname t_runningmean
//' @export
// [[Rcpp::export]]
SEXP t_running_mean(SEXP v, 
                    Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                    Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                    SEXP window = R_NilValue, 
                    SEXP wts = R_NilValue,
                    Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                    bool na_rm=false, int min_df=0, int restart_period=10000,
                    bool variable_win=false, bool wts_as_delta=true, bool check_wts=false) {

    // 2FIX: accept NULL restart period and propagate that forward:
    double wins = get_double_wins(window);
    // turn min_df = 0 into min_df = 1 for later?
    // 2FIX: if you have weights, want a different value. ugh.
    if (min_df < 1) { min_df = 1; }
    return t_runningSumishCurryFour<ret_mean>(v,time,time_deltas,wins,wts,lb_time,
                                              na_rm,min_df,restart_period,
                                              variable_win,wts_as_delta,check_wts);
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
