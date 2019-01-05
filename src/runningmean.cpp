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

  Created: 2018.12.27
  Copyright: Steven E. Pav, 2016-2018
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#include "runningmean.h"

#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' Compute sums or means over a sliding window.
//'
//' @description
//' Compute the mean or sum over 
//' an infinite or finite sliding window, returning a vector the same size as the input.
//' 
//' @param v a vector.
//' @param window the window size. if given as finite integer or double, passed through.
//' If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, equivalent
//' to an infinite window size. If negative, an error will be thrown.
//' @param restart_period the recompute period. because subtraction of elements can cause
//' loss of precision, the computation of moments is restarted periodically based on 
//' this parameter. Larger values mean fewer restarts and faster, though potentially less 
//' accurate results. Unlike in the computation of even order moments, loss of precision
//' is unlikely to be disastrous, so the default value is rather large.
//' @param na_rm whether to remove NA, false by default.
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
//' of \eqn{x_{i-window+1},x_{i-window+2},...,x_{i}}{x_(i-window+1),x_(i-window+2),...,x_i}.
//' Barring \code{NA} or \code{NaN}, this is over a window of size \code{window}.
//' During the 'burn-in' phase, we take fewer elements. If fewer than \code{min_df} for
//' \code{running_mean}, returns \code{NA}.
//'
//' @return A vector the same size as the input.
//' @examples
//' x <- rnorm(1e5)
//' xs <- running_sum(x,10)
//' xm <- running_mean(x,100)
//'
//' @template etc
//' @template ref-romo
//' @template ref-kahan
//' @inheritParams sd3
//' @rdname runningmean 
//' @export
// [[Rcpp::export]]
SEXP running_sum(SEXP v, 
                 SEXP window = R_NilValue, 
                 SEXP wts = R_NilValue,
                 bool na_rm=false, int restart_period=10000,
                 bool check_wts=false) {
    int wins=get_wins(window);
    return runningSumishCurryFour<ret_sum>(v,wts,wins,0,restart_period,na_rm,check_wts);
}

//' @rdname runningmean
//' @export
// [[Rcpp::export]]
SEXP running_mean(SEXP v, 
                  SEXP window = R_NilValue, 
                  SEXP wts = R_NilValue,
                  bool na_rm=false, int min_df=0, int restart_period=10000,
                  bool check_wts=false) {

    // 2FIX: accept NULL restart period and propagate that forward:
    int wins=get_wins(window);
    // turn min_df = 0 into min_df = 1 for later?
    // 2FIX: if you have weights, want a different value. ugh.
    if (min_df < 1) { min_df = 1; }
    return runningSumishCurryFour<ret_mean>(v,wts,wins,min_df,restart_period,na_rm,check_wts);
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
