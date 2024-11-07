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

  time-based running operations.

  Created: 2019.01.03
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_RUNCODE__
#define __DEF_RUNCODE__

#include "common.h"
#include "kahan.h"
#include "welford.h"
#include "running.h"
#include "t_running.h"

#endif /* __DEF_RUNCODE__ */

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// running sums, moments, cumulants, approximate quantiles

//' @title
//' Compute first K moments over a sliding time-based window
//' @description
//' Compute the (standardized) 2nd through kth moments, the mean, and the number of elements over
//' an infinite or finite sliding time based window, returning a matrix. 
//' 
//' @param v a vector of data.
//' @param time  an optional vector of the timestamps of \code{v}. If given, must be
//'  the same length as \code{v}. If not given, we try to infer it by summing the
//'  \code{time_deltas}.
//' @param time_deltas  an optional vector of the deltas of timestamps. If given, must be
//'  the same length as \code{v}. If not given, and \code{wts} are given and \code{wts_as_delta} is true,
//'  we take the \code{wts} as the time deltas.  The deltas must be positive. We sum them to arrive
//'  at the times.
//' @param window the window size, in time units. if given as finite integer or double, passed through.
//' If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, 
//'  and \code{variable_win} is true, then we infer the window from the lookback times: the
//'  first window is infinite, but the remaining is the deltas between lookback times.
//'  If \code{variable_win} is false, then these undefined values are equivalent to an
//'  infinite window.
//'  If negative, an error will be thrown.
//' @param lb_time  a vector of the times from which lookback will be performed. The output should
//'  be the same size as this vector. If not given, defaults to \code{time}.
//' @param na_rm whether to remove NA, false by default.
//' @param max_order the maximum order of the centered moment to be computed.
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent moments from being computed on too few observations.
//' Defaults to zero, meaning no restriction.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the centered moments computation. These are subtracted from the number of
//' observations.
//' @param variable_win  if true, and the \code{window} is not a concrete number,
//'  the computation window becomes the time between lookback times.
//' @param wts_as_delta  if true and the \code{time} and \code{time_deltas} are not
//' given, but \code{wts} are given, we take \code{wts} as the \code{time_deltas}.
//'
//' @details
//'
//' Computes the number of elements, the mean, and the 2nd through kth
//' centered (and typically standardized) moments, for \eqn{k=2,3,4}{k=2,3,4}. These
//' are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//'
//' Given the length \eqn{n} vector \eqn{x}, we output matrix \eqn{M} where
//' \eqn{M_{i,j}}{M_i,j} is the \eqn{order - j + 1} moment (\emph{i.e.}
//' excess kurtosis, skewness, standard deviation, mean or number of elements)
//' of some elements \eqn{x_i} defined by the sliding time window.
//' Barring \code{NA} or \code{NaN}, this is over a window of time width \code{window}.
//'
//' @template sec-t-win
//'
//' @return Typically a matrix, where the first columns are the kth, k-1th through 2nd standardized, 
//' centered moments, then a column of the mean, then a column of the number of (non-nan) elements in the input,
//' with the following exceptions:
//' \describe{
//' \item{t_running_cent_moments}{Computes arbitrary order centered moments. When \code{max_order_only} is set,
//' only a column of the maximum order centered moment is returned.}
//' \item{t_running_std_moments}{Computes arbitrary order standardized moments, then the standard deviation, the mean,
//' and the count. There is not yet an option for \code{max_order_only}, but probably should be.}
//' \item{t_running_cumulants}{Computes arbitrary order cumulants, and returns the kth, k-1th, through the second 
//' (which is the variance) cumulant, then the mean, and the count.}
//' }
//'
//' @note
//' the kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @examples
//' x <- rnorm(1e5)
//' xs3 <- t_running_sd3(x,time=seq_along(x),window=10)
//' xs4 <- t_running_skew4(x,time=seq_along(x),window=10)
//' # but what if you only cared about some middle values?
//' xs4 <- t_running_skew4(x,time=seq_along(x),lb_time=(length(x) / 2) + 0:10,window=20)
//'
//' @inheritParams running_sd3
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @template param-heywood
//' @template note-wts
//' @template note-heywood
//' @seealso \code{\link{running_sd3}}.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_sd3(SEXP v, 
                            Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                            SEXP window = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                            bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                            bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true, 
                            bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    NumericMatrix preval = t_runQMCurryThree<ret_sd3>(v, wts, time, time_deltas, lb_time,
                                                      2, wins, restart_period, lookahead, min_df, used_df, 
                                                      na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
    return preval;
}
// return the skew, the standard deviation, the mean, and the dof
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_skew4(SEXP v, 
                              Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                              SEXP window = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                              bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                              bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                              bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    NumericMatrix preval = t_runQMCurryThree<ret_skew4>(v, wts, time, time_deltas, lb_time,
                                                        3, wins, restart_period, lookahead, min_df, used_df, 
                                                        na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);

    return preval;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_kurt5(SEXP v, 
                              Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                              SEXP window = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                              bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                              bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                              bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    NumericMatrix preval = t_runQMCurryThree<ret_exkurt5>(v, wts, time, time_deltas, lb_time,
                                                          4, wins, restart_period, lookahead, min_df, used_df, 
                                                          na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
    return preval;
}
// just the sd nothing else.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_sd(SEXP v, 
                           Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                           Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                           SEXP window = R_NilValue, 
                           Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                           Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                           bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                           bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                           bool check_negative_moments=true) {
//2FIX: introduce used_df ... 
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    return t_runQMCurryThree<ret_stdev>(v, wts, time, time_deltas, lb_time,
                                        2, wins, restart_period, lookahead, min_df, used_df, 
                                        na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// just the skew nothing else.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_skew(SEXP v, 
                             Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                             SEXP window = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                             bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                             bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                             bool check_negative_moments=true) {
//2FIX: introduce used_df ... 
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    return t_runQMCurryThree<ret_skew>(v, wts, time, time_deltas, lb_time,
                                       3, wins, restart_period, lookahead, min_df, used_df, 
                                       na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// just the kurtosis nothing else.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_kurt(SEXP v, 
                             Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                             SEXP window = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                             bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                             bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                             bool check_negative_moments=true) {
//2FIX: introduce used_df ... 
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    return t_runQMCurryThree<ret_exkurt>(v, wts, time, time_deltas, lb_time,
                                         4, wins, restart_period, lookahead, min_df, used_df, 
                                         na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}

// return the centered moments down to the 2nd, then the mean, and the dof.
//' @param max_order_only for \code{running_cent_moments}, if this flag is set, only compute
//' the maximum order centered moment, and return in a vector.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_cent_moments(SEXP v, 
                                     Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                     Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                     SEXP window = R_NilValue, 
                                     Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                     Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                                     int max_order=5, bool na_rm=false, bool max_order_only=false, 
                                     int min_df=0, double used_df=0.0, int restart_period=100, 
                                     bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                     bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    if (max_order_only) {
        return t_runQMCurryThree<ret_centmaxonly>(v, wts, time, time_deltas, lb_time,
                                         max_order, wins, restart_period, lookahead, min_df, used_df, 
                                         na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
    } 
    return t_runQMCurryThree<ret_centmoments>(v, wts, time, time_deltas, lb_time,
                                              max_order, wins, restart_period, lookahead, min_df, used_df, 
                                              na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}

// return the standardized moments down to the 3rd, then the standard deviation, the mean, and the dof.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_std_moments(SEXP v, 
                                    Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                    Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                    SEXP window = R_NilValue, 
                                    Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                    Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                                    int max_order=5, bool na_rm=false, 
                                    int min_df=0, double used_df=0.0, int restart_period=100, 
                                    bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                    bool check_negative_moments=true) {


    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    return t_runQMCurryThree<ret_stdmoments>(v, wts, time, time_deltas, lb_time,
                                             max_order, wins, restart_period, lookahead, min_df, used_df, 
                                             na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}

// return the cumulants down to the 2nd, then the standard deviation, the mean, and the dof.
//' @rdname t_runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_cumulants(SEXP v, 
                                  Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                  Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                  SEXP window = R_NilValue, 
                                  Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                  Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                                  int max_order=5, bool na_rm=false, 
                                  int min_df=0, double used_df=0.0, int restart_period=100, 
                                  bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                  bool check_negative_moments=true) {

    NumericMatrix cumulants = t_running_cent_moments(v, time, time_deltas, window, wts, lb_time, max_order, na_rm, false,
                                                     min_df, used_df, restart_period, variable_win, wts_as_delta, check_wts, normalize_wts, check_negative_moments);


    // changes the cumulants in the background
    centmom2cumulants(cumulants, max_order);
    return cumulants;
}
//' @title
//' Compute approximate quantiles over a sliding time window
//' @description
//' Computes cumulants up to some given order, then employs the Cornish-Fisher approximation
//' to compute approximate quantiles using a Gaussian basis.
//' 
//' @param p the probability points at which to compute the quantiles. Should be in the range (0,1).
//' @inheritParams t_running_cumulants
//'
//' @details
//'
//' Computes the cumulants, then approximates quantiles using AS269 of Lee & Lin.
//' @template sec-t-win
//'
//' @return A matrix, with one row for each element of \code{x}, and one column for each element of \code{q}.
//'
//' @note
//' The current implementation is not as space-efficient as it could be, as it first computes
//' the cumulants for each row, then performs the Cornish-Fisher approximation on a row-by-row
//' basis. In the future, this computation may be moved earlier into the pipeline to be more
//' space efficient. File an issue if the memory footprint is an issue for you.
//'
//' @examples
//' x <- rnorm(1e5)
//' xq <- t_running_apx_quantiles(x,c(0.1,0.25,0.5,0.75,0.9),
//'        time=seq_along(x),window=200,lb_time=c(100,200,400))
//'
//' xq <- t_running_apx_median(x,time=seq_along(x),window=200,lb_time=c(100,200,400))
//' xq <- t_running_apx_median(x,time=cumsum(runif(length(x),min=0.5,max=1.5)),
//'       window=200,lb_time=c(100,200,400))
//'
//' # weighted median?
//' wts <- runif(length(x),min=1,max=5)
//' xq <- t_running_apx_median(x,wts=wts,wts_as_delta=TRUE,window=1000,lb_time=seq(1000,10000,by=1000))
//'
//' # these should give the same answer:
//' xr <- running_apx_median(x,window=200);
//' xt <- t_running_apx_median(x,time=seq_along(x),window=199.99)
//'
//' @seealso \code{\link{running_apx_quantiles}}, \code{\link{t_running_cumulants}}, \code{PDQutils::qapx_cf}, \code{PDQutils::AS269}.
//' @template etc
//' @template ref-cf
//' @template ref-romo
//' @template param-wts
//' @template param-heywood
//' @template note-wts
//' @rdname t_runningquantiles
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_apx_quantiles(SEXP v, NumericVector p, 
                                      Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                      Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                      SEXP window = R_NilValue, 
                                      Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                      Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                                      int max_order=5, bool na_rm=false, 
                                      int min_df=0, double used_df=0.0, int restart_period=100, 
                                      bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                      bool check_negative_moments=true) {
                                    
    NumericMatrix cumulants = t_running_cumulants(v, time, time_deltas, window, wts, lb_time,
                                                  max_order, na_rm, min_df, used_df, restart_period, 
                                                  variable_win, wts_as_delta, check_wts, normalize_wts, check_negative_moments);

    return cumulants2quantiles(cumulants, p, max_order);
}
//' @rdname t_runningquantiles
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_apx_median(SEXP v, 
                                   Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                   Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                   SEXP window = R_NilValue, 
                                   Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                   Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                                   int max_order=5, bool na_rm=false, 
                                   int min_df=0, double used_df=0.0, int restart_period=100, 
                                   bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                   bool check_negative_moments=true) {
                                   
    NumericVector p(1);
    p(0) = 0.5;
    return t_running_apx_quantiles(v,p, time, time_deltas, window, wts, lb_time,
                                   max_order, na_rm, min_df, used_df, restart_period, 
                                   variable_win, wts_as_delta, check_wts, normalize_wts, check_negative_moments);
}

//' @title
//' Compare data to moments computed over a time sliding window.
//' @description
//' Computes moments over a sliding window, then adjusts the data accordingly, centering, or scaling,
//' or z-scoring, and so on.
//' 
//' @inheritParams t_running_cent_moments
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent \emph{e.g.} Z-scores from being computed on only 3
//' observations. Defaults to zero, meaning no restriction, which can result in 
//' infinite Z-scores during the burn-in period.
//' @param lookahead for some of the operations, the value is compared to 
//' mean and standard deviation possibly using 'future' or 'past' information
//' by means of a non-zero lookahead. Positive values mean data are taken from
//' the future. This is in time units, and so should be a real.
//' @param compute_se for \code{running_sharpe}, return an extra column of the
//' standard error, as computed by Mertens' correction.
//'
//' @details
//'
//' Given the length \eqn{n} vector \eqn{x}, for
//' a given index \eqn{i}, define \eqn{x^{(i)}}{x^(i)}
//' as the elements of \eqn{x} defined by the sliding time window (see the section
//' on time windowing).
//' Then define \eqn{\mu_i}{mu_i}, \eqn{\sigma_i}{sigma_i}
//' and \eqn{n_i}{n_i} as, respectively, the sample mean, standard deviation and number of
//' non-NA elements in \eqn{x^{(i)}}{x^(i)}. 
//'
//' We compute output vector \eqn{m} the same size as \eqn{x}. 
//' For the 'centered' version of \eqn{x}, we have \eqn{m_i = x_i - \mu_i}{m_i = x_i - mu_i}.
//' For the 'scaled' version of \eqn{x}, we have \eqn{m_i = x_i / \sigma_i}{m_i = x_i / sigma_i}.
//' For the 'z-scored' version of \eqn{x}, we have \eqn{m_i = (x_i - \mu_i) / \sigma_i}{m_i = (x_i - mu_i) / sigma_i}.
//' For the 't-scored' version of \eqn{x}, we have \eqn{m_i = \sqrt{n_i} \mu_i / \sigma_i}{m_i = sqrt(n_i) mu_i / sigma_i}.
//'
//' We also allow a 'lookahead' for some of these operations.
//' If positive, the moments are computed using data from larger indices;
//' if negative, from smaller indices. 
//'
//' @template sec-t-win
//'
//' @return a vector the same size as the input consisting of the adjusted version of the input.
//' When there are not sufficient (non-nan) elements for the computation, \code{NaN} are returned.
//'
//' @seealso \code{\link{running_centered}}, \code{\link{scale}}
//' @template etc
//' @template ref-romo
//' @template note-wts
//' @rdname t_runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_centered(SEXP v, 
                                 Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                 Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                 SEXP window = R_NilValue, 
                                 Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                 bool na_rm=false, int min_df=0, double used_df=1.0, double lookahead=0.0, int restart_period=100,
                                 bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                 bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue;
    return t_runQMCurryThree<ret_centered>(v, wts, time, time_deltas, lb_time,
                                           1, wins, restart_period, lookahead, min_df, used_df, 
                                           na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// scale the input
//' @rdname t_runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_scaled(SEXP v, 
                               Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                               SEXP window = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                               bool na_rm=false, int min_df=0, double used_df=1.0, double lookahead=0.0, int restart_period=100,
                               bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                               bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue;
    return t_runQMCurryThree<ret_scaled>(v, wts, time, time_deltas, lb_time,
                                         2, wins, restart_period, lookahead, min_df, used_df, 
                                         na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// zscore the input
//' @rdname t_runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_zscored(SEXP v, 
                                Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                                Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                                SEXP window = R_NilValue, 
                                Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                bool na_rm=false, int min_df=0, double used_df=1.0, double lookahead=0.0, int restart_period=100,
                                bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                                bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue;
    return t_runQMCurryThree<ret_zscore>(v, wts, time, time_deltas, lb_time,
                                         2, wins, restart_period, lookahead, min_df, used_df, 
                                         na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// sharpe on the input
//' @rdname t_runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_sharpe(SEXP v, 
                               Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                               SEXP window = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                               bool na_rm=false, bool compute_se=false, int min_df=0, double used_df=1.0, int restart_period=100,
                               bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                               bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    if (compute_se) {
        return t_runQMCurryThree<ret_sharpese>(v, wts, time, time_deltas, lb_time,
                                               4, wins, restart_period, lookahead, min_df, used_df, 
                                               na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
    } 
    return t_runQMCurryThree<ret_sharpe>(v, wts, time, time_deltas, lb_time,
                                         2, wins, restart_period, lookahead, min_df, used_df, 
                                         na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}
// t stat of the input
//' @rdname t_runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix t_running_tstat(SEXP v, 
                              Rcpp::Nullable< Rcpp::NumericVector > time = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas = R_NilValue, 
                              SEXP window = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time = R_NilValue, 
                              bool na_rm=false, bool compute_se=false, int min_df=0, double used_df=1.0, int restart_period=100,
                              bool variable_win=false, bool wts_as_delta=true, bool check_wts=false, bool normalize_wts=true,
                              bool check_negative_moments=true) {
    double wins = get_double_wins(window);
    const double lookahead=0.0; // 2FIX: allow user to set?
    return t_runQMCurryThree<ret_tstat>(v, wts, time, time_deltas, lb_time,
                                        2, wins, restart_period, lookahead, min_df, used_df, 
                                        na_rm, check_wts, variable_win, wts_as_delta, normalize_wts, check_negative_moments);
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
