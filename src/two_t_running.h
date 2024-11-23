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

  Created: 2024.11.23
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_TWO_T_RUNNING__
#define __DEF_TWO_T_RUNNING__

#include "common.h"
#include "kahan.h"
#include "two_welford.h"
#include "runningmean.h"

#include <Rcpp.h>
using namespace Rcpp;

// running covariance, correlation, regression and so on.

// some notes:
// I believe we can just use
// Vector<REALSXP> instead of NumericVector
// Vector<INTSXP> instead of IntegerVector
// Vector<LGLSXP> instead of LogicalVector
// which may simplify some of the W vs oneW nonsense
//
// also there is a std::enable_if thingy which we can
// use, I think, to do boolean template magic.

// ok here is how we compute these;
//
// first figure out the time;
// if time is given, use it. 
// if time is null, figure out time_deltas and sum that to get times.
// if time deltas are null, throw an error.
//
// if time_deltas are given, and needed, use them.
// if time deltas are given, and not needed, throw a warning.
// if time deltas are not given and neede, and wts_as_delta is true,
// and wts are given, then use the weights as the time deltas,
// otherwise throw an error
//
// if lb_time is given, use it as the lookback period.
// if not, use time.

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool na_rm>
NumericMatrix two_t_runQM(T v,
                          T vv,
                          W wts,
                          Rcpp::Nullable< Rcpp::NumericVector > opt_time,
                          Rcpp::Nullable< Rcpp::NumericVector > opt_time_deltas,
                          Rcpp::Nullable< Rcpp::NumericVector > opt_lb_time,
                          // const int ord,
                          const double window,
                          const int recom_period,
                          const int min_df,
                          const double used_df,
                          const bool check_wts,
                          const bool variable_win,
                          const bool wts_as_delta,
                          const bool renormalize,
                          const bool check_negative_moments) {

    // a bit of a hack here, but you must have ord >= 2 for Welford
    // objects, otherwise it hits a memory leak. I know that previously
    // I had hard coded Welford-like objects for the ord=1 case,
    // but that makes huge libraries. instead, just hack this MAX in here
    //const int fake_ord = MAX(ord,2);

    TwoWelford<oneW,has_wts,na_rm> frets = TwoWelford<oneW,has_wts,na_rm>();
    frets.tare();

    NumericVector time, time_deltas, lb_time;
    if (opt_time.isNotNull()) {
        time = opt_time.get();
        if (opt_time_deltas.isNotNull()) {
            Rcpp::warning("time deltas given, but not needed; ignoring."); // #nocov
        }
        if (has_decrease<NumericVector>(time)) { stop("decreasing time detected"); } // #nocov
    } else {
        if (opt_time_deltas.isNotNull()) {
            time_deltas = opt_time_deltas.get();
        } else {
            if (wts_as_delta) {
                if (has_wts) {
                    time_deltas = wts;
                } else {
                    stop("cannot infer times, as time, time_deltas and weights not given."); // #nocov
                }
            } else {
                stop("cannot infer times, as time and time_deltas not given, and wts_as_delta is FALSE."); // #nocov
            }
        }
        // to be sure, check again; this might be redundant in the case where deltas are weights, but whatever.
        if (bad_weights<NumericVector>(time_deltas)) { stop("negative time deltas detected"); }
        // just going to use the sugar function here;
        //time = Rcpp::cumsum(time_deltas);
        // 2FIX: is this all good?
        time = (runningSumishCurryFour<ret_sum>(time_deltas,R_NilValue,NA_INTEGER,0,100000,false,false));
    }
    if (opt_lb_time.isNotNull()) {
        lb_time = opt_lb_time.get();
        if (has_decrease<NumericVector>(lb_time)) { stop("decreasing lb_time detected"); }
    } else {
        // read only so relax.
        lb_time = time;
    }

    const int numel = v.size();
    if (time.size() != numel) {
        stop("size of time does not match v"); // #nocov
    }
    const int numlb = lb_time.size();

    double nextv, prevv, nextvv, prevvv, nextw, prevw;

    if (has_wts) {
        if (wts.size() < numel) { stop("size of wts does not match v"); }
    }
    const int ord = 2;

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = NumericVector::is_na(window);
    if ((window <= 0) && (!infwin)) { stop("must give positive window"); } // #nocov
    if (variable_win && !infwin) { Rcpp::warning("variable_win specified, but not being used as a non-na window is given."); } // #nocov

    // whether to use the gap between lb_time as the effective window
    const bool gapwin = variable_win && infwin;

    if (min_df < 0) { stop("require positive min_df"); }
    // min_df is now on # of observations, but window is a time delta
    // so this makes no sense:
    // if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

    if (!(
          (retwhat==ret_correlation) ||
          (retwhat==ret_covariance) ||
          (retwhat==ret_covariance_matrix) ||
          (retwhat==ret_regression_slope) ||
          (retwhat==ret_regression_intercept) ||
          (retwhat==ret_regression_fit) ||
          (retwhat==ret_regression_diagnostics)
        )) {
        stop("NYI: only understand correlation, covariance, etc");  // #nocov
    }
    int iii,jjj,lll,tr_iii,tr_jjj;
    // these define the time window for any lll; 
    // we aim to perform computations over the half-open interval
    // (t0,tf]
    // generally we will have
    // t0 = lb_time[lll] + lookahead - window
    // tf = lb_time[lll] + lookahead
    // though here we have lookahead = 0 as it is not well defined in the x&y case.
    // note that the the window can depend on lll in a simple way:
    // if gapwin is true, then
    // t0 = lb_time[lll-1] + lookahead
    // tf = lb_time[lll] + lookahead
    // where lb_time[-1] is understood to be -inf
    double tf,t0;
    double prev_tf;
    const double tminf = time[0] - 1.0;  // effectively -inf?
    if (!gapwin && infwin) { t0 = tminf; }

    // super gross; I need these for the include later.
    double denom;
    int mmm;

    // preallocated with zeros; should
    // probably be NA?
    int ncols;
    if ((retwhat==ret_covariance_matrix)) {
        ncols = 3;
    } else if (retwhat==ret_regression_fit) {
        ncols = 2;
    } else if (retwhat==ret_regression_diagnostics) {
        ncols = 5;
    } else {
        ncols = 1; 
    }
        
    NumericMatrix xret(numlb,ncols);

    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
    }
    // set them once only
    if (!has_wts) {
        nextw = 1.0;
        prevw = 1.0;
    }

    tr_jjj = 0;
    tr_iii = -1;
    if (infwin) {
        prev_tf = tminf;
    } else {
        prev_tf = MIN(tminf,lb_time[0] - window - 1.0);  // make it less than t0 will be.
    }

    // now run through lll index//FOLDUP
    for (lll=0;lll < numlb;++lll) {
        tf = lb_time[lll];
        if (gapwin) {
            if (lll==0) {
                t0 = tminf;  // effectively -inf? 
            } else {
                t0 = lb_time[lll-1];
            }
        } else if (!infwin) {
            t0 = tf - window;
        }
        // otherwise t0 was set as tminf previously.

        // if there is no overlap, then just restart the whole thingy.
        if ((prev_tf <= t0) || (frets.subcount() >= recom_period)) {
            // could bisect, but lets not get fancy
            if (!infwin || gapwin) { while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { tr_jjj++; } }
            tr_iii = tr_jjj;
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { tr_iii++; }
            //zero it out
            frets.tare();
            add_many<T,W,oneW,has_wts,na_rm>(frets,
                                             v,vv,wts,
                                             tr_jjj,  //bottom
                                             tr_iii,  //top
                                             false);  //no need to check weights as we have done it once above.
        } else {
            if (!infwin || gapwin) {
                while ((tr_iii < numel) && (time[tr_iii] <= tf) && (time[tr_jjj] <= t0)) { 
                    nextv = double(v[tr_iii]);
                    nextvv = double(vv[tr_iii]);
                    prevv = double(v[tr_jjj]);
                    prevvv = double(vv[tr_jjj]);
                    if (has_wts) { 
                        nextw = double(wts[tr_iii]); 
                        prevw = double(wts[tr_jjj]); 
                    } 
                    frets.swap_one(nextv,nextvv,nextw,
                                   prevv,prevvv,prevw); 
                    tr_iii++; 
                    tr_jjj++; 
                }
            }
            // 2FIX: check for subcount? 
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { 
                nextv = double(v[tr_iii]);
                nextvv = double(vv[tr_iii]);
                if (has_wts) { nextw = double(wts[tr_iii]); } 
                frets.add_one(nextv,nextvv,nextw); 
                tr_iii++; 
            }
            if (!infwin || gapwin) {
                while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { 
                    prevv = double(v[tr_jjj]);
                    prevvv = double(vv[tr_jjj]);
                    if (has_wts) { prevw = double(wts[tr_jjj]); }
                    frets.rem_one(prevv,prevvv,prevw); 
                    tr_jjj++; 
                }
            }
            // may need to recompute based on the number of subtractions. bummer.
            if ((frets.subcount() >= recom_period) || (check_negative_moments && frets.has_heywood())) {
                //zero it out
                frets.tare();
                add_many<T,W,oneW,has_wts,na_rm>(frets,
                                                 v,vv,wts,
                                                 tr_jjj,  //bottom
                                                 tr_iii,  //top
                                                 false);  //no need to check weights as we have done it once above.
            }
        }

        // fill in the value in the output.
        // 2FIX: give access to v, not v[lll]...
//yuck!!
#include "two_moment_interp.h"

        prev_tf = tf;
    }//UNFOLD
    return xret;
}

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts>
NumericMatrix two_t_runQMCurryZero(T v, T vv,
                                   W wts,
                                   Rcpp::Nullable< Rcpp::NumericVector > time,
                                   Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                                   Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                                   // const int ord,
                                   const double window,
                                   const int recom_period,
                                   const int min_df,
                                   const double used_df,
                                   const bool na_rm,
                                   const bool check_wts,
                                   const bool variable_win,
                                   const bool wts_as_delta,
                                   const bool normalize_wts,
                                   const bool check_negative_moments) {
    if (na_rm) {
        return two_t_runQM<T,retwhat,W,oneW,has_wts,true>(v, vv, wts, 
                                                          time, time_deltas, lb_time,
                                                          window, recom_period, min_df, used_df, check_wts, 
                                                          variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
    } 
    return two_t_runQM<T,retwhat,W,oneW,has_wts,false>(v, vv, wts, 
                                                       time, time_deltas, lb_time,
                                                       window, recom_period, min_df, used_df, check_wts, 
                                                       variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
}

template <typename T,ReturnWhat retwhat>
NumericMatrix two_t_runQMCurryOne(T v, T vv,
                                  Rcpp::Nullable< Rcpp::NumericVector > wts,
                                  Rcpp::Nullable< Rcpp::NumericVector > time,
                                  Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                                  Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                                  // const int ord,
                                  const double window,
                                  const int recom_period,
                                  const int min_df,
                                  const double used_df,
                                  const bool na_rm,
                                  const bool check_wts,
                                  const bool variable_win,
                                  const bool wts_as_delta,
                                  const bool normalize_wts,
                                  const bool check_negative_moments) {

    //2FIX: typeof wts?
    if (wts.isNotNull()) {
        return two_t_runQMCurryZero<T,retwhat,NumericVector,double,true>(v, vv, wts.get(), 
                                                                         time, time_deltas, lb_time,
                                                                         window, recom_period, 
                                                                         min_df, used_df, na_rm, check_wts, 
                                                                         variable_win, wts_as_delta, normalize_wts, 
                                                                         check_negative_moments); 
    }
    NumericVector dummy_wts;
    return two_t_runQMCurryZero<T,retwhat,NumericVector,double,false>(v, vv, dummy_wts, 
                                                                      time, time_deltas, lb_time,
                                                                      window, recom_period, 
                                                                      min_df, used_df, na_rm, check_wts, 
                                                                      variable_win, wts_as_delta, normalize_wts,
                                                                      check_negative_moments);
}


template <ReturnWhat retwhat>
NumericMatrix two_t_runQMCurryTwo(SEXP v, SEXP vv,
                                  Rcpp::Nullable< Rcpp::NumericVector > wts,
                                  Rcpp::Nullable< Rcpp::NumericVector > time,
                                  Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                                  Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                                  // const int ord,
                                  const double window,
                                  const int recom_period,
                                  const int min_df,
                                  const double used_df,
                                  const bool na_rm,
                                  const bool check_wts,
                                  const bool variable_win,
                                  const bool wts_as_delta,
                                  const bool normalize_wts,
                                  const bool check_negative_moments) {
    switch (TYPEOF(v)) {
        // 2FIX: we might have to cast vv to the same type as v?
        case  INTSXP: 
            { 
                switch (TYPEOF(vv)) {
                    case  INTSXP: 
                        { 
                            return two_t_runQMCurryOne<IntegerVector,retwhat>(v, vv, wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    case REALSXP: // cast int to numeric
                        {  
                            return two_t_runQMCurryOne<NumericVector,retwhat>(as<NumericVector>(v), vv, wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to int
                        {
                            return two_t_runQMCurryOne<IntegerVector,retwhat>(v, as<IntegerVector>(vv), wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    default: stop("Unsupported data type for vv"); // #nocov
                }
            }
        case REALSXP: 
            {
                switch (TYPEOF(vv)) {
                    case  INTSXP: // cast int to numeric
                        { 
                            return two_t_runQMCurryOne<NumericVector,retwhat>(v, as<NumericVector>(vv), wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        } 
                    case REALSXP: 
                        {  
                            return two_t_runQMCurryOne<NumericVector,retwhat>(v, vv, wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to numeric
                        {
                            return two_t_runQMCurryOne<NumericVector,retwhat>(v, as<NumericVector>(vv), wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments);
                        }
                    default: stop("Unsupported data type for vv"); // #nocov
                }
            }
        case  LGLSXP: 
            {
                switch (TYPEOF(vv)) {
                    case  INTSXP: // cast logical to int 
                        { 
                            return two_t_runQMCurryOne<IntegerVector,retwhat>(as<IntegerVector>(v), vv, wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        } 
                    case REALSXP: // cast logical to numeric
                        {  
                            return two_t_runQMCurryOne<NumericVector,retwhat>(as<NumericVector>(v), vv, wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to integer
                        {
                            return two_t_runQMCurryOne<IntegerVector,retwhat>(as<IntegerVector>(v), as<IntegerVector>(vv), wts, 
                                                                              time, time_deltas, lb_time,
                                                                              window, recom_period, min_df, used_df, na_rm, check_wts, 
                                                                              variable_win, wts_as_delta, normalize_wts, check_negative_moments); 
                        }
                    default: stop("Unsupported data type for vv"); // #nocov
                }
            }
        default: stop("Unsupported data type for v"); // #nocov
    }
    // have to have fallthrough for CRAN check.
    NumericMatrix retv;
    return retv;
}

#endif /* __DEF_TWO_T_RUNNING__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
