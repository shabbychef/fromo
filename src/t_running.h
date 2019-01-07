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

  Created: 2019.01.03
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_T_RUNNING__
#define __DEF_T_RUNNING__

#include "common.h"
#include "kahan.h"
#include "welford.h"
#include "runningmean.h"

#include <Rcpp.h>
using namespace Rcpp;

// running sums, moments, cumulants, approximate quantiles

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

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool ord_beyond,bool renormalize,bool na_rm>
NumericMatrix t_runQM(T v,
                      W wts,
                      Rcpp::Nullable< Rcpp::NumericVector > opt_time,
                      Rcpp::Nullable< Rcpp::NumericVector > opt_time_deltas,
                      Rcpp::Nullable< Rcpp::NumericVector > opt_lb_time,
                      const int ord,
                      const double window,
                      const int recom_period,
                      const double lookahead,
                      const int min_df,
                      const double used_df,
                      const bool check_wts,
                      const bool variable_win,
                      const bool wts_as_delta,
                      const bool normalize_wts) {

    Welford<oneW,has_wts,ord_beyond,na_rm> frets = Welford<oneW,has_wts,ord_beyond,na_rm>(ord);

    NumericVector time, time_deltas, lb_time;
    if (opt_time.isNotNull()) {
        time = opt_time.get();
        if (opt_time_deltas.isNotNull()) {
            Rcpp::warning("time deltas given, but not needed; ignoring.");
        }
        if (has_decrease<NumericVector>(time)) { stop("decreasing time detected"); }
    } else {
        if (opt_time_deltas.isNotNull()) {
            time_deltas = opt_time_deltas.get();
        } else {
            if (wts_as_delta) {
                if (has_wts) {
                    time_deltas = wts;
                } else {
                    stop("cannot infer times, as time, time_deltas and weights not given."); // nocov
                }
            } else {
                stop("cannot infer times, as time and time_deltas not given, and wts_as_delta is FALSE."); // nocov
            }
        }
        // to be sure, check again; this might be redundant in the case where deltas are weights, but whatever.
        if (bad_weights<NumericVector>(time_deltas)) { stop("negative time deltas detected"); }
        // just going to use the sugar function here;
        //time = Rcpp::cumsum(time_deltas);
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
        stop("size of time does not match v"); // nocov
    }
    const int numlb = lb_time.size();

    double nextv, prevv, nextw, prevw;

    if (has_wts) {
        if (wts.size() < numel) { stop("size of wts does not match v"); }
    }

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = NumericVector::is_na(window);
    if ((window <= 0) && (!infwin)) { stop("must give positive window"); }
    if (variable_win && !infwin) { Rcpp::warning("variable_win specified, but not being used as a non-na window is given."); }

    // whether to use the gap between lb_time as the effective window
    const bool gapwin = variable_win && infwin;

    if (min_df < 0) { stop("require positive min_df"); }
    // min_df is now on # of observations, but window is a time delta
    // so this makes no sense:
    // if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

    if ((((retwhat==ret_scaled) || 
          (retwhat==ret_zscore) || 
          (retwhat==ret_sharpe) || 
          (retwhat==ret_tstat) || 
          (retwhat==ret_stdev) || 
          (retwhat==ret_sd3)) && (ord < 2)) ||
        (((retwhat==ret_skew) ||
          (retwhat==ret_skew4)) && (ord < 3)) ||
        (((retwhat==ret_sharpese) || 
          (retwhat==ret_exkurt) ||
          (retwhat==ret_exkurt5)) && (ord < 4))) { 
        stop("bad code: order too small to support this computation"); 
    }
    int iii,jjj,lll,tr_iii,tr_jjj;
    // these define the time window for any lll; 
    // we aim to perform computations over the half-open interval
    // (t0,tf]
    // generally we will have
    // t0 = lb_time[lll] + lookahead - window
    // tf = lb_time[lll] + lookahead
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
    double sg_denom,renorm,denom,sigmasq,sigma,sigmapow,mydf,dwsum,skew,exkurt,sr;
    int mmm;

    // preallocated with zeros; should
    // probably be NA?
    int ncols;
    if ((retwhat==ret_centmoments) ||
        (retwhat==ret_stdmoments) ||
        (retwhat==ret_sd3) || 
        (retwhat==ret_skew4) ||
        (retwhat==ret_exkurt5)) {
        ncols = 1+ord; 
    } else if ((retwhat==ret_sharpese)) {
        ncols = 2; 
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
        prev_tf = MIN(tminf,lb_time[0] + lookahead - window - 1.0);  // make it less than t0 will be.
    }

    // now run through lll index//FOLDUP
    for (lll=0;lll < numlb;++lll) {
        tf = lb_time[lll] + lookahead;
        if (gapwin) {
            if (lll==0) {
                t0 = tminf;  // effectively -inf? should there be lookahead? ack.
            } else {
                t0 = lb_time[lll-1] + lookahead;
            }
        } else if (!infwin) {
            t0 = tf - window;
        }
        // otherwise t0 was set as tminf previously.

        // if there is no overlap, then just restart the whole thingy.
        if ((prev_tf <= t0) || (frets.subcount() >= recom_period)) {
            // could bisect, but lets not get fancy
            if (!infwin) { while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { tr_jjj++; } }
            tr_iii = tr_jjj;
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { tr_iii++; }

            frets = quasiWeightedThing<T,W,oneW,has_wts,ord_beyond,na_rm>(v,wts,ord,
                                                                          tr_jjj,       //bottom
                                                                          tr_iii,       //top
                                                                          false);    //no need to check weights as we have done it once above.
        } else {
            if (!infwin) {
                while ((tr_iii < numel) && (time[tr_iii] <= tf) && (time[tr_jjj] <= t0)) { 
                    nextv = double(v[tr_iii]);
                    prevv = double(v[tr_jjj]);
                    if (has_wts) { 
                        nextw = double(wts[tr_iii]); 
                        prevw = double(wts[tr_jjj]); 
                    } 
                    frets.swap_one(nextv,nextw,prevv,prevw); 
                    tr_iii++; 
                    tr_jjj++; 
                }
            }
            // 2FIX: check for subcount? 
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { 
                nextv = double(v[tr_iii]);
                if (has_wts) { nextw = double(wts[tr_iii]); } 
                frets.add_one(nextv,nextw); 
                tr_iii++; 
            }
            if (!infwin) {
                while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { 
                    prevv = double(v[tr_jjj]);
                    if (has_wts) { prevw = double(wts[tr_jjj]); }
                    frets.rem_one(prevv,prevw); 
                    tr_jjj++; 
                }
            }
            // may need to recompute based on the number of subtractions. bummer.
            if (frets.subcount() >= recom_period) {
                frets = quasiWeightedThing<T,W,oneW,has_wts,ord_beyond,na_rm>(v,wts,ord,
                                                                              tr_jjj,       //bottom
                                                                              tr_iii,       //top
                                                                              false);    //no need to check weights as we have done it once above.
            }
        }

        // fill in the value in the output.
        // 2FIX: give access to v, not v[lll]...
        // moment_converter<retwhat, Welford<oneW,has_wts,ord_beyond,na_rm> ,T,renormalize>::mom_interp(xret,lll,ord,frets,v,used_df,min_df);
//yuck!!
#include "moment_interp.h"

        prev_tf = tf;
    }//UNFOLD
    return xret;
}

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool ord_beyond>
NumericMatrix t_runQMCurryZero(T v, 
                               W wts,
                               Rcpp::Nullable< Rcpp::NumericVector > time,
                               Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                               Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                               const int ord,
                               const double window,
                               const int recom_period,
                               const double lookahead,
                               const int min_df,
                               const double used_df,
                               const bool na_rm,
                               const bool check_wts,
                               const bool variable_win,
                               const bool wts_as_delta,
                               const bool normalize_wts) {
    if (has_wts && normalize_wts) {
        if (na_rm) {
            return t_runQM<T,retwhat,W,oneW,has_wts,ord_beyond,true,true>(v, wts, 
                                                                          time, time_deltas, lb_time,
                                                                          ord, window, recom_period, lookahead, min_df, used_df, check_wts, 
                                                                          variable_win, wts_as_delta, normalize_wts); 
        } else {
            return t_runQM<T,retwhat,W,oneW,has_wts,ord_beyond,true,false>(v, wts, 
                                                                          time, time_deltas, lb_time,
                                                                          ord, window, recom_period, lookahead, min_df, used_df, check_wts, 
                                                                          variable_win, wts_as_delta, normalize_wts); 
        }
    } 
    if (na_rm) {
        return t_runQM<T,retwhat,W,oneW,has_wts,ord_beyond,false,true>(v, wts, 
                                                                       time, time_deltas, lb_time,
                                                                       ord, window, recom_period, lookahead, min_df, used_df, check_wts, 
                                                                       variable_win, wts_as_delta, normalize_wts); 
    } 
    return t_runQM<T,retwhat,W,oneW,has_wts,ord_beyond,false,false>(v, wts, 
                                                                       time, time_deltas, lb_time,
                                                                       ord, window, recom_period, lookahead, min_df, used_df, check_wts, 
                                                                       variable_win, wts_as_delta, normalize_wts); 
}

template <typename T,ReturnWhat retwhat,bool ord_beyond>
NumericMatrix t_runQMCurryOne(T v, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts,
                              Rcpp::Nullable< Rcpp::NumericVector > time,
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                              const int ord,
                              const double window,
                              const int recom_period,
                              const double lookahead,
                              const int min_df,
                              const double used_df,
                              const bool na_rm,
                              const bool check_wts,
                              const bool variable_win,
                              const bool wts_as_delta,
                              const bool normalize_wts) {

    //2FIX: typeof wts?
    if (wts.isNotNull()) {
        return t_runQMCurryZero<T,retwhat,NumericVector,double,true,ord_beyond>(v, wts.get(), 
                                                                                time, time_deltas, lb_time,
                                                                                ord, window, recom_period, lookahead, 
                                                                                min_df, used_df, na_rm, check_wts, 
                                                                                variable_win, wts_as_delta, normalize_wts); 
    }
    NumericVector dummy_wts;
    return t_runQMCurryZero<T,retwhat,NumericVector,double,false,ord_beyond>(v, dummy_wts, 
                                                                             time, time_deltas, lb_time,
                                                                             ord, window, recom_period, lookahead, 
                                                                             min_df, used_df, na_rm, check_wts, 
                                                                             variable_win, wts_as_delta, normalize_wts); 
}



template <typename T,ReturnWhat retwhat>
NumericMatrix t_runQMCurryTwo(T v, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts,
                              Rcpp::Nullable< Rcpp::NumericVector > time,
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                              const int ord,
                              const double window,
                              const int recom_period,
                              const double lookahead,
                              const int min_df,
                              const double used_df,
                              const bool na_rm,
                              const bool check_wts,
                              const bool variable_win,
                              const bool wts_as_delta,
                              const bool normalize_wts) {

    if (ord==2) {
        return t_runQMCurryOne<T,retwhat,false>(v, wts, 
                                                time, time_deltas, lb_time,
                                                ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, 
                                                variable_win, wts_as_delta, normalize_wts); 
    }
    return t_runQMCurryOne<T,retwhat,true>(v, wts, 
                                           time, time_deltas, lb_time,
                                           ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, 
                                           variable_win, wts_as_delta, normalize_wts); 
}

template <ReturnWhat retwhat>
NumericMatrix t_runQMCurryThree(SEXP v, 
                                Rcpp::Nullable< Rcpp::NumericVector > wts,
                                Rcpp::Nullable< Rcpp::NumericVector > time,
                                Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                                Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                                const int ord,
                                const double window,
                                const int recom_period,
                                const double lookahead,
                                const int min_df,
                                const double used_df,
                                const bool na_rm,
                                const bool check_wts,
                                const bool variable_win,
                                const bool wts_as_delta,
                                const bool normalize_wts) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return t_runQMCurryTwo<IntegerVector,retwhat>(v, wts, 
                                                                      time, time_deltas, lb_time,
                                                                      ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, 
                                                                      variable_win, wts_as_delta, normalize_wts); }
        case REALSXP: { return t_runQMCurryTwo<NumericVector,retwhat>(v, wts, 
                                                                      time, time_deltas, lb_time,
                                                                      ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, 
                                                                      variable_win, wts_as_delta, normalize_wts); }
        // to make smaller binaries, and because who cares about logicals, I convert them to integers here...
        case  LGLSXP: { return t_runQMCurryTwo<IntegerVector,retwhat>(as<IntegerVector>(v), wts,  
                                                                      time, time_deltas, lb_time,
                                                                      ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, 
                                                                      variable_win, wts_as_delta, normalize_wts); }
        default: stop("Unsupported weight type"); // nocov
    }
    // have to have fallthrough for CRAN check.
    NumericMatrix retv;
    return retv;
}

#endif /* __DEF_T_RUNNING__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
