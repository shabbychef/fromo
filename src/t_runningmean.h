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

#ifndef __DEF_T_RUNNINGMEAN__
#define __DEF_T_RUNNINGMEAN__

#include "common.h"
#include "kahan.h"
#include "welford.h"
#include "runningmean.h"

#include <Rcpp.h>
using namespace Rcpp;

// running sums and means 

// running (weighted) sum or mean;
// and optimized by class of input?
//
// 2FIX: oneT and oneW should be the accumulator types, not one type.
// that is sum of logicals should go to int !?

template <typename RET,typename T,typename oneT,bool v_robustly,typename W,typename oneW,bool w_robustly,ReturnWhat retwhat,bool has_wts,bool do_recompute,bool na_rm>
RET t_runningSumish(T v,
                    Rcpp::Nullable< Rcpp::NumericVector > opt_time,
                    Rcpp::Nullable< Rcpp::NumericVector > opt_time_deltas,
                    double window,
                    W wts,
                    Rcpp::Nullable< Rcpp::NumericVector > opt_lb_time,
                    const int min_df, const int restart_period,
                    const bool variable_win, const bool wts_as_delta, const bool check_wts) {

    if (min_df < 0) { stop("BAD CODE: must give positive min_df"); }

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
                    time_deltas = as<NumericVector>(wts);
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

    oneT nextv, prevv;

    Kahan<oneT> fvsum;
    Kahan<oneW> fwsum;

    int nel;
    // subtraction count
    int subcount;

    oneW nextw, prevw;
    if (has_wts) {
        if (wts.size() < v.size()) { stop("size of wts does not match v"); }
    } else {
        nextw = oneW(1); 
        prevw=nextw;
    }

    fvsum = oneT(0);
    if (has_wts) { fwsum = oneW(0); }
    nel = 0;
    subcount = 0;

    //2FIX: use restart_period and do_recompute ... 
    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = NumericVector::is_na(window);
    if ((window <= 0) && (!infwin)) { stop("must give positive window"); }
    if (variable_win && !infwin) { Rcpp::warning("variable_win specified, but not being used as a non-na window is given."); }

    // whether to use the gap between lb_time as the effective window
    const bool gapwin = variable_win && infwin;

    double tf,t0;
    double prev_tf;
    const double tminf = time[0] - 1.0;  // effectively -inf?
    if (!gapwin && infwin) { t0 = tminf; }


    int iii,jjj,lll,tr_iii,tr_jjj;
    int numel = v.size();
    if (time.size() != numel) {
        stop("size of time does not match v"); // nocov
    }
    const int numlb = lb_time.size();

    RET xret(numlb);

    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
    }
    // 2FIX: this all has to be recouched in terms of a time window...
    
    tr_jjj = 0;
    tr_iii = -1;
    if (infwin) {
        prev_tf = tminf;
    } else {
        prev_tf = MIN(tminf,lb_time[0] - window - 1.0);  // make it less than t0 will be.
    }

    subcount = 0;
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
        if ((prev_tf <= t0) || (subcount >= restart_period)) {
            // could bisect, but lets not get fancy
            if (!infwin) { while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { tr_jjj++; } }
            tr_iii = tr_jjj;
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { tr_iii++; }

            // now sum up all elements from tr_jjj on the bottom, to (tr_iii-1) on the top, inclusive.

            // init//FOLDUP
            fvsum = oneT(0);
            if (has_wts) { fwsum = oneW(0); }
            nel = 0; //UNFOLD

            for (iii=tr_jjj;iii < tr_iii;++iii) {
                // add one
                if (has_wts) { //FOLDUP
                    nextw = wts[iii];
                } 
                nextv = v[iii];
                if (!na_rm) {
                    if (has_wts) {
                        fvsum += oneT(nextv * nextw);
                        fwsum += oneW(nextw);
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                } else if (!ISNAN(nextv)) {
                    if (has_wts) {
                        if (!(ISNAN(nextw) || (nextw <= 0))) {
                            fvsum += oneT(nextv * nextw);
                            fwsum += oneW(nextw);
                        }
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                }//UNFOLD
            }
            subcount = 0;
        } else {
            while ((tr_iii < numel) && (time[tr_iii] <= tf)) { 
                // add one //FOLDUP
                if (has_wts) { 
                    nextw = wts[tr_iii];
                } 
                nextv = v[tr_iii];
                if (! na_rm) {
                    if (has_wts) {
                        fvsum += oneT(nextv * nextw);
                        fwsum += oneW(nextw);
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                } else if (! ISNAN(nextv)) {
                    if (has_wts) {
                        if (! ((ISNAN(nextw) || (nextw <= 0)))) {
                            fvsum += oneT(nextv * nextw);
                            fwsum += oneW(nextw);
                        }
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                }//UNFOLD
                tr_iii++; 
            }
            if (!infwin) {
                while ((tr_jjj < numel) && (time[tr_jjj] <= t0)) { 
                    // remove one//FOLDUP
                    if (has_wts) { 
                        prevw = wts[tr_jjj];
                    } 
                    prevv = v[tr_jjj];
                    if (! na_rm) {
                        if (has_wts) {
                            fvsum -= oneT(prevv * prevw);
                            fwsum -= oneW(prevw);
                        } else {
                            fvsum -= oneT(prevv);
                            --nel;
                        }
                        if (do_recompute) { ++subcount; }
                    } else if (! ISNAN(prevv)) {
                        if (has_wts) {
                            if (! ((ISNAN(prevw) || (prevw <= 0)))) {
                                fvsum -= oneT(prevv * prevw);
                                fwsum -= oneW(prevw);
                            if (do_recompute) { ++subcount; }
                            }
                        } else {
                            fvsum -= oneT(prevv);
                            --nel;
                            if (do_recompute) { ++subcount; }
                        }
                    }//UNFOLD
                    tr_jjj++; 
                }
            }
            // may need to recompute based on the number of subtractions. bummer.
            if (subcount >= restart_period) {
                // now sum up all elements from tr_jjj on the bottom, to (tr_iii-1) on the top, inclusive.

                // init//FOLDUP
                fvsum = oneT(0);
                if (has_wts) { fwsum = oneW(0); }
                nel = 0; //UNFOLD

                for (iii=tr_jjj;iii < tr_iii;++iii) {
                    // add one
                    if (has_wts) { //FOLDUP
                        nextw = wts[iii];
                    } 
                    nextv = v[iii];
                    if (!na_rm) {
                        if (has_wts) {
                            fvsum += oneT(nextv * nextw);
                            fwsum += oneW(nextw);
                        } else {
                            fvsum += oneT(nextv);
                            ++nel;
                        }
                    } else if (!ISNAN(nextv)) {
                        if (has_wts) {
                            if (!(ISNAN(nextw) || (nextw <= 0))) {
                                fvsum += oneT(nextv * nextw);
                                fwsum += oneW(nextw);
                            }
                        } else {
                            fvsum += oneT(nextv);
                            ++nel;
                        }
                    }//UNFOLD
                }
                subcount = 0;
            }
        }
        // store em//FOLDUP
        if (has_wts) {
            if (fwsum.as() < min_df) {
                xret[lll] = oneT(NA_REAL); 
            } else { 
                if (retwhat==ret_sum) {
                    xret[lll] = fvsum.as(); 
                } else {
                    xret[lll] = fvsum.as() /  double(fwsum.as());
                }
            }
        } else {
            if (nel < min_df) {
                xret[lll] = oneT(NA_REAL); 
            } else { 
                if (retwhat==ret_sum) {
                    xret[lll] = fvsum.as(); 
                } else {
                    xret[lll] = fvsum.as()/nel; 
                }
            }
        }//UNFOLD
        prev_tf = tf;
    }//UNFOLD
    return xret;
}


template <typename T,typename oneT,bool v_robustly,typename W,typename oneW,bool w_robustly,ReturnWhat retwhat,bool has_wts,bool do_recompute>
SEXP t_runningSumishCurryOne(T v,
                             Rcpp::Nullable< Rcpp::NumericVector > time,
                             Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                             double window,
                             W wts,
                             Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                             const bool na_rm, const int min_df, const int restart_period,
                             const bool variable_win, const bool wts_as_delta, const bool check_wts,
                             const bool return_int) {
    if (return_int) {
        if (na_rm) {
            return wrap(t_runningSumish<IntegerVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,true>(v,time,time_deltas,window,wts,lb_time,min_df,restart_period,
                                                                                                                             variable_win,wts_as_delta,check_wts)); 
        } else {
            return wrap(t_runningSumish<IntegerVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,false>(v,time,time_deltas,window,wts,lb_time,min_df,restart_period,
                                                                                                                              variable_win,wts_as_delta,check_wts)); 
        }
    }
    if (na_rm) {
        return wrap(t_runningSumish<NumericVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,true>(v,time,time_deltas,window,wts,lb_time,min_df,restart_period,
                                                                                                                         variable_win,wts_as_delta,check_wts)); 
    }
    return wrap(t_runningSumish<NumericVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,false>(v,time,time_deltas,window,wts,lb_time,min_df,restart_period,
                                                                                                                      variable_win,wts_as_delta,check_wts)); 
}

template <typename T,typename oneT,bool v_robustly,ReturnWhat retwhat,bool do_recompute>
SEXP t_runningSumishCurryTwo(T v,
                             Rcpp::Nullable< Rcpp::NumericVector > time,
                             Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                             double window,
                             SEXP wts,
                             Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                             const bool na_rm, const int min_df, const int restart_period,
                             const bool variable_win, const bool wts_as_delta, const bool check_wts,
                             const bool return_int) {
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return t_runningSumishCurryOne<T,oneT,v_robustly,IntegerVector,int,false,retwhat,true,do_recompute>(v,time,time_deltas,window,wts,lb_time,
                                                                                                                                na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                                                return_int); }
            case REALSXP: { return t_runningSumishCurryOne<T,oneT,v_robustly,NumericVector,double,true,retwhat,true,do_recompute>(v,time,time_deltas,window,wts,lb_time,
                                                                                                                                  na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                                                  false); } // SIC: when double weights, cannot return int
            // to make smaller binaries, and because who cares about logicals, I convert them to integers here...
            case  LGLSXP: { return t_runningSumishCurryOne<T,oneT,v_robustly,IntegerVector,int,false,retwhat,true,do_recompute>(v,time,time_deltas,window,as<IntegerVector>(wts),lb_time,
                                                                                                                                 na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                                                 return_int); }
            default: stop("Unsupported weight type"); // nocov
        }
    }
    NumericVector dummy_wts;
    return t_runningSumishCurryOne<T,oneT,v_robustly,NumericVector,double,true,retwhat,false,do_recompute>(v,time,time_deltas,window,dummy_wts,lb_time,
                                                                                                           na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                           return_int); 
}

template <ReturnWhat retwhat,bool do_recompute>
SEXP t_runningSumishCurryThree(SEXP v,
                               Rcpp::Nullable< Rcpp::NumericVector > time,
                               Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                               double window,
                               SEXP wts,
                               Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                               const bool na_rm, const int min_df, const int restart_period,
                               const bool variable_win, const bool wts_as_delta, const bool check_wts,
                               const bool return_int) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return t_runningSumishCurryTwo<IntegerVector, int, false, retwhat, do_recompute>(v,time,time_deltas,window,wts,lb_time,
                                                                                                         na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                         return_int); }
        case REALSXP: { return t_runningSumishCurryTwo<NumericVector, double, true, retwhat, do_recompute>(v,time,time_deltas,window,wts,lb_time,
                                                                                                         na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                         return_int); }
        // to make smaller binaries, and because who cares about logicals, I convert them to integers here...
        case  LGLSXP: { return t_runningSumishCurryTwo<IntegerVector, int, false, retwhat, do_recompute>(as<IntegerVector>(v),time,time_deltas,window,wts,lb_time,
                                                                                                         na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                                                                         return_int); }
        default: stop("Unsupported input type"); // nocov
    }
    // CRAN checks are broken: 'warning: control reaches end of non-void function'
    // ... only for a crappy automated warning.
    return wrap(NumericMatrix(1,1));
}

template <ReturnWhat retwhat>
SEXP t_runningSumishCurryFour(SEXP v,
                              Rcpp::Nullable< Rcpp::NumericVector > time,
                              Rcpp::Nullable< Rcpp::NumericVector > time_deltas,
                              double window,
                              SEXP wts,
                              Rcpp::Nullable< Rcpp::NumericVector > lb_time,
                              const bool na_rm, const int min_df, const int restart_period,
                              const bool variable_win, const bool wts_as_delta, const bool check_wts) {

    const bool return_int=( ((TYPEOF(v)==INTSXP) || (TYPEOF(v)==LGLSXP)) &&
                            (retwhat==ret_sum) );
    const bool do_recompute = !IntegerVector::is_na(restart_period);
    if (do_recompute) {
        return t_runningSumishCurryThree<retwhat,true>(v,time,time_deltas,window,wts,lb_time,
                                                       na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                       return_int);
    }
    return t_runningSumishCurryThree<retwhat,false>(v,time,time_deltas,window,wts,lb_time,
                                                    na_rm,min_df,restart_period,variable_win,wts_as_delta,check_wts,
                                                    return_int);
}

#endif /* __DEF_T_RUNNINGMEAN__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
