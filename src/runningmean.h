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

#ifndef __DEF_RUNNINGMEAN__
#define __DEF_RUNNINGMEAN__

#include "common.h"
#include "kahan.h"
#include "welford.h"

#include <Rcpp.h>
using namespace Rcpp;

// running sums and means 

//template <T>
//Rcpp::Nullable<T> get_wts(SEXP wts) {
    //if (!Rf_isNull(window)) {
        //return Rcpp::Nullable<T>(wts);
    //}
    //return Rcpp::Nullable<T>();
//}

// running (weighted) sum or mean;
// and optimized by class of input?
//
// 2FIX: oneT and oneW should be the accumulator types, not one type.
// that is sum of logicals should go to int !?

template <typename RET,typename T,typename oneT,bool v_robustly,typename W,typename oneW,bool w_robustly,ReturnWhat retwhat,bool has_wts,bool do_recompute,bool na_rm>
RET runningSumish(T v,
                  W wts,
                  int window,
                  const int min_df,
                  int recom_period,
                  const bool check_wts) {
    if (min_df < 0) { stop("BAD CODE: must give positive min_df"); }

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

    //2FIX: use recom_period and do_recompute ... 
    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(window);
    if ((window < 1) && (!infwin)) { stop("must give positive window"); }

    int iii,jjj,lll;
    int numel = v.size();

    RET xret(numel);

    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
    }

    jjj = 0;
    // now run through iii
    for (iii=0;iii < numel;++iii) {
        if (!do_recompute || (subcount < recom_period)) {
            // add one
            if (has_wts) { //FOLDUP
                nextw = wts[iii];
            } 
            nextv = v[iii];
            if (! na_rm) {
                if (has_wts) {
                    fvsum += oneT(nextv * nextw);
                    fwsum += oneW(nextw);
                } else {
                    fvsum += oneT(nextv);
                    ++nel;
                }
            } else {
                if (! ISNAN(nextv)) {
                    if (has_wts) {
                        if (! ((ISNAN(nextw) || (nextw <= 0)))) {
                            fvsum += oneT(nextv * nextw);
                            fwsum += oneW(nextw);
                        }
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                }
            }//UNFOLD
            // remove one
            if (!infwin && (iii >= window)) {
                if (has_wts) { //FOLDUP
                    prevw = wts[jjj];
                } 
                prevv = v[jjj];
                if (! na_rm) {
                    if (has_wts) {
                        fvsum -= oneT(prevv * prevw);
                        fwsum -= oneW(prevw);
                    } else {
                        fvsum -= oneT(prevv);
                        --nel;
                    }
                    if (do_recompute) { ++subcount; }
                } else {
                    if (! ISNAN(prevv)) {
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
                    }
                }//UNFOLD
                ++jjj;
            }
        } else {
            // flat out recompute;//FOLDUP
            // this seems a little odd, but note we are positively adding here,
            // not subtracting, so increment the jjj first.
            ++jjj;
            // init//FOLDUP
            fvsum = oneT(0);
            if (has_wts) { fwsum = oneW(0); }
            nel = 0; //UNFOLD
            for (lll=jjj;lll <= iii;++lll) {
                // add one
                if (has_wts) { //FOLDUP
                    nextw = wts[lll];
                } 
                nextv = v[lll];
                if (!na_rm) {
                    if (has_wts) {
                        fvsum += oneT(nextv * nextw);
                        fwsum += oneW(nextw);
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                } else {
                    // in this case no need to increment nel, we know it will be window? 
                    if (!ISNAN(nextv)) {
                        if (has_wts) {
                            if (!(ISNAN(nextw) || (nextw <= 0))) {
                                fvsum += oneT(nextv * nextw);
                                fwsum += oneW(nextw);
                            }
                        } else {
                            fvsum += oneT(nextv);
                            ++nel;
                        }
                    }
                }//UNFOLD
            }
            subcount = 0;//UNFOLD
        }
        // store em
        if (has_wts) {
            if (fwsum.as() < min_df) {
                xret[iii] = oneT(NA_REAL); 
            } else { 
                if (retwhat==ret_sum) {
                    xret[iii] = fvsum.as(); 
                } else {
                    if (has_wts) {
                        xret[iii] = fvsum.as() / double(fwsum.as());
                    } else {
                        xret[iii] = fvsum.as() / double(nel); 
                    }
                }
            }
        } else {
            if (nel < min_df) {
                xret[iii] = oneT(NA_REAL); 
            } else { 
                if (retwhat==ret_sum) {
                    xret[iii] = fvsum.as(); 
                } else {
                    if (has_wts) {
                        xret[iii] = fvsum.as() / double(fwsum.as());
                    } else {
                        xret[iii] = fvsum.as() / double(nel); 
                    }
                }
            }
        }
    }
    return xret;
}


template <typename T,typename oneT,bool v_robustly,typename W,typename oneW,bool w_robustly,ReturnWhat retwhat,bool has_wts,bool do_recompute>
SEXP runningSumishCurryOne(T v,
                           W wts,
                           int window,
                           const int min_df,
                           int recom_period,
                           const bool na_rm,
                           const bool check_wts,
                           const bool return_int) {
   if (return_int) {
       if (na_rm) {
           return wrap(runningSumish<IntegerVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,true>(v,wts,window,min_df,recom_period,check_wts));
       } else {
           return wrap(runningSumish<IntegerVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,false>(v,wts,window,min_df,recom_period,check_wts));
       }
   }
   if (na_rm) {
       return wrap(runningSumish<NumericVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,true>(v,wts,window,min_df,recom_period,check_wts));
   }
   return wrap(runningSumish<NumericVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute,false>(v,wts,window,min_df,recom_period,check_wts));
}

template <typename T,typename oneT,bool v_robustly,ReturnWhat retwhat,bool do_recompute>
SEXP runningSumishCurryTwo(T v,
                           SEXP wts,
                           int window,
                           const int min_df,
                           int recom_period,
                           const bool na_rm,
                           const bool check_wts,
                           const bool return_int) {
    // 2FIX: to get smaller images, use IntegerVector instead of logicals and convert to 0/1
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return runningSumishCurryOne<T,oneT,v_robustly,IntegerVector,int,false,retwhat,true,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts,return_int); }
            case REALSXP: { return runningSumishCurryOne<T,oneT,v_robustly,NumericVector,double,true,retwhat,true,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts,false); } // SIC: when double weights, cannot return int
            // to make smaller binaries, and because who cares about logicals, I convert them to integers here...
            case  LGLSXP: { return runningSumishCurryOne<T,oneT,v_robustly,IntegerVector,int,false,retwhat,true,do_recompute>(v,as<IntegerVector>(wts),window,min_df,recom_period,na_rm,check_wts,return_int); }
            default: stop("Unsupported weight type"); // #nocov
        }
    }
    NumericVector dummy_wts;
    return runningSumishCurryOne<T,oneT,v_robustly,NumericVector,double,true,retwhat,false,do_recompute>(v,dummy_wts,window,min_df,recom_period,na_rm,check_wts,return_int);
}

template <ReturnWhat retwhat,bool do_recompute>
SEXP runningSumishCurryThree(SEXP v,
                             SEXP wts,
                             int window,
                             const int min_df,
                             int recom_period,
                             const bool na_rm,
                             const bool check_wts,
                             const bool return_int) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return runningSumishCurryTwo<IntegerVector, int, false, retwhat, do_recompute>(v, wts, window, min_df, recom_period, na_rm, check_wts, return_int); }
        case REALSXP: { return runningSumishCurryTwo<NumericVector, double, true, retwhat, do_recompute>(v, wts, window, min_df, recom_period, na_rm, check_wts, return_int); }
        // to make smaller binaries, and because who cares about logicals, I convert them to integers here...
        case  LGLSXP: { return runningSumishCurryTwo<IntegerVector, int, false, retwhat, do_recompute>(as<IntegerVector>(v), wts, window, min_df, recom_period, na_rm, check_wts, return_int); }
        default: stop("Unsupported input type"); // #nocov
    }
    // CRAN checks are broken: 'warning: control reaches end of non-void function'
    // ... only for a crappy automated warning.
    return wrap(NumericMatrix(1,1));
}

template <ReturnWhat retwhat>
SEXP runningSumishCurryFour(SEXP v,
                            SEXP wts,
                            int window,
                            const int min_df,
                            int recom_period,
                            const bool na_rm,
                            const bool check_wts) {
    const bool return_int=( ((TYPEOF(v)==INTSXP) || (TYPEOF(v)==LGLSXP)) &&
                              (retwhat==ret_sum) );
    const bool do_recompute = !IntegerVector::is_na(recom_period);
    if (do_recompute) {
        return runningSumishCurryThree<retwhat,true>(v,wts,window,min_df,recom_period,na_rm,check_wts,return_int);
    }
    return runningSumishCurryThree<retwhat,false>(v,wts,window,min_df,recom_period,na_rm,check_wts,return_int);
}

#endif /* __DEF_RUNNINGMEAN__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
