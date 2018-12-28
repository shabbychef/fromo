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

#ifndef __DEF_FROMO__
#define __DEF_FROMO__

#include <math.h>
#include "common.h"

#endif /* __DEF_FROMO__ */

#include <Rcpp.h>
using namespace Rcpp;

// yeah yeah, I know.
#include "kahan.cpp"


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

template <typename RET,typename T,typename oneT,bool v_robustly,typename W,typename oneW,bool w_robustly,ReturnWhat retwhat,bool has_wts,bool do_recompute>
RET runningSumish(T v,
                  W wts,
                  int window,
                  const int min_df,
                  int recom_period,
                  const bool na_rm,
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

    if (has_wts && check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }

    jjj = 0;
    // now run through iii
    for (iii=0;iii < numel;++iii) {
        if (!do_recompute || (subcount < recom_period)) {
            // add one
            if (has_wts) { //FOLDUP
                nextw = wts[iii];
            } 
            nextv = v[iii];
            if (! (na_rm && (ISNAN(nextv) || (has_wts && (ISNAN(nextw) || (nextw <= 0)))))) { 
                if (has_wts) {
                    fvsum += oneT(nextv * nextw);
                    fwsum += oneW(nextw);
                } else {
                    fvsum += oneT(nextv);
                    ++nel;
                }
            }//UNFOLD
            // remove one
            if (!infwin && (iii >= window)) {
                if (has_wts) { //FOLDUP
                    prevw = wts[jjj];
                } 
                prevv = v[jjj];
                if (! (na_rm && (ISNAN(prevv) || (has_wts && (ISNAN(prevw) || (prevw <= 0)))))) { 
                    if (do_recompute) { ++subcount; }
                    if (has_wts) {
                        fvsum -= oneT(prevv * prevw);
                        fwsum -= oneW(prevw);
                    } else {
                        fvsum -= oneT(prevv);
                        --nel;
                    }
                }//UNFOLD
                ++jjj;
            }
        } else {
            // flat out recompute;//FOLDUP
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
                if (! (na_rm && (ISNAN(nextv) || (has_wts && (ISNAN(nextw) || (nextw <= 0)))))) { 
                    if (has_wts) {
                        fvsum += oneT(nextv * nextw);
                        fwsum += oneW(nextw);
                    } else {
                        fvsum += oneT(nextv);
                        ++nel;
                    }
                }//UNFOLD
            }
            subcount = 0;//UNFOLD
        }
        // store em
        if ((has_wts && (fwsum.as() < min_df)) || (!has_wts && (nel < min_df))) { 
            xret[iii] = oneT(NA_REAL); 
        } else { 
            if (retwhat==ret_sum) {
                xret[iii] = fvsum.as(); 
            } else {
                if (has_wts) {
                    xret[iii] = fvsum.as() /  double(fwsum.as());
                } else {
                    xret[iii] = fvsum.as()/nel; 
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
       return wrap(runningSumish<IntegerVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts));
   }
   return wrap(runningSumish<NumericVector,T,oneT,v_robustly,W,oneW,w_robustly,retwhat,has_wts,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts));
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
    NumericVector dummy_wts;
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return runningSumishCurryOne<T,oneT,v_robustly,IntegerVector,int,false,retwhat,true,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts,return_int); }
            case REALSXP: { return runningSumishCurryOne<T,oneT,v_robustly,NumericVector,double,true,retwhat,true,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts,false); } // SIC: when double weights, cannot return int
            case  LGLSXP: { return runningSumishCurryOne<T,oneT,v_robustly,LogicalVector,bool,false,retwhat,true,do_recompute>(v,wts,window,min_df,recom_period,na_rm,check_wts,return_int); }
            default: stop("Unsupported weight type"); // nocov
        }
    }
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
        case  LGLSXP: { return runningSumishCurryTwo<LogicalVector, bool, false, retwhat, do_recompute>(v, wts, window, min_df, recom_period, na_rm, check_wts, return_int); }
        default: stop("Unsupported input type");
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
