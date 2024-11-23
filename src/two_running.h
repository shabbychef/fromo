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

  Created: 2024.11.19
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_TWO_RUNNING__
#define __DEF_TWO_RUNNING__

#include "common.h"
#include "two_welford.h"

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

// start here. will have to add new RetWhat thingys for correlation, regression and so on.
// 2FIX: do we need the renormalize thingy?

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool na_rm>
NumericMatrix two_runQM(T v,
                        T vv,
                        W wts,
                        const int window,
                        const int recom_period,
                        const int min_df,
                        const double used_df,   // remove this?
                        const bool check_wts,
                        const bool renormalize,  // confusing to have had two versions of this.
                        const bool check_negative_moments) {

    // a bit of a hack here, but you must have ord >= 2 for Welford
    // objects, otherwise it hits a memory leak. I know that previously
    // I had hard coded Welford-like objects for the ord=1 case,
    // but that makes huge libraries. instead, just hack this MAX in here
    //const int fake_ord = MAX(ord,2);

    TwoWelford<oneW,has_wts,na_rm> frets = TwoWelford<oneW,has_wts,na_rm>();
    frets.tare();

    const int numel = v.size();
    if (v.size() != numel) { stop("size of v and vv do not match"); }

    double nextv, prevv, nextvv, prevvv, nextw, prevw;

    if (has_wts) {
        if (wts.size() < numel) { stop("size of wts does not match v"); }
    }
    const int ord = 2;

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(window);
    if ((window < 1) && (!infwin)) { stop("must give positive window"); }

    const int quasiwin = (infwin)? (numel):(window);

    if (min_df < 0) { stop("require positive min_df"); }
    if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

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

    // super gross; I need these for the include later.
    double denom;

    const int firstpart = MIN(numel,quasiwin);

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
        
    NumericMatrix xret(numel,ncols);

    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
    }
    // set them once only
    if (!has_wts) {
        nextw = 1.0;
        prevw = 1.0;
    }
    // aligned case
    // as an invariant, we will start the computation
    // with frets, which is initialized as the summed
    // means on [jjj,lll]

    // now run through lll index//FOLDUP
    for (lll=0;lll < firstpart;++lll) {
        // check subcount first and just recompute if needed.
        if (frets.subcount() >= recom_period) {
            //zero it out
            frets.tare();
            add_many<T,W,oneW,has_wts,na_rm>(frets,
                                             v,vv,
                                             wts,
                                             0,      //bottom
                                             lll+1,  //top
                                             false); //no need to check weights as we have done it once above.
        } else {
            // add on nextv:
            nextv = double(v[lll]);
            nextvv = double(vv[lll]);
            if (has_wts) { nextw = double(wts[lll]); }  
            frets.add_one(nextv,nextvv,nextw); 
            if (check_negative_moments && frets.has_heywood()) {
                //zero it out
                frets.tare();
                add_many<T,W,oneW,has_wts,na_rm>(frets,
                                                 v,vv,wts,
                                                 0,      //bottom
                                                 lll+1,  //top
                                                 false); //no need to check weights as we have done it once above.
            }
        }

        // fill in the value in the output.
        // 2FIX: give access to v, not v[lll]...
//yuck!!
#include "two_moment_interp.h"
    }//UNFOLD
    if (firstpart < numel) {
        tr_jjj = 0;
        // now run through lll index//FOLDUP
        for (lll=firstpart;lll < numel;++lll) {
            // check subcount first and just recompute if needed.
            if (frets.subcount() >= recom_period) {
                // fix this
                jjj = tr_jjj+1;
                //zero it out
                frets.tare();
                add_many<T,W,oneW,has_wts,na_rm>(frets,
                                                 v,vv,wts,
                                                 jjj,    //bottom
                                                 lll+1,  //top
                                                 false); //no need to check weights as we have done it once above.
            } else {
                // add on nextv:
                nextv = double(v[lll]);
                nextvv = double(vv[lll]);
                // remove prevv:
                prevv = double(v[tr_jjj]);
                prevvv = double(vv[tr_jjj]);
                if (has_wts) { 
                    nextw = double(wts[lll]); 
                    prevw = double(wts[tr_jjj]); 
                }
                frets.swap_one(nextv,nextvv,nextw,
                               prevv,prevvv,prevw); 
                if (check_negative_moments && frets.has_heywood()) {
                    // fix this
                    jjj = tr_jjj+1;
                    //zero it out
                    frets.tare();
                    add_many<T,W,oneW,has_wts,na_rm>(frets,
                                                     v,vv,wts,
                                                     jjj,    //bottom
                                                     lll+1,  //top
                                                     false); //no need to check weights as we have done it once above.
                }
            }
            tr_jjj++;

            // fill in the value in the output.
            // 2FIX: give access to v, not v[lll]...
//yuck!!
#include "two_moment_interp.h"
        }//UNFOLD
    }
    return xret;
}

// 2FIX: do we need the renormalize thingy?
template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts>
NumericMatrix two_runQMCurryZero(T v, T vv, 
                             W wts,
                             const int window,
                             const int recom_period,
                             const int min_df,
                             const double used_df,
                             const bool na_rm,
                             const bool check_wts,
                             const bool normalize_wts,
                             const bool check_negative_moments) {
    if (na_rm) {
        return two_runQM<T,retwhat,W,oneW,has_wts,true>(v, vv, wts, window, recom_period, 
                                                        min_df, used_df, check_wts, normalize_wts, check_negative_moments); 
    } 
    return two_runQM<T,retwhat,W,oneW,has_wts,false>(v, vv, wts,  window, recom_period, 
                                                     min_df, used_df, check_wts, normalize_wts, check_negative_moments); 
}

template <typename T,ReturnWhat retwhat>
NumericMatrix two_runQMCurryOne(T v, T vv,
                                Rcpp::Nullable< Rcpp::NumericVector > wts,
                                const int window,
                                const int recom_period,
                                const int min_df,
                                const double used_df,
                                const bool na_rm,
                                const bool check_wts,
                                const bool normalize_wts,
                                const bool check_negative_moments) {

    //2FIX: typeof wts?
    if (wts.isNotNull()) {
        return two_runQMCurryZero<T,retwhat,NumericVector,double,true>(v, vv, wts.get(), window, recom_period, 
                                                                   min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
    }
    NumericVector dummy_wts;
    return two_runQMCurryZero<T,retwhat,NumericVector,double,false>(v, vv, dummy_wts,  window, recom_period,
                                                                min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
}


template <ReturnWhat retwhat>
NumericMatrix two_runQMCurryTwo(SEXP v, SEXP vv,
                                Rcpp::Nullable< Rcpp::NumericVector > wts,
                                const int window,
                                const int recom_period,
                                const int min_df,
                                const double used_df,
                                const bool na_rm,
                                const bool check_wts,
                                const bool normalize_wts,
                                const bool check_negative_moments) {
    switch (TYPEOF(v)) {
        // 2FIX: we might have to cast vv to the same type as v?
        case  INTSXP: 
            { 
                switch (TYPEOF(vv)) {
                    case  INTSXP: 
                        { 
                            return two_runQMCurryOne<IntegerVector,retwhat>(v, vv, wts, window, recom_period, 
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        } 
                    case REALSXP: // cast int to numeric
                        {  
                            return two_runQMCurryOne<NumericVector,retwhat>(as<NumericVector>(v), vv, wts, window, recom_period, 
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to int
                        {
                            return two_runQMCurryOne<IntegerVector,retwhat>(v, as<IntegerVector>(vv), wts, window, recom_period, 
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        }
                    default: stop("Unsupported data type for vv"); // #nocov
                }
            }
        case REALSXP:
            {
                switch (TYPEOF(vv)) {
                    case  INTSXP: // cast int to numeric
                        { 
                            return two_runQMCurryOne<NumericVector, retwhat>(v, as<NumericVector>(vv), wts, window, recom_period, 
                                                                         min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        } 
                    case REALSXP: 
                        {  
                            return two_runQMCurryOne<NumericVector,retwhat>(v, vv, wts, window, recom_period, 
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to numeric
                        {
                            return two_runQMCurryOne<NumericVector,retwhat>(v, as<NumericVector>(vv), wts, window, recom_period,
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        }
                    default: stop("Unsupported data type for vv"); // #nocov
                }
            }
        case  LGLSXP: 
            {
                switch (TYPEOF(vv)) {
                    case  INTSXP: // cast logical to int 
                        { 
                            return two_runQMCurryOne<IntegerVector, retwhat>(as<IntegerVector>(v), vv, wts, window, recom_period,
                                                                         min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        } 
                    case REALSXP: // cast logical to numeric
                        {  
                            return two_runQMCurryOne<NumericVector,retwhat>(as<NumericVector>(v), vv, wts, window, recom_period,
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
                        }
                    case  LGLSXP: // cast logical to integer
                        {
                            return two_runQMCurryOne<IntegerVector,retwhat>(as<IntegerVector>(v), as<IntegerVector>(vv), wts, window, recom_period,
                                                                        min_df, used_df, na_rm, check_wts, normalize_wts, check_negative_moments); 
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

#endif /* __DEF_TWO_RUNNING__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
