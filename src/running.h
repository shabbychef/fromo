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

#ifndef __DEF_RUNNING__
#define __DEF_RUNNING__

#include "common.h"
#include "kahan.h"
#include "welford.h"

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

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool ord_beyond,bool renormalize,bool na_rm>
NumericMatrix runQM(T v,
                    W wts,
                    const int ord,
                    const int window,
                    const int recom_period,
                    const int lookahead,
                    const int min_df,
                    const double used_df,
                    const bool check_wts,
                    const bool normalize_wts) {

    Welford<oneW,has_wts,ord_beyond,na_rm> frets = Welford<oneW,has_wts,ord_beyond,na_rm>(ord);

    const int numel = v.size();

    double nextv, prevv, nextw, prevw;

    if (has_wts) {
        if (wts.size() < numel) { stop("size of wts does not match v"); }
    }

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(window);
    if ((window < 1) && (!infwin)) { stop("must give positive window"); }

    const int quasiwin = (infwin)? (numel):(window);

    if (min_df < 0) { stop("require positive min_df"); }
    if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

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
    bool aligned = (lookahead == 0);

    // super gross; I need these for the include later.
    double sg_denom,renorm,denom,sigmasq,sigma,sigmapow,mydf,dwsum,skew,exkurt,sr;
    int mmm;

    const int firstpart = MIN(numel,quasiwin);

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
        
    NumericMatrix xret(numel,ncols);

    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
    }
    // set them once only
    if (!has_wts) {
        nextw = 1.0;
        prevw = 1.0;
    }

    if (aligned) {
        // aligned case
        // sigh. broken.
        // as an invariant, we will start the computation
        // with frets, which is initialized as the summed
        // means on [jjj,lll]
        //

        // now run through lll index//FOLDUP
        for (lll=0;lll < firstpart;++lll) {
            // check subcount first and just recompute if needed.
            if (frets.subcount() >= recom_period) {
                // fix this
                frets = quasiWeightedThing<T,W,oneW,has_wts,ord_beyond,na_rm>(v,wts,ord,
                                                                              0,         //bottom
                                                                              lll+1,     //top
                                                                              false);    //no need to check weights as we have done it once above.
            } else {
                // add on nextv:
                nextv = double(v[lll]);
                if (has_wts) { nextw = double(wts[lll]); }  
                frets.add_one(nextv,nextw); 
            }

            // fill in the value in the output.
            // 2FIX: give access to v, not v[lll]...
            // moment_converter<retwhat, Welford<oneW,has_wts,ord_beyond,na_rm> ,T,renormalize>::mom_interp(xret,lll,ord,frets,v,used_df,min_df);
//yuck!!
#include "moment_interp.h"
        }//UNFOLD
        if (firstpart < numel) {
            tr_jjj = 0;
            // now run through lll index//FOLDUP
            for (lll=firstpart;lll < numel;++lll) {
                // check subcount first and just recompute if needed.
                if (frets.subcount() >= recom_period) {
                    // fix this
                    jjj = tr_jjj+1;
                    frets = quasiWeightedThing<T,W,oneW,has_wts,ord_beyond,na_rm>(v,wts,ord,
                                                                                  jjj,       //bottom
                                                                                  lll+1,     //top
                                                                                  false);    //no need to check weights as we have done it once above.
                } else {
                    // add on nextv:
                    nextv = double(v[lll]);
                    // remove prevv:
                    prevv = double(v[tr_jjj]);
                    if (has_wts) { 
                        nextw = double(wts[lll]); 
                        prevw = double(wts[tr_jjj]); 
                    }
                    frets.swap_one(nextv,nextw,prevv,prevw); 
                }
                tr_jjj++;

                // fill in the value in the output.
                // 2FIX: give access to v, not v[lll]...
                // moment_converter<retwhat, Welford<oneW,has_wts,ord_beyond,na_rm> ,T,renormalize>::mom_interp(xret,lll,ord,frets,v,used_df,min_df);
    //yuck!!
#include "moment_interp.h"
            }//UNFOLD
        }
    } else {
        // nonaligned case
        // as an invariant, we will start the computation
        // with frets, which is initialized as the summed
        // means on [jjj,iii]
        tr_iii = lookahead - 1;
        tr_jjj = lookahead - quasiwin;

        // now run through lll index//FOLDUP
        for (lll=0;lll < numel;++lll) {
            tr_iii++;
            // 2FIX: I feel like there is a corner case here waiting to bite me.
            // perhaps no computation is performed for lll=0, but
            // needs to be done for lll=1? 
            
            // check subcount first and just recompute if needed.
            if ((lll==0) || (frets.subcount() >= recom_period)) {
                // fix this
                iii = MIN(numel-1,tr_iii);
                jjj = MAX(0,tr_jjj+1);
                if (jjj <= iii) {
                    frets = quasiWeightedThing<T,W,oneW,has_wts,ord_beyond,na_rm>(v,wts,ord,
                                                                                  jjj,       //bottom
                                                                                  iii+1,     //top
                                                                                  false);    //no need to check weights as we have done it once above.
                }
            } else {
                if ((tr_iii < numel) && (tr_iii >= 0)) {
                    // add on nextv:
                    nextv = double(v[tr_iii]);
                    if (has_wts) { nextw = double(wts[tr_iii]); } 
                    frets.add_one(nextv,nextw); 
                }
                // remove prevv:
                if ((tr_jjj < numel) && (tr_jjj >= 0)) {
                    prevv = double(v[tr_jjj]);
                    if (has_wts) { prevw = double(wts[tr_jjj]); }
                    frets.rem_one(prevv,prevw); 
                }
            }
            tr_jjj++;

            // fill in the value in the output.
            // 2FIX: give access to v, not v[lll]...
            // moment_converter<retwhat, Welford<oneW,has_wts,ord_beyond,na_rm> ,T,renormalize>::mom_interp(xret,lll,ord,frets,v,used_df,min_df);
//yuck!!
#include "moment_interp.h"
        }//UNFOLD
    }
    return xret;
}

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool ord_beyond>
NumericMatrix runQMCurryZero(T v, 
                             W wts,
                             const int ord,
                             const int window,
                             const int recom_period,
                             const int lookahead,
                             const int min_df,
                             const double used_df,
                             const bool na_rm,
                             const bool check_wts,
                             const bool normalize_wts) {
    if (has_wts && normalize_wts) {
        if (na_rm) {
            return runQM<T,retwhat,W,oneW,has_wts,ord_beyond,true,true>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, check_wts, normalize_wts); 
        } else {
            return runQM<T,retwhat,W,oneW,has_wts,ord_beyond,true,false>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, check_wts, normalize_wts); 
        }
    } 
    if (na_rm) {
        return runQM<T,retwhat,W,oneW,has_wts,ord_beyond,false,true>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, check_wts, normalize_wts); 
    } 
    return runQM<T,retwhat,W,oneW,has_wts,ord_beyond,false,false>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, check_wts, normalize_wts); 
}

template <typename T,ReturnWhat retwhat,bool ord_beyond>
NumericMatrix runQMCurryOne(T v, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts,
                            const int ord,
                            const int window,
                            const int recom_period,
                            const int lookahead,
                            const int min_df,
                            const double used_df,
                            const bool na_rm,
                            const bool check_wts,
                            const bool normalize_wts) {

    //2FIX: typeof wts?
    if (wts.isNotNull()) {
        return runQMCurryZero<T,retwhat,NumericVector,double,true,ord_beyond>(v, wts.get(), ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
    }
    NumericVector dummy_wts;
    return runQMCurryZero<T,retwhat,NumericVector,double,false,ord_beyond>(v, dummy_wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
}



template <typename T,ReturnWhat retwhat>
NumericMatrix runQMCurryTwo(T v, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts,
                            const int ord,
                            const int window,
                            const int recom_period,
                            const int lookahead,
                            const int min_df,
                            const double used_df,
                            const bool na_rm,
                            const bool check_wts,
                            const bool normalize_wts) {

    if (ord==2) {
        return runQMCurryOne<T,retwhat,false>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
    }
    return runQMCurryOne<T,retwhat,true>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
}

template <ReturnWhat retwhat>
NumericMatrix runQMCurryThree(SEXP v, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts,
                              const int ord,
                              const int window,
                              const int recom_period,
                              const int lookahead,
                              const int min_df,
                              const double used_df,
                              const bool na_rm,
                              const bool check_wts,
                              const bool normalize_wts) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return runQMCurryTwo<IntegerVector,retwhat>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); } 
        case REALSXP: { return runQMCurryTwo<NumericVector,retwhat>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); } 
        case  LGLSXP: { return runQMCurryTwo<IntegerVector,retwhat>(as<IntegerVector>(v), wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); }  // bools can be upcast to save build size.
        default: stop("Unsupported weight type"); // nocov
    }
    // have to have fallthrough for CRAN check.
    NumericMatrix retv;
    return retv;
}

#endif /* __DEF_RUNNING__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
