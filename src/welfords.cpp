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

  numerically stable computation of sample moments.
 
  see also: 
  * http://www.johndcook.com/blog/2008/09/26/comparing-three-methods-of-computing-standard-deviation/
  * http://www.johndcook.com/standard_deviation.html
  * J. Bennett, et. al., 'Numerically Stable, Single-Pass,
    Parallel Statistics Algorithms,' Proceedings of IEEE
    International Conference on Cluster Computing, 2009.
    http://www.janinebennett.org/index_files/ParallelStatisticsAlgorithms.pdf
    https://www.semanticscholar.org/paper/Numerically-stable-single-pass-parallel-statistics-Bennett-Grout/a83ed72a5ba86622d5eb6395299b46d51c901265
  * T. Terriberry, 'Computing Higher-Order Moments Online,' 
    http://people.xiph.org/~tterribe/notes/homs.html

  Created: 2016.03.25
  Copyright: Steven E. Pav, 2016
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

// general note: this is *not* just a header file. 
// it has code in it. 
// I am using CPP 'magic' to keep it DRY.

#ifndef __DEF_WELFORD__
#define __DEF_WELFORD__

#include <math.h>
#include <Rcpp.h>
#include "common.h"
#include "kahan.cpp"
using namespace Rcpp;

template<class W,bool has_wts,bool na_rm,bool check_wts,bool normalize,MaxOrder maxord>
class Welford {
    public:
        int m_ord;
        NumericVector m_xx;
    public:
        inline Welford(const int &ord);
        inline Welford(const int &ord, 
                       const NumericVector &xx);

    public:
        // getters 
        inline int nel() const;
        inline NumericVector as() const;
        inline NumericVector asvec() const;
        inline NumericVector vecpart() const;
    public:
        // reset to zero
        inline Welford& tare();
    public:
        // interpreters
        inline double var(const double used_df) const;
        inline double mean() const;
        inline double sd(const double used_df) const;
        inline double skew() const;
        inline double exkurt() const;

        inline double sharpe(const double used_df) const { return double(m_xx[1]) / sd(used_df); }
        inline double centered(const double xval) const { return (xval - m_xx[1]); }
        inline double scaled(const double xval,const double used_df) const { return (xval/sd(used_df)); }
        inline double zscored(const double xval,const double used_df) const { return ((xval - m_xx[1])/sd(used_df)); }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt);
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt);
        // add one, remove one
        inline Welford& swap_one (const double addxval, const W addwt, const double remxval, const W remwt);
        // join two Welford objects together
        inline Welford& join(const Welford& rhs);
        // remove one from another
        inline Welford& unjoin(const Welford& rhs);
};

#endif /* __DEF_WELFORD__ */

using namespace Rcpp;


// Welford-Terriberry object.
// holds the following:
// m_nel : the number of distinct observations going into this sum.
// note that weighted observations are treated as one here.
// m_wsum : a Kahan<W> object of the total weight. When unweighted
// observations are used, this is just a duplicate of m_nel.
// m_ord : the integer order of the summation.
// a vector xx of length m_ord;
// xx[0] is ignored. keep it zero.
// xx[1] is weighted mean of x.
// xx[2] is weighted sum of (x - xx[1])^2
// xx[3] is weighted sum of (x - xx[1])^3
// xx[4] is weighted sum of (x - xx[1])^4
// ...

// Welford Terriberry
// ord_beyond must be used for (ord > 2)
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
//
//

template<class W>
class Welford<W,
#ifdef HAS_WTS
      true,
#else
      false,
#endif
#ifdef NA_RM
      true,
#else
      false,
#endif
#ifdef CHECK_WT
      true,
#else
      false,
#endif
#ifdef NORMALIZE
      true,
#else
      false,
#endif
      MAX_ORDER > {
    public:
        int m_ord;
        NumericVector m_xx;
    private:
#ifdef HAS_WTS
        Kahan<W> m_wsum;
#endif
        int m_nel;
    public:
#ifdef HAS_WTS
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_wsum(Kahan<W>(0)), m_xx(NumericVector(ord+1)) {}
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {}
        inline Welford(const int &ord, 
                       const int &nel,
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {}
#else
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_xx(NumericVector(ord+1)) {}
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_xx(NumericVector(xx)) {}

#endif
    public:
        // reset to zero
        inline Welford& tare() {
            m_nel = 0; 
#ifdef HAS_WTS
            m_wsum = W(0);
#endif
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const double used_df) const {
#ifdef HAS_WTS
#ifdef NORMALIZE
            double renorm;
            renorm = double(m_nel) / double(m_wsum.as());
            return ((renorm * m_xx[2]) / (double(m_nel) - used_df));
#else
            return ((m_xx[2]) / (double(m_wsum.as()) - used_df));
#endif
#else
            return ((m_xx[2]) / (double(m_nel) - used_df));
#endif
        }
        inline double mean() const { return m_xx[1]; }
        inline double sd(const double used_df) const { return sqrt(var(used_df)); }
        inline double skew() const {
#ifdef HAS_WTS
            return (sqrt(double(m_wsum.as())) * m_xx[3] / pow(m_xx[2],1.5));
#else
            return (sqrt(double(m_nel)) * m_xx[3] / pow(m_xx[2],1.5));
#endif
        }
        inline double exkurt() const {
#ifdef HAS_WTS
            return ((double(m_wsum.as()) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
#else
            return ((double(m_nel) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
#endif
        }
        // getters 
        inline int nel() const { return m_nel; }
#ifdef HAS_WTS
        inline W wsum() const { return m_wsum.as(); }
#endif
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
#ifdef HAS_WTS
            resu[0] = double(wsum());
#else
            resu[0] = double(m_nel);
#endif
            return resu;
        }
        inline NumericVector vecpart() const { 
            return m_xx;
        }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt) {
            double xb_les_muA, pre_del_mu, muD_les_muA, wtD, wtA;
            double term_left, div_left, rem_right, div_right, inner_term;
            // xval = x_b
            // wt = w_b
            // xb_les_muA = x_b - mu_A
#if(defined(HAS_WTS) && defined(CHECK_WT))
            if (wt < 0) { stop("negative weight detected"); }
#endif

#ifdef NA_RM
#ifdef HAS_WTS
            if (! (ISNAN(xval) || ISNAN(wt))) {
#else
            if (! (ISNAN(xval))) {
#endif
#endif
#ifdef HAS_WTS
                m_nel++; 
                wtA = double(m_wsum.as());
                m_wsum += wt;
                wtD = double(m_wsum.as());
#else
                wtA = double(m_nel);
                m_nel++;
                wtD = double(m_nel);
#endif
                xb_les_muA = xval - m_xx[1];
#ifdef HAS_WTS
                pre_del_mu  = xb_les_muA * double(wt);
                muD_les_muA = pre_del_mu / wtD;
#else
                muD_les_muA = xb_les_muA / wtD;
#endif
                m_xx[1] += muD_les_muA;
                // the mean is computed. drop out if ord==1
#if MAX_ORDER != ORDER_ONE

#if MAX_ORDER == ORDER_TWO 
#ifdef HAS_WTS
                m_xx[2] += pre_del_mu * (xval - m_xx[1]);
#else
                m_xx[2] += xb_les_muA * (xval - m_xx[1]);
#endif
#else

                div_left = -muD_les_muA;
                term_left = pow(div_left,m_ord) * wtA;
#ifdef HAS_WTS
                div_right = -wtA / double(wt);
#else
                div_right = -wtA;
#endif
                rem_right = pow(div_right,m_ord - 1);

                for (int ppp=m_ord;ppp >= 2;ppp--) {
                    m_xx[ppp] += term_left * (1.0 - rem_right);
                    if (ppp > 2) {
                        term_left /= div_left;
                        rem_right /= div_right;
                        inner_term = div_left;
                        for (int qqq=1;qqq <= ppp-2;qqq++) {
                            m_xx[ppp] += bincoef[ppp][qqq] * inner_term * m_xx[ppp-qqq];
                            if (qqq < ppp - 2) { inner_term *= div_left; }
                        }
                    }
                }
#endif

#ifdef NA_RM
            }
#endif
#endif
            return *this;
        }

        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            double xc_les_muA, pre_del_mu, muD_les_muA, wtD, wtA;
            double term_left, div_left, rem_right, div_right, inner_term;
            // xval = x_c
            // wt = w_c
            // xc_les_muA = x_c - mu_A
#if(defined(HAS_WTS) && defined(CHECK_WT))
            if (wt < 0) { stop("negative weight detected"); }
#endif

#ifdef NA_RM
#ifdef HAS_WTS
            if (! (ISNAN(xval) || ISNAN(wt))) {
#else
            if (! (ISNAN(xval))) {
#endif
#endif
#ifdef HAS_WTS
                m_nel--; 
                wtA = double(m_wsum.as());
                m_wsum -= wt;
                wtD = double(m_wsum.as());
#else
                wtA = double(m_nel);
                m_nel--;
                wtD = double(m_nel);
#endif
                xc_les_muA = xval - m_xx[1];
#ifdef HAS_WTS
                pre_del_mu  = xc_les_muA * double(wt);
                muD_les_muA = - pre_del_mu / wtD;
#else
                muD_les_muA = - xc_les_muA / wtD;
#endif
                m_xx[1] += muD_les_muA;
                // the mean is computed. drop out if ord==1
#if MAX_ORDER != ORDER_ONE

#if MAX_ORDER == ORDER_TWO 
#ifdef HAS_WTS
                m_xx[2] -= pre_del_mu * (xval - m_xx[1]);
#else
                m_xx[2] -= xc_les_muA * (xval - m_xx[1]);
#endif
#else

                div_left = -muD_les_muA;
                term_left = pow(div_left,m_ord) * wtA;
#ifdef HAS_WTS
                div_right = wtA / double(wt);
#else
                div_right = wtA;
#endif
                rem_right = pow(div_right,m_ord - 1);

                for (int ppp=m_ord;ppp >= 2;ppp--) {
                    m_xx[ppp] += term_left * (1.0 - rem_right);
                    if (ppp > 2) {
                        term_left /= div_left;
                        rem_right /= div_right;
                        inner_term = div_left;
                        for (int qqq=1;qqq <= ppp-2;qqq++) {
                            m_xx[ppp] += bincoef[ppp][qqq] * inner_term * m_xx[ppp-qqq];
                            if (qqq < ppp - 2) { inner_term *= div_left; }
                        }
                    }
                }
#endif

#ifdef NA_RM
            }
#endif
#endif
            return *this;
        }


        // remove one (weighted) and add one (weighted) observation from our set of x
        inline Welford& swap_one (const double addxval, const W addwt,
                                  const double remxval, const W remwt) {
            // you only get speedups when the following occur
            // no weights (or rather equal weights, which we do not want to check), and
            // order == 2, and
            // neither is NA.
#if(DEFINED(HAS_WTS) || (MAX_ORDER != ORDER_TWO))
                add_one(addxval,addwt);
                rem_one(remxval,remwt);
#else

#ifdef NA_RM
            if (! (ISNAN(addxval))) {
                rem_one(remxval,remwt);
            } else if (! (ISNAN(remxval))) {
                add_one(addxval,addwt);
            } else {
#endif
                double xb_les_xc, wtD, muA;
                xb_les_xc = addxval - remxval;
                wtD = double(m_nel);
                muA = m_xx[1];
                m_xx[1] += xb_les_xc / wtD;
                m_xx[2] += xb_les_xc * (addxval + remxval - (m_xx[1] + muA));
#ifdef NA_RM
            }
#endif

#endif
            return *this;
        }

        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double wtA, wtB, wtD,
                   muB_les_muA, muD_les_muA, muB_les_muD,
                   aterm_div, aterm_first, bterm_div, bterm_first,
                   cterm, dterm;
            int ppp,qqq;

#ifdef HAS_WTS
            wtA = double(m_wsum.as());
#else
            wtA = double(m_nel);
#endif
            if (wtA <= 0) {
                m_nel = rhs.m_nel;
                // clone?
#ifdef HAS_WTS
                m_wsum = rhs.m_wsum;
#endif
                // clone!
                m_xx = Rcpp::clone(rhs.m_xx);
                return *this;
            }
            // else onboard the observations
            m_nel += rhs.m_nel;
#ifdef HAS_WTS
            wtB = double(rhs.m_wsum.as());
            m_wsum += rhs.m_wsum;
            wtD = double(m_wsum.as());
#else
            wtB = double(rhs.m_nel);
            wtD = double(m_nel);
#endif

            if (wtB <= 0) {
                return *this;
            }

            muB_les_muA = rhs.m_xx[1] - m_xx[1];
            muD_les_muA = (wtB / wtD) * muB_les_muA;
            m_xx[1] += muD_les_muA;
#if MAX_ORDER != ORDER_ONE
            muB_les_muD = rhs.m_xx[1] - m_xx[1];

            aterm_div = -muD_les_muA;
            aterm_first = wtA * pow(aterm_div,m_ord);
            bterm_div = muB_les_muD;
            bterm_first = wtB * pow(bterm_div,m_ord);

            for (ppp=m_ord;ppp >= 2;ppp--) {
                m_xx[ppp] += rhs.m_xx[ppp] + aterm_first + bterm_first;
#if MAX_ORDER == ORDER_BEYOND
                if (ppp > 2) {
                    aterm_first /= aterm_div;
                    bterm_first /= bterm_div;
                    cterm = aterm_div;
                    dterm = bterm_div;
                    for (int qqq=1;qqq <= (ppp-2); qqq++) {
                        m_xx[ppp] += bincoef[ppp][qqq] * (cterm * m_xx[ppp-qqq] + dterm * rhs.m_xx[ppp-qqq]);
                        if (qqq < (ppp-2)) {
                            cterm *= aterm_div;
                            dterm *= bterm_div;
                        }
                    }
                }
#endif
            }
#endif
            return *this;
        }

        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double wtA, wtC, wtD,
                   muC_les_muA, muD_les_muA, muC_les_muD,
                   aterm_div, aterm_first, bterm_div, bterm_first,
                   cterm, dterm;
            int ppp,qqq;

#ifdef HAS_WTS
            wtC = double(rhs.m_wsum.as());
#else
            wtC = double(rhs.m_nel);
#endif
            if (wtC <= 0) {
                return *this;
            }


#ifdef HAS_WTS
            wtA = double(m_wsum.as());
#else
            wtA = double(m_nel);
#endif
            // else onboard the observations
            m_nel -= rhs.m_nel;
#ifdef HAS_WTS
            m_wsum -= rhs.m_wsum;
            wtD = double(m_wsum.as());
#else
            wtD = double(m_nel);
#endif
            if (wtD <= 0) { stop("cannot subtract more observations than were seen."); }

            muC_les_muA = rhs.m_xx[1] - m_xx[1];
            muD_les_muA = (- wtC / wtD) * muC_les_muA;
            m_xx[1] += muD_les_muA;
#if MAX_ORDER != ORDER_ONE
            muC_les_muD = rhs.m_xx[1] - m_xx[1];

            aterm_div = -muD_les_muA;
            aterm_first = wtA * pow(aterm_div,m_ord);
            bterm_div = muC_les_muD;
            bterm_first = - wtC * pow(bterm_div,m_ord);

            for (ppp=m_ord;ppp >= 2;ppp--) {
                m_xx[ppp] += aterm_first + bterm_first - rhs.m_xx[ppp];
#if MAX_ORDER == ORDER_BEYOND
                if (ppp > 2) {
                    aterm_first /= aterm_div;
                    bterm_first /= bterm_div;
                    cterm = aterm_div;
                    dterm = bterm_div;
                    for (int qqq=1;qqq <= (ppp-2); qqq++) {
                        m_xx[ppp] += bincoef[ppp][qqq] * (cterm * m_xx[ppp-qqq] - dterm * rhs.m_xx[ppp-qqq]);
                        if (qqq < (ppp-2)) {
                            cterm *= aterm_div;
                            dterm *= bterm_div;
                        }
                    }
                }
#endif
            }
#endif
            return *this;
        }

};

// initialize.

template<class T,class Wvec,class W>
Welford<W, 
#ifdef HAS_WTS
      true,
#else
      false,
#endif
#ifdef NA_RM
      true,
#else
      false,
#endif
#ifdef CHECK_WT
      true,
#else
      false,
#endif
#ifdef NORMALIZE
      true,
#else
      false,
#endif
      MAX_ORDER > init(T x,Wvec wt,int ord) {
         // first compute the mean robustly;
         int top;
         int nel;
         int iii, jjj;
         double nextv,nextdiff,powdiff;
         double mu;
         Kahan<double> wtsum;
         NumericVector xx = NumericVector(ord+1);
         NumericVector errs = NumericVector(ord+1);
         double tmp1, tmp2;
#ifdef HAS_WTS
         Kahan<W> totwt;
#endif


#ifdef HAS_WTS
         W nextwt;
#endif
         top = x.size();
#ifdef NA_RM
         LogicalVector isok=LogicalVector(top);
#endif
         nel = 0;
         for (int iii=0;iii<top;++iii) {
             nextv = double(x[iii]);
#ifdef HAS_WTS
             nextw = W(wt[iii]);
#endif
#ifdef NA_RM
#ifdef HAS_WTS
             isok[iii] = !ISNAN(nextv) && !ISNAN(nextw) && !(nextw < 0);
#else
             isok[iii] = !ISNAN(nextv);
#endif
             if (isok[iii]) {
#endif
#ifdef HAS_WTS
                 wtsum += nextw * nextv;
                 totwt += nextw;
#else
                 wtsum += nextv;
#endif
                 nel++;
#ifdef NA_RM
             }
#endif
         }

#if MAX_ORDER != ORDER_ONE
         if (nel > 0) {
#ifdef HAS_WTS
             mu = double(wtsum.as()) / double(totwt.as());
#else
             mu = double(wtsum.as()) / double(nel);
#endif

             for (int iii=0;iii<top;++iii) {
#ifdef NA_RM
                 if (isok[iii]) {
#endif
                     nextdiff = double(x[iii]) - mu;
#ifdef HAS_WTS
                     powdiff = double(wt[iii]) * nextdiff * nextdiff;
#else
                     powdiff = nextdiff * nextdiff;
#endif
                     for (int jjj=2;jjj < ord;++jjj) {
                         KAHAN_ADD(xx[jjj],errs[jjj],powdiff,tmp1,tmp2);
                         powdiff *= nextdiff;
                     }
                     KAHAN_ADD(xx[jjj],errs[jjj],powdiff,tmp1,tmp2);

#ifdef NA_RM
                 }
#endif
             }
         }
#endif

         xx[0] = double(nel);
         xx[1] = mu;
         
         Welford<W,
#ifdef HAS_WTS
      true,
#else
      false,
#endif
#ifdef NA_RM
      true,
#else
      false,
#endif
#ifdef CHECK_WT
      true,
#else
      false,
#endif
#ifdef NORMALIZE
      true,
#else
      false,
#endif
      MAX_ORDER > retv = 
         Welford<W,
#ifdef HAS_WTS
      true,
#else
      false,
#endif
#ifdef NA_RM
      true,
#else
      false,
#endif
#ifdef CHECK_WT
      true,
#else
      false,
#endif
#ifdef NORMALIZE
      true,
#else
      false,
#endif
      MAX_ORDER >(ord,nel,xx);
         return retv;
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
