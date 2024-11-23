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
 

  Created: 2024.11.13
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_TWO_WELFORD__
#define __DEF_TWO_WELFORD__

#include "common.h"
#include "kahan.h"

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// Welford object on a pair of observations, storing the first and second (co)sums
//
// holds the following:
//
// m_nel : the (integer) number of distinct observations going into this sum.
// we only track this in the case of weighted observations, as in the
// unweighted case it is redundant to m_wsum.
//
// m_subc : the (integer) number of subtractions performed on this object.
// not accurate for (un)joins.
//
// m_wsum : a Kahan<W> object of the total weight. When unweighted
// observations are used, this is just a duplicate of m_nel.
//
// m_ord : the integer order of the summation, can only equal 2.
//
// a vector m_xx of length m_ord;
// m_xx[0] is ignored. keep it zero.
// m_xx[1] is weighted mean of x.
// m_xx[2] is weighted mean of y.
// m_xx[3] is weighted sum of (x - m_xx[1])^2
// m_xx[4] is weighted sum of (x - m_xx[1]) * (y - m_xx[2])
// m_xx[5] is weighted sum of (y - m_xx[2])^2

// TwoWelford 
// generic//FOLDUP
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
template<class W,bool has_wts,bool na_rm>
class TwoWelford {
    private: 
        int m_nel;
        int m_subc;
    private:
        Kahan<W> m_wsum;
    public:
        NumericVector m_xx;
    public:
        inline TwoWelford() : m_nel(0), m_subc(0), m_wsum(Kahan<W>(0)), m_xx(NumericVector(6)) {
        }
        inline TwoWelford(const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_nel(nel), m_subc(0), m_wsum(Kahan<W>(sumwt)), m_xx(NumericVector(xx)) {
            if (m_xx.size() != 6) { stop("wrong sized object, must give 6"); } // #nocov
        }
        inline TwoWelford(const NumericVector &xx) :  m_nel(int(xx[0])), m_subc(0), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {
            if (m_xx.size() != 6) { stop("wrong sized object, must give 6"); } // #nocov
        }
    public:
        // reset to zero
        inline TwoWelford& tare() {
            m_nel = 0;
            m_subc = 0;
            m_wsum = W(0);
            for (int iii=0;iii < 6;++iii) { m_xx[iii] = 0; }
            return *this;
        }
        inline double var_denominator(const bool normalize,const double used_df) const {
            double renorm;
            if (has_wts) {
                if (normalize) {
                    if (used_df == 0.0) {
                        return (double(m_wsum.as()));
                    } else {
                        renorm = double(m_nel) / double(m_wsum.as());
                        return (double(m_nel) - used_df) / renorm;
                    }
                } else {
                    return (double(m_wsum.as()) - used_df);
                }
            } else {
                return (double(m_nel) - used_df);
            }
        }
        inline double x_var(const bool normalize,const double used_df) const {
            return m_xx[3] / var_denominator(normalize, used_df);
        }
        inline double y_var(const bool normalize,const double used_df) const {
            return m_xx[5] / var_denominator(normalize, used_df);
        }
        inline double x_mean() const {
            return m_xx[1];
        }
        inline double y_mean() const {
            return m_xx[2];
        }
        inline double x_sd(const bool normalize,const double used_df) const {
            return sqrt(x_var(normalize,used_df));
        }
        inline double y_sd(const bool normalize,const double used_df) const {
            return sqrt(y_var(normalize,used_df));
        }
        inline double correlation() const {
            return m_xx[4] / sqrt(m_xx[3] * m_xx[5]);
        }
        inline double covariance(const bool normalize,const double used_df) const {
            return m_xx[4] / var_denominator(normalize, used_df);
        }
        inline double regression_slope() const {
            // assumes that v = x and vv = y and we want the slope for y = mx + b
            return m_xx[4] / m_xx[3];
        }
        inline double regression_intercept() const {
            // assumes that v = x and vv = y and we want the intercept for y = mx + b
            return m_xx[2] - m_xx[1] * m_xx[4] / m_xx[3];
        }
        // I really hate C++. here I will assign values into a matrix at a given position because
        // returning two values is just too difficult. 
        // assigns the intercept and slope into the row of xret.
        inline void assign_regression_fit(NumericMatrix xret, const int rownum) const {
            const double slope = m_xx[4] / m_xx[3];
            xret(rownum, 1) = slope;
            xret(rownum, 0) = m_xx[2] - m_xx[1] * slope;
        }
        // assigns the intercept and slope, standard error, s.e. of the intercept and slope into the row of xret.
        inline void assign_regression_diagnostics(NumericMatrix xret, const int rownum, const bool normalize,const double used_df) const {
            // assign_regression_diagnostics(xret, rownum, normalize, used_df);
            const double slope = m_xx[4] / m_xx[3];
            // slope
            xret(rownum, 1) = slope;
            // intercept
            xret(rownum, 0) = m_xx[2] - m_xx[1] * slope;
            const double denom = var_denominator(normalize, used_df);
            const double reg_se = sqrt((m_xx[5] - m_xx[4] * slope)/denom);
            const double slope_se = reg_se / sqrt(m_xx[3]);
            xret(rownum, 2) = reg_se;
            // slope se
            xret(rownum, 4) = slope_se;
            double na;
            if (has_wts) {
                na = double(m_wsum.as());
            } else {
                na = double(m_nel);
            }
            // intercept se
            xret(rownum, 3) = slope_se * sqrt(m_xx[3]/na + m_xx[1]*m_xx[1]);
        }
        // getters 
        inline int nel() const { if (has_wts) { return m_nel; } else { return int(wsum()); } }  // not sure I understand this...
        inline int subcount() const { return m_subc; }

        // return true if any even order sums are negative, or if the matrix has a negative eigenvalue.
        inline bool has_heywood() const {
            return ((m_xx[3] < 0) || (m_xx[5] < 0) || (m_xx[3] * m_xx[5] < m_xx[4] * m_xx[4]));
        }
        inline W wsum() const { 
            if (has_wts) { return m_wsum.as(); }
            return W(m_nel);
        }
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
            resu[0] = double(wsum());
            return resu;
        }
        inline NumericVector vecpart() const { 
            return m_xx;
        }
    public:
        // add another (weighted) observation to our set of x
        inline TwoWelford& add_one (const double xval, const double yval, const W wt) {
            if (na_rm) {
                if (ISNAN(xval) || ISNAN(yval)) { return *this; }
                if (has_wts) {
                    if (ISNAN(wt) || (wt <= 0)) {
                        return *this;
                    }
                }
            }

            double xb_les_muA, x_pre_del_mu, x_muD_les_muA, 
                   yb_les_muA, y_pre_del_mu, y_muD_les_muA, 
                   y_post_diff,
                   wtD, wtA;

            if (has_wts) {
                m_nel++; 
                wtA = double(m_wsum.as());
                m_wsum += wt;
                wtD = double(m_wsum.as());
            } else {
                wtA = double(m_nel);
                m_nel++;
                wtD = double(m_nel);
            }
            xb_les_muA = xval - m_xx[1];
            yb_les_muA = yval - m_xx[2];
            if (has_wts) {
                x_pre_del_mu  = xb_les_muA * double(wt);
                x_muD_les_muA = x_pre_del_mu / wtD;
                y_pre_del_mu  = yb_les_muA * double(wt);
                y_muD_les_muA = y_pre_del_mu / wtD;
            } else {
                x_muD_les_muA = xb_les_muA / wtD;
                y_muD_les_muA = yb_les_muA / wtD;
            }
            m_xx[1] += x_muD_les_muA;
            m_xx[2] += y_muD_les_muA;
            // the mean is computed. drop out if ord==1
            y_post_diff = (yval - m_xx[2]);
            if (has_wts) {
                m_xx[3] += x_pre_del_mu * (xval - m_xx[1]);
                m_xx[4] += x_pre_del_mu * y_post_diff;
                m_xx[5] += y_pre_del_mu * y_post_diff;
            } else {
                m_xx[3] += xb_les_muA * (xval - m_xx[1]);
                m_xx[4] += xb_les_muA * y_post_diff;
                m_xx[5] += yb_les_muA * y_post_diff;
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline TwoWelford& rem_one (const double xval, const double yval, const W wt) {
            if (na_rm) {
                if (ISNAN(xval) || ISNAN(yval)) { return *this; }
                if (has_wts) {
                    if (ISNAN(wt) || (wt <= 0)) {
                        return *this;
                    }
                }
            }
            m_subc++;

            double xc_les_muA, x_pre_del_mu, x_muD_les_muA, 
                   yc_les_muA, y_pre_del_mu, y_muD_les_muA, 
                   y_post_diff,
                   wtD, wtA;
            // xval = x_c
            // wt = w_c
            // xc_les_muA = x_c - mu_A
            if (has_wts) {
                m_nel--; 
                wtA = double(m_wsum.as());
                m_wsum -= wt;
                wtD = double(m_wsum.as());
            } else {
                wtA = double(m_nel);
                m_nel--;
                wtD = double(m_nel);
            }
            if (wtD > 0) {
                xc_les_muA = xval - m_xx[1];
                yc_les_muA = yval - m_xx[2];

                if (has_wts) {
                    x_pre_del_mu  = xc_les_muA * double(wt);
                    x_muD_les_muA = - x_pre_del_mu / wtD;
                    y_pre_del_mu  = yc_les_muA * double(wt);
                    y_muD_les_muA = - y_pre_del_mu / wtD;
                } else {
                    x_muD_les_muA = - xc_les_muA / wtD;
                    y_muD_les_muA = - yc_les_muA / wtD;
                }
                m_xx[1] += x_muD_les_muA;
                m_xx[2] += y_muD_les_muA;
                y_post_diff = (yval - m_xx[2]);
                // the mean is computed. drop out if ord==1
                if (has_wts) {
                    m_xx[3] -= x_pre_del_mu * (xval - m_xx[1]);
                    m_xx[4] -= x_pre_del_mu * y_post_diff;
                    m_xx[5] -= y_pre_del_mu * y_post_diff;
                } else {
                    m_xx[3] -= xc_les_muA * (xval - m_xx[1]);
                    m_xx[4] -= xc_les_muA * y_post_diff;
                    m_xx[5] -= yc_les_muA * y_post_diff;
                }
            } else {
                // zero it out?
                m_wsum = W(0);
                m_nel = 0;
                for (int iii=0;iii < 6;++iii) { m_xx[iii] = 0; }
            }
            return *this;
        }
        inline TwoWelford& swap_one (const double addxval, const double addyval, const W addwt,
                                  const double remxval, const double remyval, const W remwt) {
            add_one(addxval,addyval,addwt);
            rem_one(remxval,remyval,remwt);
            return *this;
        }
        // join two TwoWelford objects together
        inline TwoWelford& join(const TwoWelford& rhs) {
            double n1, n2, ntot, 
            x_del21, x_mupart, x_nfoo, 
            y_del21, y_mupart, y_nfoo, 
            n1rat, n2rat;
            double ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
            int ppp,qqq;
            if (has_wts) {
                n1 = double(m_wsum.as());
            } else {
                n1 = double(m_nel);
            }
            if (n1 <= 0) {
                // lhs is empty; just copy the rhs.
                m_nel = rhs.m_nel;
                m_wsum = rhs.m_wsum;
                m_subc = rhs.m_subc;
                // clone it?
                for (int iii=1;iii < 6;++iii) { m_xx[iii] = rhs.m_xx[iii]; }
                return *this;
            }
            // problem: can we join a weighted and unweighted welford object together? ack.
            if (has_wts) {
                n2 = double(rhs.m_wsum.as());
            } else {
                n2 = double(rhs.m_nel);
            }
            if (n2 <= 0) {
                // rhs is empty; just return the lhs.
                return *this;
            }
            // else onboard the observations
            m_nel += rhs.m_nel; 
            m_wsum += rhs.m_wsum;
            m_subc += rhs.m_subc;

            //ntot = double(m_wsum.as());
            ntot = n1 + n2;
            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            ac_n2 = 1.0/n2;
            ac_mn1 = 1.0/n1;

            x_del21 = rhs.m_xx[1] - m_xx[1];
            x_mupart = x_del21 * n2rat;
            y_del21 = rhs.m_xx[2] - m_xx[2];
            y_mupart = y_del21 * n2rat;

            m_xx[1] += x_mupart;
            m_xx[2] += y_mupart;

            x_nfoo = n1 * x_mupart;
            y_nfoo = n1 * y_mupart;

            m_xx[3] += rhs.m_xx[3] + (x_nfoo * x_nfoo * (ac_n2 + ac_mn1));
            // spit balling here, not sure.
            m_xx[4] += rhs.m_xx[4] + (x_nfoo * y_nfoo * (ac_n2 + ac_mn1));
            m_xx[5] += rhs.m_xx[5] + (y_nfoo * y_nfoo * (ac_n2 + ac_mn1));

            return *this;
        }
        // remove one from another
        inline TwoWelford& unjoin(const TwoWelford& rhs) {
            double n1, n2, ntot, 
                   x_del21, x_mupart, x_nfoo, 
                   y_del21, y_mupart, y_nfoo, 
                   n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1;
            double ac_del,ac_mn2,ac_n1;
            int ppp,qqq;

            // problem: can we join a weighted and unweighted welford object together? ack.
            if (has_wts) {
                ntot = double(m_wsum.as());
                n2 = double(rhs.m_wsum.as());
            } else {
                ntot = double(m_nel);
                n2 = double(rhs.m_nel);
            }

            if (n2 <= 0) { 
                // rhs is empty; return lhs.
                return *this; 
            }
            if (n2 > ntot) { stop("cannot subtract more observations than were seen."); } // #nocov

            x_mupart = rhs.m_xx[1] - m_xx[1];
            y_mupart = rhs.m_xx[2] - m_xx[2];

            m_nel -= rhs.m_nel; 
            if (has_wts) {
                m_wsum -= rhs.m_wsum;
                n1 = double(m_wsum.as());
            } else {
                n1 = double(m_nel);
            }
            // crap. this is just wrong...
            m_subc += rhs.m_subc;

            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            ac_n2 = 1.0 / n2;
            ac_mn1 = 1.0 / n1;

            m_xx[1] -= (n2/n1) * x_mupart;
            m_xx[2] -= (n2/n1) * y_mupart;

            x_del21 = x_mupart / n1rat;
            y_del21 = y_mupart / n1rat;
            x_nfoo = x_mupart * n2;
            y_nfoo = y_mupart * n2;

            m_xx[3] -= rhs.m_xx[3] + (x_nfoo * x_nfoo * (ac_n2 + ac_mn1));
            // spit balling here, not sure.
            m_xx[4] -= rhs.m_xx[4] + (y_nfoo * x_nfoo * (ac_n2 + ac_mn1));
            m_xx[5] -= rhs.m_xx[5] + (y_nfoo * y_nfoo * (ac_n2 + ac_mn1));
            return *this;
        }

};
//UNFOLD

// univariate sums, moments, cumulants//FOLDUP

template <typename T,typename W,typename oneW,bool has_wts,bool na_rm>
void add_many(TwoWelford<oneW,has_wts,na_rm> & frets,
              T v,
              T vv,
              W wts,
              int bottom,
              int top,
              const bool check_wts) {
    double nextwt;

    if (!has_wts) { nextwt = 1.0; }
    if ((top < 0) || (top > v.size())) { top = v.size(); }
    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); } // #nocov
        if (wts.size() < top) { stop("size of wts does not match v"); } // #nocov
    }
    for (int iii=bottom;iii < top;++iii) {
        if (has_wts) { 
            nextwt = double(wts[iii]); 
        }
        frets.add_one(v[iii],vv[iii],nextwt);
    }
}

template <typename T,typename W,typename oneW,bool has_wts,bool na_rm>
TwoWelford<oneW,has_wts,na_rm> quasiWeightedThing(T v,
                                                  T vv,
                                                  W wts,
                                                  int bottom,
                                                  int top,
                                                  const bool check_wts) {
    TwoWelford<oneW,has_wts,na_rm> frets = TwoWelford<oneW,has_wts,na_rm>();
    add_many<T,W,oneW,has_wts,na_rm>(frets,v,vv,wts,bottom,top,check_wts);
    return frets;
}

// this function returns a NumericVector of:
//   the number of elements or the sum of wts, 
//   the mean, and 
//   an (ord - 1)-vector consisting of the
//   2nd through ord'th centered sum, defined
//   as sum_j wts[j] * (v[j] - mean)^i
// if top < 0, take the length of v.
//
// if normalize_wts and wts are non-null, then
// we essentially renormalize the weights to
// have mean 1. in that case, the zeroth element
// returned is the 

template <typename T,typename W,typename oneW,bool has_wts,bool na_rm>
NumericVector quasiWeightedMoments(T v,
                                   T vv,
                                   W wts,
                                   int bottom,
                                   int top,
                                   const bool check_wts,
                                   const bool normalize_wts) {
    double nextv, nextw, renorm, nok;
    NumericVector xret;

    TwoWelford<oneW,has_wts,na_rm> irets = quasiWeightedThing<T,W,oneW,has_wts,na_rm>(v,vv,wts,bottom,top,check_wts);
    xret = irets.asvec();
    xret[0] = double(irets.wsum());
    nok = double(irets.nel());

    if (has_wts && normalize_wts) {//FOLDUP
        renorm = nok / xret[0];
        xret[0] = nok;
        for (int ppp=3;ppp <= 5;ppp++) {
            xret[ppp] *= renorm;
        }
    }//UNFOLD
    return xret;
}

// wrap one level
template <typename T,typename W,typename oneW,bool has_wts>
NumericVector quasiWeightedMomentsCurryZero(T v, 
                                            T vv,
                                            W wts,
                                            int bottom,
                                            int top,
                                            const bool na_rm,
                                            const bool check_wts,
                                            const bool normalize_wts) {

    if (na_rm) {
        return quasiWeightedMoments<T,W,oneW,has_wts,true>(v, vv, wts, bottom, top, check_wts, normalize_wts); 
    } 
    // have to have fallthrough for CRAN check.
    return quasiWeightedMoments<T,W,oneW,has_wts,false>(v, vv, wts, bottom, top, check_wts, normalize_wts); 
}

// wrap one level
// fix ord_beyond
template <typename T>
NumericVector quasiWeightedMomentsCurryOne(T v, 
                                           T vv,
                                           SEXP wts, 
                                           const bool na_rm, 
                                           const bool check_wts, 
                                           const bool normalize_wts) {
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return quasiWeightedMomentsCurryZero<T,IntegerVector,int,true>(v, vv, wts, 0, -1, na_rm, check_wts, normalize_wts); }
            case REALSXP: { return quasiWeightedMomentsCurryZero<T,NumericVector,double,true>(v, vv, wts, 0, -1, na_rm, check_wts, normalize_wts); }
            case  LGLSXP: { return quasiWeightedMomentsCurryZero<T,IntegerVector,int,true>(v, vv, as<IntegerVector>(wts), 0, -1, na_rm, check_wts, normalize_wts); } // bools can be upcast to save build size.
            default: stop("Unsupported weight type"); // #nocov
        }
    }
    NumericVector dummy_wts;
    // have to have fallthrough for CRAN check.
    return quasiWeightedMomentsCurryZero<T,NumericVector,int,false>(v, vv, dummy_wts, 0, -1, na_rm, check_wts, normalize_wts); 
}

// wrap one level
NumericVector inline quasiWeightedMomentsCurryTwo(SEXP v, 
                                                  SEXP vv,
                                                  SEXP wts, 
                                                  const bool na_rm, 
                                                  const bool check_wts, 
                                                  const bool normalize_wts) {
    if (!Rf_isNull(v)) {  
        switch (TYPEOF(v)) {
            case  INTSXP: { return quasiWeightedMomentsCurryOne<IntegerVector>(v, vv, wts, na_rm, check_wts, normalize_wts); }
            case REALSXP: { return quasiWeightedMomentsCurryOne<NumericVector>(v, vv, wts, na_rm, check_wts, normalize_wts); }
            case  LGLSXP: { return quasiWeightedMomentsCurryOne<IntegerVector>(as<IntegerVector>(v), as<IntegerVector>(vv), wts, na_rm, check_wts, normalize_wts); }  // bools can be upcast to save build size.
            default: stop("Unsupported data type"); // #nocov
        }
    }
    // have to have fallthrough for CRAN check.
    NumericVector retv(6);
    return retv;
}

//UNFOLD

#endif /* __DEF_TWO_WELFORD__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
