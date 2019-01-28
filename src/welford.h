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
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_WELFORD__
#define __DEF_WELFORD__

#include "common.h"
#include "kahan.h"

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// Welford-Terriberry object.
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
// m_ord : the integer order of the summation.
//
// a vector xx of length m_ord;
// xx[0] is ignored. keep it zero.
// xx[1] is weighted mean of x.
// xx[2] is weighted sum of (x - xx[1])^2
// xx[3] is weighted sum of (x - xx[1])^3
// xx[4] is weighted sum of (x - xx[1])^4
// ...

// Welford Terriberry
// generic//FOLDUP
// ord_beyond must be used for (ord > 2)
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
template<class W,bool has_wts,bool ord_beyond,bool na_rm>
class Welford {
    public:
        int m_ord;
    private: 
        int m_nel;
        int m_subc;
    private:
        Kahan<W> m_wsum;
    public:
        NumericVector m_xx;
    public:
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_subc(0), m_wsum(Kahan<W>(0)), m_xx(NumericVector(ord+1)) {
            if (!ord_beyond) {
                if (ord < 2) { stop("must use ord >= 2"); } // #nocov
            } else {
                if (ord < 1) { stop("must use ord >= 1"); } // #nocov
            }
        }
        inline Welford(const int &ord, 
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_subc(0), m_wsum(Kahan<W>(sumwt)), m_xx(NumericVector(xx)) {
            if (!ord_beyond) {
                if (ord < 2) { stop("must use ord >= 2"); } // #nocov
            } else {
                if (ord < 1) { stop("must use ord >= 1"); } // #nocov
            }
        }
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_subc(0), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {
            if (!ord_beyond) {
                if (ord < 2) { stop("must use ord >= 2"); } // #nocov
            } else {
                if (ord < 1) { stop("must use ord >= 1"); } // #nocov
            }
        }
    public:
        // reset to zero
        inline Welford& tare() {
            m_nel = 0;
            m_subc = 0;
            m_wsum = W(0);
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const bool normalize,const double used_df) const {
            double renorm;
            if (has_wts) {
                if (normalize) {
                    renorm = double(m_nel) / double(m_wsum.as());
                    return ((renorm * m_xx[2]) / (double(m_nel) - used_df));
                } else {
                    return ((m_xx[2]) / (double(m_wsum.as()) - used_df));
                }
            } else {
                return ((m_xx[2]) / (double(m_nel) - used_df));
            }
        }
        inline double mean() const {
            return m_xx[1];
        }
        inline double sd(const bool normalize,const double used_df) const {
            return sqrt(var(normalize,used_df));
        }
        // one centered moment, not standardized.
        // this should be the approximate weighted mean of
        // (x - weighted_mean)^ord
        inline double a_cent_mom(const int ord,const bool normalize,const double used_df) const {
            if (ord==2) {
                return sqrt(var(normalize,used_df));
            } else if (ord==1) {
                // Rcpp::warning("first centered moment is zero."); // #nocov
                return 0.0;
            } else {
                if (ord > m_ord) { stop("cannot compute this high of a moment."); } // #nocov
                if (has_wts) {
                    return m_xx[ord] / double(m_wsum.as());
                } else {
                    return m_xx[ord] / double(m_nel);
                }
            }
        }
        inline double skew() const {
            if (has_wts) {
                return (sqrt(double(m_wsum.as())) * m_xx[3] / pow(m_xx[2],1.5));
            } else {
                return (sqrt(double(m_nel)) * m_xx[3] / pow(m_xx[2],1.5));
            }
        }
        inline double exkurt() const {
            if (has_wts) {
                return ((double(m_wsum.as()) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
            } else {
                return ((double(m_nel) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
            }
        }
        inline double sharpe(const bool normalize,const double used_df) const {
            return double(m_xx[1]) / sd(normalize,used_df);
        }
        inline double centered(const double xval) const {
            return (xval - m_xx[1]);
        }
        inline double scaled(const double xval,const bool normalize,const double used_df) const {
            return (xval/sd(normalize,used_df));
        }
        inline double zscored(const double xval,const bool normalize,const double used_df) const {
            return ((xval - m_xx[1])/sd(normalize,used_df));
        }

        // getters 
        inline int nel() const { if (has_wts) { return m_nel; } else { return int(wsum()); } }  // not sure I understand this...
        inline int subcount() const { return m_subc; }

        // return true if any even order sums are negative
        inline bool has_heywoord() const {
            if (!ord_beyond) {
                return (m_xx[2] < 0);
            } else {
                for (int ppp=2;ppp <= m_ord;ppp += 2) {
                    if (m_xx[ppp] < 0) { return true; }
                }
                return false;
            }
            return false;
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
        inline Welford& add_one (const double xval, const W wt) {
            if (na_rm) {
                if (ISNAN(xval)) { return *this; }
                if (has_wts) {
                    if (ISNAN(wt) || (wt <= 0)) {
                        return *this;
                    }
                }
            }

            double xb_les_muA, pre_del_mu, muD_les_muA, wtD, wtA;
            double term_left, div_left, rem_right, div_right, inner_term;
            // xval = x_b
            // wt = w_b
            // xb_les_muA = x_b - mu_A

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
            if (has_wts) {
                pre_del_mu  = xb_les_muA * double(wt);
                muD_les_muA = pre_del_mu / wtD;
            } else {
                muD_les_muA = xb_les_muA / wtD;
            }
            m_xx[1] += muD_les_muA;
            // the mean is computed. drop out if ord==1
            if (!ord_beyond) {
                if (has_wts) {
                    m_xx[2] += pre_del_mu * (xval - m_xx[1]);
                } else {
                    m_xx[2] += xb_les_muA * (xval - m_xx[1]);
                }
            } else {
                if ((wtA > 0) && (m_ord > 1)) {
                    div_left = -muD_les_muA;
                    term_left = pow(div_left,m_ord) * wtA;
                    if (has_wts) {
                        div_right = -wtA / double(wt);
                    } else {
                        div_right = -wtA;
                    }
                    rem_right = pow(div_right,m_ord - 1);

                    for (int ppp=m_ord;ppp > 2;ppp--) {
                        m_xx[ppp] += term_left * (1.0 - rem_right);
                        // could hit division by zero here ? 
                        term_left /= div_left;
                        rem_right /= div_right;
                        inner_term = div_left;
                        for (int qqq=1;qqq <= ppp-2;qqq++) {
                            m_xx[ppp] += bincoef[ppp][qqq] * inner_term * m_xx[ppp-qqq];
                            if (qqq < ppp - 2) { inner_term *= div_left; }
                        }
                    }
                    m_xx[2] += term_left * (1.0 - rem_right);
                }
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            if (na_rm) {
                if (ISNAN(xval)) { return *this; }
                if (has_wts) {
                    if (ISNAN(wt) || (wt <= 0)) {
                        return *this;
                    }
                }
            }
            m_subc++;

            double xc_les_muA, pre_del_mu, muD_les_muA, wtD, wtA;
            double term_left, div_left, rem_right, div_right, inner_term;

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

                if (has_wts) {
                    pre_del_mu  = xc_les_muA * double(wt);
                    muD_les_muA = - pre_del_mu / wtD;
                } else {
                    muD_les_muA = - xc_les_muA / wtD;
                }
                m_xx[1] += muD_les_muA;
                // the mean is computed. drop out if ord==1
                if (!ord_beyond) {
                    if (has_wts) {
                        m_xx[2] -= pre_del_mu * (xval - m_xx[1]);
                    } else {
                        m_xx[2] -= xc_les_muA * (xval - m_xx[1]);
                    }
                } else {
                    div_left = -muD_les_muA;
                    term_left = pow(div_left,m_ord) * wtA;
                    if (has_wts) {
                        div_right = wtA / double(wt);
                    } else {
                        div_right = wtA;
                    }
                    rem_right = pow(div_right,m_ord - 1);

                    for (int ppp=m_ord;ppp >= 2;ppp--) {
                        m_xx[ppp] += term_left * (1.0 - rem_right);
                        if (ppp > 2) {
                            // could hit division by zero here ? 
                            // in fact, you will for the first value encountered.
                            term_left /= div_left;
                            rem_right /= div_right;
                            inner_term = div_left;
                            for (int qqq=1;qqq <= ppp-2;qqq++) {
                                m_xx[ppp] += bincoef[ppp][qqq] * inner_term * m_xx[ppp-qqq];
                                if (qqq < ppp - 2) { inner_term *= div_left; }
                            }
                        }
                    }
                }
            } else {
                // zero it out?
                m_wsum = W(0);
                m_nel = 0;
                if (!ord_beyond) {
                    m_xx[1] = 0.0;
                    m_xx[2] = 0.0;
                } else {
                    for (int ppp=1;ppp <= m_ord;ppp++) {
                        m_xx[ppp] = 0.0;
                    }
                }
            }
            return *this;
        }
        inline Welford& swap_one (const double addxval, const W addwt,
                                  const double remxval, const W remwt) {
            // na checking
            if (na_rm) {
                if (ISNAN(addxval)) {
                    if (ISNAN(remxval)) {
                        return *this;
                    } else {
                        rem_one(remxval,remwt);
                        return *this;
                    }
                } else if (ISNAN(remxval)) {
                    add_one(addxval,addwt);
                    return *this;
                }
                if (has_wts) {
                    if (ISNAN(addwt) || (addwt <= 0)) {
                        if (ISNAN(remwt) || (remwt <= 0)) {
                            return *this;
                        } else {
                            rem_one(remxval,remwt);
                            return *this;
                        }
                    } else if (ISNAN(remwt) || (remwt <= 0)) {
                        add_one(addxval,addwt);
                        return *this;
                    }
                }
            }
            m_subc++;

            double diffmu,prevmu,nel;
            double diffw,diffx,diffxw,addxw,remxw,nelm;
            if (!ord_beyond) {
                if (has_wts) {
                    // try this instead, maybe my math is wrong below...
                    add_one(addxval,addwt);
                    rem_one(remxval,remwt);
                    // yuck; maybe just call add_one and rem_one instead?
                    
                    //nelm = double(m_wsum.as());
                    //addxw = addxval * double(addwt);
                    //remxw = remxval * double(remwt);
                    //diffw = double(addwt) - double(remwt);
                    //diffx = addxval - remxval;
                    //diffxw = addxw - remxw;
                    //diffmu = m_xx[1] * diffw + diffxw;

                    //m_wsum += diffw;
                    //nel = double(m_wsum.as());
                    //// 2FIX: check for bottoming out?
                    //prevmu = m_xx[1];
                    //m_xx[1] += (diffmu/nel);
                    //m_xx[2] += (nelm * (-prevmu * (diffmu + diffxw) + (addxw * addxval - remxw * remxval)) - double(addwt) * double(remwt) * diffx * diffx) / nel;
                } else {
                    nel = double(m_nel);
                    diffmu = addxval - remxval;
                    prevmu = m_xx[1];
                    m_xx[1] += (diffmu/nel);
                    m_xx[2] += diffmu*(addxval + remxval - prevmu - m_xx[1]);
                }
            } else {
                // too hard for ord > 2 case;
                add_one(addxval,addwt);
                rem_one(remxval,remwt);
            }
            return *this;
        }
        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
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
                if (!ord_beyond) {
                    m_xx[1] = rhs.m_xx[1];
                    m_xx[2] = rhs.m_xx[2];
                } else {
                    m_xx = Rcpp::clone(rhs.m_xx);
                }
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
            del21 = rhs.m_xx[1] - m_xx[1];
            mupart = del21 * n2rat;

            m_xx[1] += mupart;
            nfoo = n1 * mupart;
            if (!ord_beyond) {
                ac_nfoo = nfoo*nfoo;
                ac_n2 = 1.0/n2;
                ac_mn1 = 1.0/n1;
                m_xx[2] += rhs.m_xx[2] + (ac_nfoo * (ac_n2 + ac_mn1));
            } else {
                ac_nfoo = pow(nfoo,m_ord);
                ac_n2 = pow(n2,1-m_ord);
                ac_mn1 = pow(-n1,1-m_ord);

                for (ppp=m_ord;ppp >= 2;ppp--) {
                    m_xx[ppp] += rhs.m_xx[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
                    if (ord_beyond) {
                        if (ppp > 2) {
                            if (nfoo != 0) { ac_nfoo /= nfoo; }
                            ac_n2 *= n2;
                            ac_mn1 *= (-n1);

                            ac_del = del21;
                            ac_mn2 = -n2rat;
                            ac_n1 = n1rat;
                            for (int qqq=1;qqq <= (ppp-2); qqq++) {
                                m_xx[ppp] += bincoef[ppp][qqq] * ac_del * (ac_mn2 * m_xx[ppp-qqq] + ac_n1 * rhs.m_xx[ppp-qqq]);
                                if (qqq < (ppp-2)) {
                                    ac_del *= del21;
                                    ac_mn2 *= (-n2rat);
                                    ac_n1  *= (n1rat);
                                }
                            }
                        }
                    }
                }
            }
            return *this;
        }
        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
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

            mupart = rhs.m_xx[1] - m_xx[1];

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

            m_xx[1] -= (n2/n1) * mupart;

            del21 = mupart / n1rat;
            nfoo = mupart * n2;

            if (!ord_beyond) {
                ac_nfoo = nfoo * nfoo;
                ac_n2 = 1.0 / n2;
                ac_mn1 = 1.0 / n1;
                m_xx[2] -= rhs.m_xx[2] + (ac_nfoo * (ac_n2 + ac_mn1));
            } else {
                ac_nfoo = nfoo * nfoo;
                ac_n2 = 1.0 / n2;
                ac_mn1 = -1.0 / n1;
                for (ppp=2;ppp <= m_ord;ppp++) {
                    m_xx[ppp] -= rhs.m_xx[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
                    if (ord_beyond) {
                        if (ppp < m_ord) {
                            ac_nfoo *= nfoo; 
                            ac_n2 /= n2;
                            ac_mn1 /= (-n1);
                        }
                        ac_del = del21;
                        ac_mn2 = -n2rat;
                        ac_n1 = n1rat;
                        for (int qqq=1;qqq <= (ppp-2); qqq++) {
                            m_xx[ppp] -= bincoef[ppp][qqq] * ac_del * (ac_mn2 * m_xx[ppp-qqq] + ac_n1 * rhs.m_xx[ppp-qqq]);
                            if (qqq < (ppp-2)) {
                                ac_del *= del21;
                                ac_mn2 *= (-n2rat);
                                ac_n1  *= (n1rat);
                            }
                        }
                    }
                }
            }
            return *this;
        }

};
//UNFOLD

// univariate sums, moments, cumulants//FOLDUP

// quasiSumThing :
// given weights and values, computes the
// sum of weights, and the weighted mean of the values.
// computes from bottom to top; if top is negative, change to the size of v;
// note that top should be like numel. we sum over indices iii such that
// bottom <= iii < top
// if it happens to be the case that bottom == top, we should return a tared Welford object.

template <typename T,typename W,typename oneW,bool has_wts,bool na_rm>
NumericVector quasiSumThing(T v,
                            W wts,
                            int bottom,
                            int top,
                            const bool check_wts,
                            const bool normalize_wts) {
    double nextv, nextw;
    Kahan<double> fwvsum;
    Kahan<oneW> fwsum;
    double totwt;
    int nel = 0;

    if ((top < 0) || (top > v.size())) { top = v.size(); }
    if (has_wts) {
        if (wts.size() < top) { stop("size of wts does not match v"); }
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); }
        //2FIX: push na_rm into template params?
        if (na_rm) {
            for (int iii=bottom;iii < top;++iii) {
                nextv = v[iii];
                nextw = double(wts[iii]); 
                if (! (ISNAN(nextv) || ISNAN(nextw))) {
                    // 2FIX: check for zero weight??
                    fwvsum += nextv * nextw;
                    fwsum += nextw;
                    ++nel;
                }
            }
        } else {
            for (int iii=bottom;iii < top;++iii) {
                nextv = v[iii];
                nextw = double(wts[iii]); 
                fwvsum += nextv * nextw;
                fwsum += nextw;
                ++nel;
            }
        }
    } else {
        if (na_rm) {
            for (int iii=bottom;iii < top;++iii) {
                nextv = v[iii];
                if (! (ISNAN(nextv))) { 
                    fwvsum += nextv; 
                    ++fwsum;
                }
            }
        } else {
            for (int iii=bottom;iii < top;++iii) {
                nextv = v[iii];
                fwvsum += nextv; 
                ++fwsum;
            }
        }
    }
    totwt = double(fwsum.as());
    NumericVector vret = NumericVector::create(totwt,double(fwvsum.as()) / totwt);
    // the mean does not change, but the 'sum weights' becomes the number of elements
    if (has_wts && normalize_wts) {
        vret[0] = double(nel);
    }
    return vret;
}

// ord is largely ignred here;
template <typename T,typename W,typename oneW,bool has_wts,bool ord_beyond,bool na_rm>
void add_many(Welford<oneW,has_wts,ord_beyond,na_rm> & frets,
              T v,
              W wts,
              int ord,
              int bottom,
              int top,
              const bool check_wts) {
    double nextval, nextwt;

    if (!has_wts) { nextwt = 1.0; }
    if ((top < 0) || (top > v.size())) { top = v.size(); }
    if (has_wts) {
        if (check_wts && bad_weights<W>(wts)) { stop("negative weight detected"); } // #nocov
        if (wts.size() < top) { stop("size of wts does not match v"); } // #nocov
    }
    for (int iii=bottom;iii < top;++iii) {
        nextval = v[iii];
        if (has_wts) { 
            nextwt = double(wts[iii]); 
        }
        frets.add_one(nextval,nextwt);
    }
}

template <typename T,typename W,typename oneW,bool has_wts,bool ord_beyond,bool na_rm>
Welford<oneW,has_wts,ord_beyond,na_rm> quasiWeightedThing(T v,
                                                          W wts,
                                                          int ord,
                                                          int bottom,
                                                          int top,
                                                          const bool check_wts) {
    Welford<oneW,has_wts,ord_beyond,na_rm> frets = Welford<oneW,has_wts,ord_beyond,na_rm>(ord);
    add_many<T,W,oneW,has_wts,ord_beyond,na_rm>(frets,v,wts,ord,bottom,top,check_wts);
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
                                   W wts,
                                   int ord,
                                   int bottom,
                                   int top,
                                   const bool check_wts,
                                   const bool normalize_wts) {
    double nextv, nextw, renorm, nok;
    if (ord < 1) { stop("require positive order"); } // #nocov
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); } // #nocov
    NumericVector xret;

    if (ord == 1) {
        //2FIX: no normalization??
        xret = quasiSumThing<T,W,oneW,has_wts,na_rm>(v,wts,bottom,top,check_wts,normalize_wts);
        return xret;
    } else if (ord > 2) {
        Welford<oneW,has_wts,true,na_rm> frets = quasiWeightedThing<T,W,oneW,has_wts,true,na_rm>(v,wts,ord,bottom,top,check_wts);
        xret = frets.asvec();
        nok = double(frets.nel());
    } else {
        Welford<oneW,has_wts,false,na_rm> irets = quasiWeightedThing<T,W,oneW,has_wts,false,na_rm>(v,wts,ord,bottom,top,check_wts);
        xret = irets.asvec();
        xret[0] = double(irets.wsum());
        nok = double(irets.nel());
    }

    if (has_wts && normalize_wts) {//FOLDUP
        renorm = nok / xret[0];
        xret[0] = nok;
        for (int ppp=2;ppp <= ord;ppp++) {
            xret[ppp] *= renorm;
        }
    }//UNFOLD
    return xret;
}

// wrap one level
template <typename T,typename W,typename oneW,bool has_wts>
NumericVector quasiWeightedMomentsCurryZero(T v, 
                                            W wts,
                                            int ord,
                                            int bottom,
                                            int top,
                                            const bool na_rm,
                                            const bool check_wts,
                                            const bool normalize_wts) {

    if (na_rm) {
        return quasiWeightedMoments<T,W,oneW,has_wts,true>(v, wts, ord, bottom, top, check_wts, normalize_wts); 
    } 
    // have to have fallthrough for CRAN check.
    return quasiWeightedMoments<T,W,oneW,has_wts,false>(v, wts, ord, bottom, top, check_wts, normalize_wts); 
}

// wrap one level
// fix ord_beyond
template <typename T>
NumericVector quasiWeightedMomentsCurryOne(T v, 
                                           SEXP wts, 
                                           int ord, 
                                           const bool na_rm, 
                                           const bool check_wts, 
                                           const bool normalize_wts) {
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return quasiWeightedMomentsCurryZero<T,IntegerVector,int,true>(v, wts, ord, 0, -1, na_rm, check_wts, normalize_wts); }
            case REALSXP: { return quasiWeightedMomentsCurryZero<T,NumericVector,double,true>(v, wts, ord, 0, -1, na_rm, check_wts, normalize_wts); }
            case  LGLSXP: { return quasiWeightedMomentsCurryZero<T,IntegerVector,int,true>(v, as<IntegerVector>(wts), ord, 0, -1, na_rm, check_wts, normalize_wts); } // bools can be upcast to save build size.
            default: stop("Unsupported weight type"); // #nocov
        }
    }
    NumericVector dummy_wts;
    // have to have fallthrough for CRAN check.
    return quasiWeightedMomentsCurryZero<T,NumericVector,int,false>(v, dummy_wts, ord, 0, -1, na_rm, check_wts, normalize_wts); 
}

// wrap one level
NumericVector inline quasiWeightedMomentsCurryTwo(SEXP v, 
                                                  SEXP wts, 
                                                  int ord, 
                                                  const bool na_rm, 
                                                  const bool check_wts, 
                                                  const bool normalize_wts) {
    if (!Rf_isNull(v)) {  
        switch (TYPEOF(v)) {
            case  INTSXP: { return quasiWeightedMomentsCurryOne<IntegerVector>(v, wts, ord, na_rm, check_wts, normalize_wts); }
            case REALSXP: { return quasiWeightedMomentsCurryOne<NumericVector>(v, wts, ord, na_rm, check_wts, normalize_wts); }
            case  LGLSXP: { return quasiWeightedMomentsCurryOne<IntegerVector>(as<IntegerVector>(v), wts, ord, na_rm, check_wts, normalize_wts); }  // bools can be upcast to save build size.
            default: stop("Unsupported data type"); // #nocov
        }
    }
    // have to have fallthrough for CRAN check.
    NumericVector retv(ord+1);
    return retv;
}

//UNFOLD

#endif /* __DEF_WELFORD__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
