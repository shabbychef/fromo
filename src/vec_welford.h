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

  Created: 2024.11.09
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_VEC_WELFORD__
#define __DEF_VEC_WELFORD__

#include "common.h"
#include "kahan.h"

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// Welford-Terriberry higher order object.
// This tracks sums and cosums of input vector x.
//
// Object holds the following:
//
// m_nx : the (integer) length of the vector x.
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
// a vector m_xx of length 
// 1 + m_nx + m_nx * (m_nx + 1)/2 + (m_nx + 2) choose 3 + ... 
//
// m_xx[0] is ignored. keep it zero.
// m_xx[1:m_nx] is weighted mean of x[1:m_nx].
// m_xx[m_nx+1:m_nx + (m_nx + 1) choose 2] is weighted sum of (x - m_xx[1:m_nx])^2
// m_xx[3] is weighted sum of (x - xx[1])^3
// m_xx[4] is weighted sum of (x - xx[1])^4
// ...
//
// Is that right? Think of the case where n = 2
// for the third order object we want to track
// x111 x112 x122 x222
// there are four of those. Is that 4 choose 3? Yes.
// How about n = 3
// for the third order object we have
// x111 x112 x113 x122 x123 x133 x222 x223 x233 x333
// there are 10 of those. Is that 5 choose 3? Yes.
//
// OK, but how do we index those?
//
// WARNING: at the moment we only support the case m_ord = 2.
// The higher order moments are a bit too complicated for me at the moment.

// Vector Welford Terriberry
//
// generic//FOLDUP
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
// going to have to fake ord_beyond and na_rm at the object level
template<class W,bool has_wts>
class VecWelford {
    public:
        int m_ord;
        int m_nx;
        bool na_rm;
    private: 
        int m_nel;
        int m_subc;
    private:
        Kahan<W> m_wsum;
    public:
        NumericVector m_xx;
    public:
        inline VecWelford(const int &ord, const int &nx) : m_ord(ord), m_nx(nx), m_nel(0), m_subc(0), m_wsum(Kahan<W>(0)) {
            if (ord < 1) { stop("must use ord >= 1"); } // #nocov
            if (ord > 2) { stop("NYI: must use ord <= 2"); } // #nocov
            int n_xx = 1 + nx;
            if (ord > 1) {
                n_xx += nx * (nx+1)/2;
                if (ord > 2) {
                    for (int iii=3;iii <= ord;++iii) {
                        // this could fail for large nx in theory.
                        n_xx += bincoef[iii + nx][iii + 1];
                    }
                }
            }
            m_xx = NumericVector(n_xx);
        }
        inline VecWelford(const int &ord, 
                         const int &nx,
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nx(nx), m_nel(nel), m_subc(0), m_wsum(Kahan<W>(sumwt)), m_xx(NumericVector(xx)) {
            if (ord < 1) { stop("must use ord >= 1"); } // #nocov
            if (ord > 2) { stop("NYI: must use ord <= 2"); } // #nocov
        }
    public:
        // reset to zero
        inline VecWelford& tare() {
            m_nel = 0;
            m_subc = 0;
            m_wsum = W(0);
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        // getters 
        inline int nel() const { if (has_wts) { return m_nel; } else { return int(wsum()); } }  // not sure I understand this...
        inline int subcount() const { return m_subc; }

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
        inline NumericVector vecpart() const { return m_xx; }
    public:
        // add another (weighted) observation to our set of x
        inline VecWelford& add_one (const NumericMatrix::Row xval, const W wt) {
            if (xval.size() != m_nx) { stop("gave row vector of inconsistent size") }
            if (na_rm) {
                for (int iii=0;iii<m_nx;++iii) {
                    if (ISNAN(xval[iii])) { return *this; }
                }
                if (has_wts) {
                    if (ISNAN(wt) || (wt <= 0)) {
                        return *this;
                    }
                }
            }

            NumericVector xb_les_muA, pre_del_mu, muD_les_muA;
            double wtD, wtA;
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
            xb_les_muA = xval - m_xx[1:m_nx];
            if (has_wts) {
                pre_del_mu  = xb_les_muA * double(wt);
                muD_les_muA = pre_del_mu / wtD;
            } else {
                muD_les_muA = xb_les_muA / wtD;
            }
            m_xx[1:m_nx] += muD_les_muA;
            // the mean is computed. drop out if ord==1
            // 2FIX: start here and do these right...
            if (has_wts) {
                m_xx[2] += pre_del_mu * (xval - m_xx[1]);
            } else {
                m_xx[2] += xb_les_muA * (xval - m_xx[1]);
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
                if (has_wts) {
                    m_xx[2] -= pre_del_mu * (xval - m_xx[1]);
                } else {
                    m_xx[2] -= xc_les_muA * (xval - m_xx[1]);
                }
            } else {
                // zero it out?
                m_wsum = W(0);
                m_nel = 0;
                m_xx[1] = 0.0;
                m_xx[2] = 0.0;
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
            // too hard for now; maybe some simplification possible later
            add_one(addxval,addwt);
            rem_one(remxval,remwt);
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

#endif /* __DEF_VEC_WELFORD__ */

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
