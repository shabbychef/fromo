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

#ifndef __DEF_FROMO__
#define __DEF_FROMO__

#include <math.h>

// preallocate the Binomial Coefficients, for efficiency

// this is R code used to generate C code. [ducks]

//MAXORD <- 30
//refv <- matrix(0,nrow=MAXORD,ncol=MAXORD)
//for (iii in seq(1,nrow(refv))) {
   //refv[iii,1] = 1; refv[iii,iii] = 1; 
   //if (iii > 2) { for (jjj in seq(2,iii-1)) { refv[iii,jjj] = refv[iii-1,jjj-1] + refv[iii-1,jjj]; } }
 //} 
//cat(sprintf('#define MAX_ORD %d\nconst int bincoef[%d][%d] = {',ncol(refv)-1,ncol(refv),ncol(refv)),
    //paste0('{\n',lapply(seq_len(nrow(refv)),function(rn) { paste0(sprintf('%9s',as.character(refv[rn,])),collapse=', ') }),'},\n'),
    //'};\n\n',file='/tmp/binc.txt')
   

#define MAX_ORD 29
const int bincoef[30][30] = {
 {        1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         2,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         3,         3,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         4,         6,         4,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         5,        10,        10,         5,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         6,        15,        20,        15,         6,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         7,        21,        35,        35,        21,         7,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         8,        28,        56,        70,        56,        28,         8,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,         9,        36,        84,       126,       126,        84,        36,         9,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        10,        45,       120,       210,       252,       210,       120,        45,        10,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        11,        55,       165,       330,       462,       462,       330,       165,        55,        11,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        12,        66,       220,       495,       792,       924,       792,       495,       220,        66,        12,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        13,        78,       286,       715,      1287,      1716,      1716,      1287,       715,       286,        78,        13,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        14,        91,       364,      1001,      2002,      3003,      3432,      3003,      2002,      1001,       364,        91,        14,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        15,       105,       455,      1365,      3003,      5005,      6435,      6435,      5005,      3003,      1365,       455,       105,        15,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        16,       120,       560,      1820,      4368,      8008,     11440,     12870,     11440,      8008,      4368,      1820,       560,       120,        16,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        17,       136,       680,      2380,      6188,     12376,     19448,     24310,     24310,     19448,     12376,      6188,      2380,       680,       136,        17,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        18,       153,       816,      3060,      8568,     18564,     31824,     43758,     48620,     43758,     31824,     18564,      8568,      3060,       816,       153,        18,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        19,       171,       969,      3876,     11628,     27132,     50388,     75582,     92378,     92378,     75582,     50388,     27132,     11628,      3876,       969,       171,        19,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        20,       190,      1140,      4845,     15504,     38760,     77520,    125970,    167960,    184756,    167960,    125970,     77520,     38760,     15504,      4845,      1140,       190,        20,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        21,       210,      1330,      5985,     20349,     54264,    116280,    203490,    293930,    352716,    352716,    293930,    203490,    116280,     54264,     20349,      5985,      1330,       210,        21,         1,         0,         0,         0,         0,         0,         0,         0,         0},
 {        1,        22,       231,      1540,      7315,     26334,     74613,    170544,    319770,    497420,    646646,    705432,    646646,    497420,    319770,    170544,     74613,     26334,      7315,      1540,       231,        22,         1,         0,         0,         0,         0,         0,         0,         0},
 {        1,        23,       253,      1771,      8855,     33649,    100947,    245157,    490314,    817190,   1144066,   1352078,   1352078,   1144066,    817190,    490314,    245157,    100947,     33649,      8855,      1771,       253,        23,         1,         0,         0,         0,         0,         0,         0},
 {        1,        24,       276,      2024,     10626,     42504,    134596,    346104,    735471,   1307504,   1961256,   2496144,   2704156,   2496144,   1961256,   1307504,    735471,    346104,    134596,     42504,     10626,      2024,       276,        24,         1,         0,         0,         0,         0,         0},
 {        1,        25,       300,      2300,     12650,     53130,    177100,    480700,   1081575,   2042975,   3268760,   4457400,   5200300,   5200300,   4457400,   3268760,   2042975,   1081575,    480700,    177100,     53130,     12650,      2300,       300,        25,         1,         0,         0,         0,         0},
 {        1,        26,       325,      2600,     14950,     65780,    230230,    657800,   1562275,   3124550,   5311735,   7726160,   9657700,  10400600,   9657700,   7726160,   5311735,   3124550,   1562275,    657800,    230230,     65780,     14950,      2600,       325,        26,         1,         0,         0,         0},
 {        1,        27,       351,      2925,     17550,     80730,    296010,    888030,   2220075,   4686825,   8436285,  13037895,  17383860,  20058300,  20058300,  17383860,  13037895,   8436285,   4686825,   2220075,    888030,    296010,     80730,     17550,      2925,       351,        27,         1,         0,         0},
 {        1,        28,       378,      3276,     20475,     98280,    376740,   1184040,   3108105,   6906900,  13123110,  21474180,  30421755,  37442160,  40116600,  37442160,  30421755,  21474180,  13123110,   6906900,   3108105,   1184040,    376740,     98280,     20475,      3276,       378,        28,         1,         0},
 {        1,        29,       406,      3654,     23751,    118755,    475020,   1560780,   4292145,  10015005,  20030010,  34597290,  51895935,  67863915,  77558760,  77558760,  67863915,  51895935,  34597290,  20030010,  10015005,   4292145,   1560780,    475020,    118755,     23751,      3654,       406,        29,         1},
 };
#define MAX(a,b) ((a>b)? (a):(b))
#define MIN(a,b) ((a<b)? (a):(b))

#endif /* __DEF_FROMO__ */


#include <Rcpp.h>
using namespace Rcpp;

#define COMP_SD_TWO(preval,used_df) (sqrt(preval[2]/(preval[0]-used_df)))
#define COMP_SD(preval) COMP_SD_TWO(preval,1.0)
#define COMP_SKEW(preval) (sqrt(preval[0]) * preval[3] / pow(preval[2],1.5))
#define COMP_EXKURT(preval) ((preval[0] * preval[4] / (pow(preval[2],2.0))) - 3.0)
#define COMP_SHARPE(preval) ((preval[1]) / (sqrt(preval[2] / (preval[0]-1.0))))

#define COMP_CENTERED(x,preval) (x - preval[1])

// Kahan compensated summation object.//FOLDUP
template<class T>
class Kahan {
    public:
        inline Kahan() : m_val(0), m_errs(0) {}

        inline Kahan(const T &t): m_val(t), m_errs(0) {}

        inline T as() const { return m_val; }

        inline Kahan& operator += (const T& rhs) {
            return add(rhs);
        }
        inline Kahan& operator -= (const T& rhs) {
            return add(-rhs);
        }
        inline Kahan& operator += (const Kahan<T>& rhs) {
            return join(rhs);
        }
        inline Kahan& operator -= (const Kahan<T>& rhs) {
            return unjoin(rhs);
        }
        inline Kahan& operator= (T rhs) {
            m_val = rhs;
            m_errs = 0;
            return *this;
        }
    public:
        T m_val;
    private:
        T m_errs;
        inline Kahan& add(const T& rhs) {
            T tmpv, nxtv;
            tmpv = rhs - m_errs;
            nxtv = m_val + tmpv;
            m_errs = (nxtv - m_val) - tmpv;
            m_val = nxtv;
            return *this;
        }
        inline Kahan& join(const Kahan<T>& rhs) {
            T tmpv, nxtv;
            tmpv = rhs.m_val - m_errs - rhs.m_errs;
            nxtv = m_val + tmpv;
            m_errs = (nxtv - m_val) - tmpv;
            m_val = nxtv;
            return *this;
        }
        inline Kahan& unjoin(const Kahan<T>& rhs) {
            T tmpv, nxtv;
            tmpv = -rhs.m_val - m_errs + rhs.m_errs;
            nxtv = m_val + tmpv;
            m_errs = (nxtv - m_val) - tmpv;
            m_val = nxtv;
            return *this;
        }
};
// specialization to int which do not require special treatment
template<>
class Kahan<int> {
    public:
        inline Kahan() : m_val(0) {}

        inline Kahan(const int &t): m_val(t) {}

        inline int as() const { return m_val; }

        inline Kahan& operator += (const int& rhs) {
            return add(rhs);
        }
        inline Kahan& operator -= (const int& rhs) {
            return add(-rhs);
        }
        inline Kahan& operator += (const Kahan<int>& rhs) {
            return add(rhs.m_val);
        }
        inline Kahan& operator -= (const Kahan<int>& rhs) {
            return add(-rhs.m_val);
        }
        inline Kahan& operator= (int rhs) {
            m_val = rhs;
            return *this;
        }
    public:
        int m_val;
    private:
        inline Kahan& add(const int& rhs) {
            m_val += rhs;
            return *this;
        }
};
//UNFOLD

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
// generic//FOLDUP
// ord_beyond must be used for (ord > 2)
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
template<class W,bool has_wts,bool ord_beyond>
class Welford {
    public:
        int m_ord;
        NumericVector m_xx;
    public:
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_wsum(Kahan<W>(0)), m_xx(NumericVector(ord+1)) {}
        inline Welford(const int &ord, 
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_wsum(Kahan<W>(sumwt)), m_xx(NumericVector(xx)) {}
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {}
    private:
        Kahan<W> m_wsum;
        int m_nel;
    public:
        // reset to zero
        inline Welford& tare() {
            if (has_wts) { m_nel = 0; }
            m_wsum = W(0);
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const bool normalize,const double used_df) const {
            double renorm;
            if (has_wts && normalize) {
                renorm = double(m_nel) / double(m_wsum.as());
                return ((renorm * m_xx[2]) / (double(m_nel) - used_df));
            } else {
                return ((m_xx[2]) / (double(m_wsum.as()) - used_df));
            }
        }
        inline double mean() const {
            return m_xx[1];
        }
        inline double sd(const bool normalize,const double used_df) const {
            return sqrt(var(normalize,used_df));
        }
        inline double skew() const {
            return (sqrt(double(m_wsum.as())) * m_xx[3] / pow(m_xx[2],1.5));
        }
        inline double exkurt() const {
            return ((double(m_wsum.as()) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
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
        inline int nel() const { if (has_wts) { return m_nel; } else { return wsum(); } }
        inline W wsum() const { return m_wsum.as(); }
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
            resu[0] = double(wsum());
            return resu;
        }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            if (has_wts) { m_nel++; }
            della = xval - m_xx[1];
            nelm = double(m_wsum.as());
            m_wsum += wt;
            nel = double(m_wsum.as());
            delnel = della * double(wt) / nel;
            m_xx[1] += delnel;
            if (nelm > 0) {
                drat = della * nelm / nel;
                nbyn = -wt / nelm;
                ac_dn = pow(drat,m_ord);
                ac_on = pow(nbyn,m_ord-1);

                for (int ppp=m_ord;ppp >= 2;ppp--) {
                    m_xx[ppp] += ac_dn * wt * (1.0 - ac_on);
                    if (ord_beyond) {
                        if (ppp > 2) {
                            if (drat != 0) { ac_dn /= drat; }
                            ac_on /= nbyn;
                            ac_de = -delnel;

                            for (int qqq=1;qqq <= ppp-2;qqq++) {
                                m_xx[ppp] += bincoef[ppp][qqq] * ac_de * m_xx[ppp-qqq];
                                if (qqq < ppp - 2) {
                                    ac_de *= -delnel;
                                }
                            }
                        }
                    }
                }
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            if (has_wts) { m_nel--; }
            della = xval - m_xx[1];
            nelm = double(m_wsum.as());
            m_wsum -= wt;
            nel = double(m_wsum.as());
            if (nel > 0) {
                delnel = della * double(wt) / nel;
                m_xx[1] -= delnel;
                drat = delnel * nel;
                nbyn = -double(wt) / nel;
                ac_dn = drat*drat;
                ac_on = nbyn;

                for (int ppp=2;ppp <= m_ord;ppp++) {
                    m_xx[ppp] -= ac_dn * (1.0 - ac_on);
                    if (ord_beyond) {
                        if (ppp < m_ord) {
                            ac_dn *= drat;
                            ac_on *= nbyn;
                        }
                        ac_de = -delnel;
                        for (int qqq=1;qqq <= ppp-2;qqq++) {
                            m_xx[ppp] -= bincoef[ppp][qqq] * ac_de * m_xx[ppp-qqq];
                            if (qqq < ppp - 2) {
                                ac_de *= -delnel;
                            }
                        }
                    }
                }
            }
            return *this;
        }
        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
            int ppp,qqq;
            n1 = double(m_wsum.as());
            if (n1 <= 0) {
                if (has_wts) { m_nel = rhs.m_nel; }
                m_wsum = rhs.m_wsum;
                m_xx = rhs.m_xx;
                return *this;
            }
            // else onboard the observations
            if (has_wts) { m_nel += rhs.m_nel; }
            n2 = double(rhs.m_wsum.as());
            m_wsum += rhs.m_wsum;
            if (n2 <= 0) {
                return *this;
            }
            ntot = double(m_wsum.as());
            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            del21 = rhs.m_xx[1] - m_xx[1];
            mupart = del21 * n2rat;

            m_xx[1] += mupart;
            nfoo = n1 * mupart;
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
            return *this;
        }
        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1;
            double ac_del,ac_mn2,ac_n1;
            int ppp,qqq;

            ntot = double(m_wsum.as());
            n2 = double(rhs.m_wsum.as());

            if (n2 <= 0) { return *this; }
            if (n2 > ntot) { stop("cannot subtract more observations than were seen."); }

            mupart = rhs.m_xx[1] - m_xx[1];

            m_wsum -= rhs.m_wsum;
            n1 = double(m_wsum.as());
            if (has_wts) { m_nel -= rhs.m_nel; }

            n1rat = n1 / ntot;
            n2rat = n2 / ntot;

            m_xx[1] -= (n2/n1) * mupart;

            del21 = mupart / n1rat;
            nfoo = mupart * n2;

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
            return *this;
        }

};
//UNFOLD
// no wts//FOLDUP
// ord_beyond must be used for (ord > 2)
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
template<class W,bool ord_beyond>
class Welford<W,false,ord_beyond> {
    public:
        int m_ord;
        NumericVector m_xx;
    public:
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_xx(NumericVector(ord+1)) {}
        inline Welford(const int &ord, 
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_xx(NumericVector(xx)) {}
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_xx(NumericVector(xx)) {}
    private:
        int m_nel;
    public:
        // reset to zero
        inline Welford& tare() {
            m_nel = 0;
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const bool normalize,const double used_df) const {
            double renorm;
            return ((m_xx[2]) / (double(m_nel) - used_df));
        }
        inline double mean() const {
            return m_xx[1];
        }
        inline double sd(const bool normalize,const double used_df) const {
            return sqrt(var(normalize,used_df));
        }
        inline double skew() const {
            return (sqrt(double(m_nel)) * m_xx[3] / pow(m_xx[2],1.5));
        }
        inline double exkurt() const {
            return ((double(m_nel) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
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
        inline int nel() const { return m_nel; }
        inline W wsum() const { return W(m_nel); }
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
            resu[0] = double(m_nel);
            return resu;
        }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            della = xval - m_xx[1];
            nelm = double(m_nel);
            m_nel++;
            nel = double(m_nel);
            delnel = della * double(wt) / nel;
            m_xx[1] += delnel;
            if (nelm > 0) {
                drat = della * nelm / nel;
                nbyn = -wt / nelm;
                ac_dn = pow(drat,m_ord);
                ac_on = pow(nbyn,m_ord-1);

                for (int ppp=m_ord;ppp >= 2;ppp--) {
                    m_xx[ppp] += ac_dn * wt * (1.0 - ac_on);
                    if (ord_beyond) {
                        if (ppp > 2) {
                            if (drat != 0) { ac_dn /= drat; }
                            ac_on /= nbyn;
                            ac_de = -delnel;

                            for (int qqq=1;qqq <= ppp-2;qqq++) {
                                m_xx[ppp] += bincoef[ppp][qqq] * ac_de * m_xx[ppp-qqq];
                                if (qqq < ppp - 2) {
                                    ac_de *= -delnel;
                                }
                            }
                        }
                    }
                }
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            della = xval - m_xx[1];
            nelm = double(m_nel);
            --m_nel;
            nel = double(m_nel);
            if (nel > 0) {
                delnel = della * double(wt) / nel;
                m_xx[1] -= delnel;
                drat = delnel * nel;
                nbyn = -double(wt) / nel;
                ac_dn = drat*drat;
                ac_on = nbyn;

                for (int ppp=2;ppp <= m_ord;ppp++) {
                    m_xx[ppp] -= ac_dn * (1.0 - ac_on);
                    if (ord_beyond) {
                        if (ppp < m_ord) {
                            ac_dn *= drat;
                            ac_on *= nbyn;
                        }
                        ac_de = -delnel;
                        for (int qqq=1;qqq <= ppp-2;qqq++) {
                            m_xx[ppp] -= bincoef[ppp][qqq] * ac_de * m_xx[ppp-qqq];
                            if (qqq < ppp - 2) {
                                ac_de *= -delnel;
                            }
                        }
                    }
                }
            }
            return *this;
        }
        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
            int ppp,qqq;
            if (m_nel <= 0) {
                m_nel = rhs.m_nel;
                // clone it?
                m_xx = Rcpp::clone(rhs.m_xx);
                return *this;
            }
            // else onboard the observations
            if (rhs.m_nel <= 0) {
                return *this;
            }
            n1 = double(m_nel);
            n2 = double(rhs.m_nel);
            m_nel += rhs.m_nel; 
            ntot = double(m_nel);
            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            del21 = rhs.m_xx[1] - m_xx[1];
            mupart = del21 * n2rat;

            m_xx[1] += mupart;
            nfoo = n1 * mupart;
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
            return *this;
        }
        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1;
            double ac_del,ac_mn2,ac_n1;
            int ppp,qqq;

            ntot = double(m_nel);
            n2 = double(rhs.m_nel);

            if (n2 <= 0) { return *this; }
            if (n2 > ntot) { stop("cannot subtract more observations than were seen."); }

            mupart = rhs.m_xx[1] - m_xx[1];

            m_nel -= rhs.m_nel;
            n1 = double(m_nel);

            n1rat = n1 / ntot;
            n2rat = n2 / ntot;

            m_xx[1] -= (n2/n1) * mupart;

            del21 = mupart / n1rat;
            nfoo = mupart * n2;

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
            return *this;
        }

};
//UNFOLD
// ord_beyond is false//FOLDUP
template<class W,bool has_wts>
class Welford<W,has_wts,false> {
    public:
        int m_ord;
        NumericVector m_xx;
    public:
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_wsum(Kahan<W>(0)), m_xx(NumericVector(ord+1)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
        inline Welford(const int &ord, 
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_wsum(Kahan<W>(sumwt)), m_xx(NumericVector(xx)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_wsum(Kahan<W>(W(xx[0]))), m_xx(NumericVector(xx)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
    private:
        Kahan<W> m_wsum;
        int m_nel;
    public:
        // reset to zero
        inline Welford& tare() {
            if (has_wts) { m_nel = 0; }
            m_wsum = W(0);
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const bool normalize,const double used_df) const {
            double renorm;
            if (has_wts && normalize) {
                renorm = double(m_nel) / double(m_wsum.as());
                return ((renorm * m_xx[2]) / (double(m_nel) - used_df));
            } else {
                return ((m_xx[2]) / (double(m_wsum.as()) - used_df));
            }
        }
        inline double mean() const {
            return m_xx[1];
        }
        inline double sd(const bool normalize,const double used_df) const {
            return sqrt(var(normalize,used_df));
        }
        inline double skew() const {
            return (sqrt(double(m_wsum.as())) * m_xx[3] / pow(m_xx[2],1.5));
        }
        inline double exkurt() const {
            return ((double(m_wsum.as()) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
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
        inline int nel() const { if (has_wts) { return m_nel; } else { return wsum(); } }
        inline W wsum() const { return m_wsum.as(); }
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
            resu[0] = double(wsum());
            return resu;
        }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            if (has_wts) { m_nel++; }
            della = xval - m_xx[1];
            nelm = double(m_wsum.as());
            m_wsum += wt;
            nel = double(m_wsum.as());
            delnel = della * double(wt) / nel;
            m_xx[1] += delnel;
            if (nelm > 0) {
                drat = della * nelm / nel;
                ac_dn = drat * drat;
                ac_on = -wt / nelm;
                m_xx[2] += ac_dn * wt * (1.0 - ac_on);
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            if (has_wts) { m_nel--; }
            della = xval - m_xx[1];
            nelm = double(m_wsum.as());
            m_wsum -= wt;
            nel = double(m_wsum.as());
            if (nel > 0) {
                delnel = della * double(wt) / nel;
                m_xx[1] -= delnel;
                drat = delnel * nel;
                ac_dn = drat*drat;
                ac_on = -double(wt) / nel;
                m_xx[2] -= ac_dn * (1.0 - ac_on);
            }
            return *this;
        }
        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
            int ppp,qqq;
            n1 = double(m_wsum.as());
            if (n1 <= 0) {
                if (has_wts) { m_nel = rhs.m_nel; }
                m_wsum = rhs.m_wsum;
                m_xx = rhs.m_xx;
                return *this;
            }
            // else onboard the observations
            if (has_wts) { m_nel += rhs.m_nel; }
            n2 = double(rhs.m_wsum.as());
            m_wsum += rhs.m_wsum;
            if (n2 <= 0) {
                return *this;
            }
            ntot = double(m_wsum.as());
            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            del21 = rhs.m_xx[1] - m_xx[1];
            mupart = del21 * n2rat;

            m_xx[1] += mupart;
            nfoo = n1 * mupart;
            ac_nfoo = nfoo*nfoo;
            ac_n2 = 1.0/n2;
            ac_mn1 = -1.0/n1;
            m_xx[2] += rhs.m_xx[2] + (ac_nfoo * (ac_n2 - ac_mn1));
            return *this;
        }
        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1;
            double ac_del,ac_mn2,ac_n1;
            int ppp,qqq;

            ntot = double(m_wsum.as());
            n2 = double(rhs.m_wsum.as());

            if (n2 <= 0) { return *this; }
            if (n2 > ntot) { stop("cannot subtract more observations than were seen."); }

            mupart = rhs.m_xx[1] - m_xx[1];

            m_wsum -= rhs.m_wsum;
            n1 = double(m_wsum.as());
            if (has_wts) { m_nel -= rhs.m_nel; }

            n1rat = n1 / ntot;
            n2rat = n2 / ntot;

            m_xx[1] -= (n2/n1) * mupart;

            del21 = mupart / n1rat;
            nfoo = mupart * n2;

            ac_nfoo = nfoo * nfoo;
            ac_n2 = 1.0 / n2;
            ac_mn1 = -1.0 / n1;
            m_xx[2] -= rhs.m_xx[2] + (ac_nfoo * (ac_n2 - ac_mn1));
            return *this;
        }

};
//UNFOLD
// no wts, no ord_beyond //FOLDUP
// ord_beyond must be used for (ord > 2)
// when has_wts is true, we accumulate the number of
// elements in m_nel; 
template<class W>
class Welford<W,false,false> {
    public:
        int m_ord;
        NumericVector m_xx;
    public:
        inline Welford(const int &ord) : m_ord(ord), m_nel(0), m_xx(NumericVector(ord+1)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
        inline Welford(const int &ord, 
                       const int &nel, 
                       const W &sumwt, 
                       const NumericVector &xx) : m_ord(ord), m_nel(nel), m_xx(NumericVector(xx)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
        inline Welford(const int &ord, 
                       const NumericVector &xx) : m_ord(ord), m_nel(int(xx[0])), m_xx(NumericVector(xx)) {
            if (ord != 2) { stop("should use ord >= 2"); }
        }
    private:
        int m_nel;
    public:
        // reset to zero
        inline Welford& tare() {
            m_nel = 0;
            for (int iii=0;iii < m_xx.length();++iii) {
                m_xx[iii] = 0;
            }
            return *this;
        }
        inline double var(const bool normalize,const double used_df) const {
            double renorm;
            return ((m_xx[2]) / (double(m_nel) - used_df));
        }
        inline double mean() const {
            return m_xx[1];
        }
        inline double sd(const bool normalize,const double used_df) const {
            return sqrt(var(normalize,used_df));
        }
        inline double skew() const {
            return (sqrt(double(m_nel)) * m_xx[3] / pow(m_xx[2],1.5));
        }
        inline double exkurt() const {
            return ((double(m_nel) * m_xx[4] / (pow(m_xx[2],2.0))) - 3.0);
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
        inline int nel() const { return m_nel; }
        inline W wsum() const { return W(m_nel); }
        inline NumericVector as() const { return m_xx; }
        inline NumericVector asvec() const { 
            // copy
            //NumericVector resu = NumericVector(m_xx);
            NumericVector resu = Rcpp::clone(m_xx);
            resu[0] = double(m_nel);
            return resu;
        }
    public:
        // add another (weighted) observation to our set of x
        inline Welford& add_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            della = xval - m_xx[1];
            nelm = double(m_nel);
            m_nel++;
            nel = double(m_nel);
            delnel = della * double(wt) / nel;
            m_xx[1] += delnel;
            if (nelm > 0) {
                drat = della * nelm / nel;
                ac_dn = drat * drat;
                ac_on = -wt / nelm;
                m_xx[2] += ac_dn * wt * (1.0 - ac_on);
            }
            return *this;
        }
        // remove one (weighted) observation from our set of x
        inline Welford& rem_one (const double xval, const W wt) {
            double della,nel,delnel,nelm,drat,nbyn,ac_dn,ac_on,ac_de;
            della = xval - m_xx[1];
            nelm = double(m_nel);
            --m_nel;
            nel = double(m_nel);
            if (nel > 0) {
                delnel = della * double(wt) / nel;
                m_xx[1] -= delnel;
                drat = delnel * nel;
                ac_dn = drat*drat;
                ac_on = -double(wt) / nel;
                m_xx[2] -= ac_dn * (1.0 - ac_on);
            }
            return *this;
        }
        // join two Welford objects together
        inline Welford& join(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1,ac_del,ac_mn2,ac_n1;
            int ppp,qqq;
            if (m_nel <= 0) {
                m_nel = rhs.m_nel;
                // clone it?
                m_xx = Rcpp::clone(rhs.m_xx);
                return *this;
            }
            // else onboard the observations
            if (rhs.m_nel <= 0) {
                return *this;
            }
            n1 = double(m_nel);
            n2 = double(rhs.m_nel);
            m_nel += rhs.m_nel; 
            ntot = double(m_nel);
            n1rat = n1 / ntot;
            n2rat = n2 / ntot;
            del21 = rhs.m_xx[1] - m_xx[1];
            mupart = del21 * n2rat;

            m_xx[1] += mupart;
            nfoo = n1 * mupart;
            ac_nfoo = nfoo*nfoo;
            ac_n2 = 1.0 / n2;
            ac_mn1 = -1.0 / n1;
            m_xx[2] += rhs.m_xx[2] + (ac_nfoo * (ac_n2 - ac_mn1));
            return *this;
        }
        // remove one from another
        inline Welford& unjoin(const Welford& rhs) {
            double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
            double ac_nfoo,ac_n2,ac_mn1;
            double ac_del,ac_mn2,ac_n1;
            int ppp,qqq;

            ntot = double(m_nel);
            n2 = double(rhs.m_nel);

            if (n2 <= 0) { return *this; }
            if (n2 > ntot) { stop("cannot subtract more observations than were seen."); }

            mupart = rhs.m_xx[1] - m_xx[1];

            m_nel -= rhs.m_nel;
            n1 = double(m_nel);

            n1rat = n1 / ntot;
            n2rat = n2 / ntot;

            m_xx[1] -= (n2/n1) * mupart;

            del21 = mupart / n1rat;
            nfoo = mupart * n2;

            ac_nfoo = nfoo * nfoo;
            ac_n2 = 1.0 / n2;
            ac_mn1 = -1.0 / n1;
            m_xx[2] -= rhs.m_xx[2] + (ac_nfoo * (ac_n2 - ac_mn1));
            return *this;
        }

};
//UNFOLD

// univariate sums, moments, cumulants//FOLDUP

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
// have mean 1.

template <typename T,typename W,typename oneW,bool has_wts,bool ord_beyond>
Welford<oneW,has_wts,ord_beyond> quasiWeightedThing(T v,
                                                    W wts,
                                                    int ord,
                                                    int bottom,
                                                    int top,
                                                    const bool na_rm,
                                                    const bool check_wts) {
    double nextv, nextw;
    Welford<oneW,has_wts,ord_beyond> frets = Welford<oneW,has_wts,ord_beyond>(ord);

    if ((top < 0) || (top > v.size())) { top = v.size(); }
    if (has_wts) {
        if (wts.size() < top) { stop("size of wts does not match v"); }
        for (int iii=bottom;iii < top;++iii) {
            nextv = v[iii];
            nextw = double(wts[iii]); 
            if (check_wts && (nextw < 0)) { stop("negative weight detected"); }

            if (! (na_rm && (ISNAN(nextv) || ISNAN(nextw)))) {
                // 2FIX: check for zero weight??
                frets.add_one(nextv,nextw);
            }
        }
    } else {
        for (int iii=bottom;iii < top;++iii) {
            nextv = v[iii];
            if (! (na_rm && ISNAN(nextv))) { frets.add_one(nextv,1); }
        }
    }
    return frets;
}

template <typename T,typename W,typename oneW,bool has_wts>
NumericVector quasiWeightedMoments(T v,
                                   W wts,
                                   int ord,
                                   int bottom,
                                   int top,
                                   const bool na_rm,
                                   const bool check_wts,
                                   const bool normalize_wts) {
    double nextv, nextw, renorm, nok;
    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }
    NumericVector xret;

    if (ord > 2) {
        Welford<oneW,has_wts,true> frets = quasiWeightedThing<T,W,oneW,has_wts,true>(v,wts,ord,bottom,top,na_rm,check_wts);
        xret = frets.asvec();
        nok = double(frets.nel());
    } else {
        Welford<oneW,has_wts,false> irets = quasiWeightedThing<T,W,oneW,has_wts,false>(v,wts,ord,bottom,top,na_rm,check_wts);
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
// fix ord_beyond
template <typename T>
NumericVector quasiWeightedMomentsCurryOne(T v, 
                                           SEXP wts, 
                                           int ord, 
                                           const bool na_rm, 
                                           const bool check_wts, 
                                           const bool normalize_wts) {
    NumericVector dummy_wts;
    if (!Rf_isNull(wts)) {  
        switch (TYPEOF(wts)) {
            case  INTSXP: { return quasiWeightedMoments<T,IntegerVector,int,true>(v, wts, ord, 0, -1, na_rm, check_wts, normalize_wts); }
            case REALSXP: { return quasiWeightedMoments<T,NumericVector,double,true>(v, wts, ord, 0, -1, na_rm, check_wts, normalize_wts); }
            case  LGLSXP: { return quasiWeightedMoments<T,LogicalVector,int,true>(v, wts, ord, 0, -1, na_rm, check_wts, normalize_wts); }
            default: stop("Unsupported weight type"); // nocov
        }
    }
    return quasiWeightedMoments<T,NumericVector,int,false>(v, dummy_wts, ord, 0, -1, na_rm, check_wts, normalize_wts); 
}

// wrap one level
NumericVector quasiWeightedMomentsCurryTwo(SEXP v, 
                                           SEXP wts, 
                                           int ord, 
                                           const bool na_rm, 
                                           const bool check_wts, 
                                           const bool normalize_wts) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return quasiWeightedMomentsCurryOne<IntegerVector>(v, wts, ord, na_rm, check_wts, normalize_wts); }
        case REALSXP: { return quasiWeightedMomentsCurryOne<NumericVector>(v, wts, ord, na_rm, check_wts, normalize_wts); }
        case  LGLSXP: { return quasiWeightedMomentsCurryOne<LogicalVector>(v, wts, ord, na_rm, check_wts, normalize_wts); }
        default: stop("Unsupported weight type"); // nocov
    }
    // have to have fallthrough for CRAN check.
    NumericVector retv;
    return retv;
}


// compute the standard deviation and mean and # df
//' @title
//' Compute first K moments
//' @description
//' Compute the (standardized) 2nd through kth moments, the mean, and the number of elements.
//' 
//' @param v a vector
//' @param na_rm whether to remove NA, false by default.
//' @param max_order the maximum order of the centered moment to be computed.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the centered moments computation. These are subtracted from the number of
//' observations. 
//' @param sg_df the number of degrees of freedom consumed in the computation of
//' the variance or standard deviation. This defaults to 1 to match the 
//' \sQuote{Bessel correction}.
//'
//' @details
//'
//' Computes the number of elements, the mean, and the 2nd through kth
//' centered standardized moment, for \eqn{k=2,3,4}{k=2,3,4}. These
//' are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//' In general they will \emph{not} match exactly with the 'standard'
//' implementations, due to differences in roundoff.
//'
//' These methods are reasonably fast, on par with the 'standard' implementations.
//' However, they will usually be faster than calling the various standard implementations
//' more than once.
//'
//' @return a vector, filled out as follows:
//' \describe{
//' \item{sd3}{A vector of the (sample) standard devation, mean, and number of elements (or the total weight when \code{wts}
//' are given).}
//' \item{skew4}{A vector of the (sample) skewness, standard devation, mean, and number of elements (or the total weight when 
//' \code{wts} are given).}
//' \item{kurt5}{A vector of the (sample) excess kurtosis, skewness, standard devation, mean, and number of elements (or the
//' total weight when \code{wts} are given).}
//' \item{cent_moments}{A vector of the (sample) \eqn{k}th centered moment, then \eqn{k-1}th centered moment, ..., 
//'  then the \emph{variance}, the mean, and number of elements (total weight when \code{wts} are given).}
//' \item{std_moments}{A vector of the (sample) \eqn{k}th standardized (and centered) moment, then 
//'  \eqn{k-1}th, ..., then standard devation, mean, and number of elements (total weight).}
//' \item{cent_cumulants}{A vector of the (sample) \eqn{k}th (centered, but this is redundant) cumulant, then the \eqn{k-1}th, ...,
//'  then the \emph{variance} (which is the second cumulant), then \emph{the mean}, then the number of elements (total weight).}
//' \item{std_cumulants}{A vector of the (sample) \eqn{k}th standardized (and centered, but this is redundant) cumulant, then the \eqn{k-1}th, ...,
//'  down to the third, then \emph{the variance}, \emph{the mean}, then the number of elements (total weight).}
//' }
//'
//' @note
//' The first centered (and standardized) moment is often defined to be identically 0. Instead \code{cent_moments}
//' and \code{std_moments} returns the mean. 
//' Similarly, the second standardized moments defined to be identically 1; \code{std_moments} instead returns the standard
//' deviation. The reason is that a user can always decide to ignore the results and fill in a 0 or 1 as they need, but 
//' could not efficiently compute the mean and standard deviation from scratch if we discard it.
//' The antepenultimate element of the output of \code{std_cumulants} is not a one, even though that \sQuote{should} be
//' the standardized second cumulant.
//' 
//' @note
//' The antepenultimate element of the output of \code{cent_moments}, \code{cent_cumulants} and \code{std_cumulants} is the \emph{variance},
//' not the standard deviation. All other code return the standard deviation in that place.
//'
//' @note
//' The kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @note
//' 'centered cumulants' is redundant. The intent was to avoid possible collision with existing code named 'cumulants'.
//'
//' @examples
//' x <- rnorm(1e5)
//' sd3(x)[1] - sd(x)
//' skew4(x)[4] - length(x)
//' skew4(x)[3] - mean(x)
//' skew4(x)[2] - sd(x)
//' if (require(moments)) {
//'   skew4(x)[1] - skewness(x)
//' }
//' # check 'robustness'; only the mean should change:
//' kurt5(x + 1e12) - kurt5(x)
//' # check speed
//' if (require(microbenchmark) && require(moments)) {
//'   dumbk <- function(x) { c(kurtosis(x) - 3.0,skewness(x),sd(x),mean(x),length(x)) }
//'   set.seed(1234)
//'   x <- rnorm(1e6)
//'   print(kurt5(x) - dumbk(x))
//'   microbenchmark(dumbk(x),kurt5(x),times=10L)
//' }
//' y <- std_moments(x,6)
//' cml <- cent_cumulants(x,6)
//' std <- std_cumulants(x,6)
//'
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector sd3(SEXP v, bool na_rm=false, SEXP wts = R_NilValue, double sg_df=1.0, bool check_wts=false, bool normalize_wts=true) {
    NumericVector preval = quasiWeightedMomentsCurryTwo(v, wts, 2, na_rm, check_wts, normalize_wts);
    NumericVector vret = NumericVector::create(COMP_SD_TWO(preval,sg_df),
                                               preval[1],
                                               preval[0]);

    return vret;
}

// return the skew, the standard deviation, the mean, and the dof
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector skew4(SEXP v, bool na_rm=false, SEXP wts = R_NilValue, double sg_df=1.0, bool check_wts=false, bool normalize_wts=true) {
    NumericVector preval = quasiWeightedMomentsCurryTwo(v, wts, 3, na_rm, check_wts, normalize_wts);
    NumericVector vret = NumericVector::create(COMP_SKEW(preval),
                                               COMP_SD_TWO(preval,sg_df),
                                               preval[1],
                                               preval[0]);
    return vret;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector kurt5(SEXP v, bool na_rm=false, SEXP wts = R_NilValue, double sg_df=1.0, bool check_wts=false, bool normalize_wts=true) {
    NumericVector preval = quasiWeightedMomentsCurryTwo(v, wts, 4, na_rm, check_wts, normalize_wts);
    NumericVector vret = NumericVector::create(COMP_EXKURT(preval),
                                               COMP_SKEW(preval),
                                               COMP_SD_TWO(preval,sg_df),
                                               preval[1],
                                               preval[0]);
    return vret;
}
// 2 functions useful for converting between forms
NumericVector sums2revm(NumericVector input,double used_df=0.0) {
    int ord = input.size() - 1;
    double nel = input[0] - used_df;
    int mmm;
    NumericVector output(ord+1);
    for (mmm=0;mmm <= MIN(1,ord);++mmm) {
        output[ord-mmm] = input[mmm];
    }
    for (mmm=2;mmm <= ord;++mmm) {
        output[ord-mmm] = input[mmm] / nel;
    }
    return output;
}
//NumericVector revm2sums(NumericVector input,double used_df=0.0) {
//    int ord = input.size() - 1;
//    double nel = input[ord] - used_df;
//    int mmm;
//    NumericVector output(ord+1);
//    for (mmm=0;mmm <= MIN(1,ord);++mmm) {
//        output[mmm] = input[ord-mmm];
//    }
//    for (mmm=2;mmm <= ord;++mmm) {
//        output[mmm] = input[ord-mmm] * nel;
//    }
//    return output;
//}
// return the centered moments up to order max_order
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector cent_moments(SEXP v, int max_order=5, int used_df=0, bool na_rm=false, SEXP wts=R_NilValue, 
                           bool check_wts=false, bool normalize_wts=true) {
    // 2FIX: add sg_df here?
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = quasiWeightedMomentsCurryTwo(v, wts, max_order, na_rm, check_wts, normalize_wts);
    NumericVector vret = sums2revm(preval,(double)used_df);
    return vret;
}
// return the standardized moments up to order max_order
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector std_moments(SEXP v, int max_order=5, int used_df=0, bool na_rm=false, SEXP wts=R_NilValue, 
                          bool check_wts=false, bool normalize_wts=true) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector cmoms = cent_moments(v, max_order, used_df, na_rm, wts, check_wts, normalize_wts);
    double sigma, adj;
    int mmm;
    if (max_order > 1) {
        adj = cmoms(max_order - 2);
        sigma = sqrt(adj);
        cmoms(max_order-2) = sigma;  // put the stdev in
        for (mmm=3;mmm <= max_order;++mmm) {
            adj *= sigma;
            cmoms(max_order-mmm) /= adj;
        }
    }
    return cmoms;
}
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector cent_cumulants(SEXP v, int max_order=5, int used_df=0, bool na_rm=false, SEXP wts=R_NilValue, 
                             bool check_wts=false, bool normalize_wts=true) {
    NumericVector cmoms = cent_moments(v,max_order,used_df,na_rm, wts, check_wts, normalize_wts);
    NumericVector cumuls(cmoms.size());
    int jjj,mmm;
    // copy over
    for (jjj=0;jjj < cumuls.size();++jjj) {
        cumuls(jjj) = cmoms(jjj);
    }
    if (max_order > 0) {
        // make it really centered? 
        cmoms(max_order - 1) = 0;
    }
    // moments to cumuls. it's a snap! (c.f. PDQutils)
    for (jjj=4;jjj <= max_order;++jjj) {
        // compute the jth order cumulant.
        for (mmm=2;mmm <= jjj-2;mmm++) {
            cumuls(max_order-jjj) -= bincoef[jjj-1][mmm-1] * cumuls(max_order-mmm) * cmoms(max_order-(jjj-mmm));
        }
    }
    return cumuls;
}
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector std_cumulants(SEXP v, int max_order=5, int used_df=0, bool na_rm=false, SEXP wts=R_NilValue, 
                            bool check_wts=false, bool normalize_wts=true) {
    NumericVector cumuls = cent_cumulants(v,max_order,used_df,na_rm,wts,check_wts,normalize_wts);
    double sigma,adj;
    int jjj;
    if (max_order > 1) {
        adj = cumuls(max_order-2);
        sigma = sqrt(adj);
        for (jjj=3;jjj <= max_order;++jjj) {
            adj *= sigma;
            cumuls(max_order - jjj) /= adj;
        }
    }
    return cumuls;
}
//UNFOLD

// monoid mumbo jumbo: add and subtract centered sums//FOLDUP

//' @title
//' Centered sums; join and unjoined.
//'
//' @description
//'
//' Compute, join, or unjoin centered sums.
//'
//' @param ret1 an \eqn{ord+1} vector as output by \code{\link{cent_sums}} consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//' @param ret2 an \eqn{ord+1} vector as output by \code{\link{cent_sums}} consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//' @param ret3 an \eqn{ord+1} vector as output by \code{\link{cent_sums}} consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//' @inheritParams cent_moments
//'
//' @return a vector the same size as the input consisting of the adjusted version of the input.
//' When there are not sufficient (non-nan) elements for the computation, \code{NaN} are returned.
//'
//' @examples
//'
//'  set.seed(1234)
//'  x1 <- rnorm(1e3,mean=1)
//'  x2 <- rnorm(1e3,mean=1)
//'  max_ord <- 6L
//'  rs1 <- cent_sums(x1,max_ord)
//'  rs2 <- cent_sums(x2,max_ord)
//'  rs3 <- cent_sums(c(x1,x2),max_ord)
//'  rs3alt <- join_cent_sums(rs1,rs2)
//'  stopifnot(max(abs(rs3 - rs3alt)) < 1e-7)
//'  rs1alt <- unjoin_cent_sums(rs3,rs2)
//'  rs2alt <- unjoin_cent_sums(rs3,rs1)
//'  stopifnot(max(abs(rs1 - rs1alt)) < 1e-7)
//'  stopifnot(max(abs(rs2 - rs2alt)) < 1e-7)
//'
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector cent_sums(SEXP v, int max_order=5, bool na_rm=false, SEXP wts=R_NilValue, bool check_wts=false, bool normalize_wts=true) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = quasiWeightedMomentsCurryTwo(v, wts, max_order, na_rm, check_wts, normalize_wts);
    return preval;
}
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector join_cent_sums(NumericVector ret1,NumericVector ret2) {
    if (ret1.size() != ret2.size()) { stop("mismatch in sizes."); }
    int ord = ret1.size() - 1;
    NumericVector cret1 = Rcpp::clone(ret1);
    NumericVector cret2 = Rcpp::clone(ret2);

    // always safe to do ord_beyond=true, as it is just an optimization.
    Welford<double,true,true> frets1 = Welford<double,true,true>(ord,cret1);
    Welford<double,true,true> frets2 = Welford<double,true,true>(ord,cret2);
    frets1.join(frets2);
    return frets1.asvec();
}
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector unjoin_cent_sums(NumericVector ret3,NumericVector ret2) {
    // subtract ret2 from ret3
    if (ret3.size() != ret2.size()) { stop("mismatch in sizes."); }
    int ord = ret3.size() - 1;
    NumericVector cret3 = Rcpp::clone(ret3);
    NumericVector cret2 = Rcpp::clone(ret2);

    // always safe to do ord_beyond=true, as it is just an optimization.
    Welford<double,true,true> frets3 = Welford<double,true,true>(ord,cret3);
    Welford<double,true,true> frets2 = Welford<double,true,true>(ord,cret2);
    frets3.unjoin(frets2);
    return frets3.asvec();
}
//UNFOLD

// 2FIX: add multivariate *marginal* operations ignoring comovement

// multivariate operations taking comovement into account//FOLDUP

// this function takes a n x p matrix of observations, and returns a (p+1) x (p+1) matrix;
// the 1,1 element is the count.
// the 1,2:(p+1) subcolumn is the mean
// the 2:(p+1),2:(p+1) submatrix is the squared sum (the covariance up to scaling by the inverse count)
// if na_omit is true, ignore rows of the input with any NA/NAN element.
template <typename T>
NumericMatrix quasiTheta(T v,bool na_omit = false) {
    const int n=v.nrow();
    const int p=v.ncol();

    double  nelm, nel;
    int iii,jjj,nnn;
    NumericVector mu(p);
    NumericVector della(p);
    NumericVector delnel(p);
    bool isok;

    // preallocated with zeros:
    NumericMatrix xret(1+p,1+p);

    for (nnn=0;nnn<n;nnn++) {
        isok = true;
        for (iii=0;iii<p;iii++) {
            della(iii) = v(nnn,iii) - xret(iii+1,0);
            if (na_omit && ISNAN(v(nnn,iii))) {
                isok = false;
                break;
            }
        }
        if (isok) {
            nelm = xret(0,0);
            nel = ++xret(0,0);
            for (iii=0;iii<p;iii++) {
                xret(iii+1,0) += della(iii) / nel;
                delnel(iii) = della(iii) * (nelm/nel);
            }
            for (iii=0;iii<p;iii++) {
                for (jjj=iii;jjj<p;jjj++) {
                    xret(1+iii,1+jjj) += della(iii) * delnel(jjj);
                }
            }
        }
    }
    // fill out the upper triangle
    for (iii=0;iii<p;iii++) {
        xret(0,1+iii) = xret(1+iii,0);
        for (jjj=iii+1;jjj<p;jjj++) {
            xret(1+jjj,1+iii) = xret(1+iii,1+jjj);
        }
    }

    return xret;
}


//' @title
//' Multivariate centered sums; join and unjoined.
//'
//' @description
//'
//' Compute, join, or unjoin multivariate centered (co-) sums.
//'
//' @param v an \eqn{m} by \eqn{n} matrix, each row an independent observation of some
//' \eqn{n} variate variable.
//' @param max_order the maximum order of cosum to compute. For now this can only be
//' 2; in the future higher order cosums should be possible.
//' @param na_omit a boolean; if \code{TRUE}, then only rows of \code{v} with complete
//' observations will be used.
//' @param ret1 a multdimensional array as output by \code{\link{cent_cosums}}.
//' @param ret2 a multdimensional array as output by \code{\link{cent_cosums}}.
//' @param ret3 a multdimensional array as output by \code{\link{cent_cosums}}.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the centered moments computation. These are subtracted from the number of
//' observations.
//'
//' @return a multidimensional arry of dimension \code{max_order}, each side of length
//' \eqn{1+n}. For the case currently implemented where \code{max_order} must be 2, the
//' output is a symmetric matrix, where the element in the \code{1,1} position is the count of 
//' complete) rows of \code{v}, the \code{2:(n+1),1} column is the mean, and the
//' \code{2:(n+1),2:(n+1)} is the co \emph{sums} matrix, which is the covariance up to scaling
//' by the count. \code{cent_comoments} performs this normalization for you.
//'
//' @seealso cent_sums
//'
//' @examples
//'
//'  set.seed(1234)
//'  x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
//'  x2 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
//'  max_ord <- 2L
//'  rs1 <- cent_cosums(x1,max_ord)
//'  rs2 <- cent_cosums(x2,max_ord)
//'  rs3 <- cent_cosums(rbind(x1,x2),max_ord)
//'  rs3alt <- join_cent_cosums(rs1,rs2)
//'  stopifnot(max(abs(rs3 - rs3alt)) < 1e-7)
//'  rs1alt <- unjoin_cent_cosums(rs3,rs2)
//'  rs2alt <- unjoin_cent_cosums(rs3,rs1)
//'  stopifnot(max(abs(rs1 - rs1alt)) < 1e-7)
//'  stopifnot(max(abs(rs2 - rs2alt)) < 1e-7)
//'
//' @template etc
//' @template ref-romo
//' @rdname centcosums 
//' @export
// [[Rcpp::export]]
NumericMatrix cent_cosums(SEXP v, int max_order=2, bool na_omit=false) {
    if (max_order != 2) { stop("only support order 2 for now"); }
    NumericMatrix retv;
    switch (TYPEOF(v)) {
        case  INTSXP: { retv = quasiTheta<IntegerMatrix>(v, na_omit); break; }
        case REALSXP: { retv = quasiTheta<NumericMatrix>(v, na_omit); break; }
        case  LGLSXP: { retv = quasiTheta<LogicalMatrix>(v, na_omit); break; }
        default: stop("Unsupported input type");
    }
    return retv;
}
//' @rdname centcosums 
//' @export
// [[Rcpp::export]]
NumericMatrix cent_comoments(SEXP v, int max_order=2, int used_df=0, bool na_omit=false) {
    NumericMatrix retv = cent_cosums(v,max_order,na_omit);
    double denom = retv(0,0) - used_df;
    int osize = retv.ncol();
    int iii,jjj;
    for (iii=1;iii<osize;++iii) {
        for (jjj=1;jjj<osize;++jjj) {
            retv(iii,jjj) /= denom;
        }
    }
    return retv;
}
//' @rdname centcosums 
//' @export
// [[Rcpp::export]]
NumericMatrix join_cent_cosums(NumericMatrix ret1,NumericMatrix ret2) {
    if ((ret1.ncol() != ret1.nrow()) ||
        (ret2.ncol() != ret2.nrow())) {
        stop("malformed input"); // nocov
    }

    const int p=ret1.ncol() - 1;
    double n1,n2,ntot,n2rat,muv;
    int iii,jjj;
    
    NumericVector della(p);
    NumericVector delnel(p);
    NumericMatrix ret3(p+1,p+1);
    
    n1 = ret1(0,0);
    if (n1 <= 0) {
        return ret2;
    }
    n2 = ret2(0,0);
    if (n2 <= 0) {
        return ret1;
    }
    ntot = n1 + n2;
    ret3(0,0) = ntot;
    n2rat = n2 / ntot;

    for (iii=0;iii<p;iii++) {
        muv = ret1(iii+1,0);
        della(iii) = ret2(iii+1,0) - muv;
        delnel(iii) = della(iii) * n2rat;
        ret3(iii+1,0) = muv + delnel(iii);
    }
    for (iii=0;iii<p;iii++) {
        for (jjj=iii;jjj<p;jjj++) {
            ret3(iii+1,jjj+1) = ret1(iii+1,jjj+1) + ret2(iii+1,jjj+1) + n1 * delnel(iii) * della(jjj);
        }
    }
    // fill out the upper triangle
    for (iii=0;iii<p;iii++) {
        ret3(0,1+iii) = ret3(1+iii,0);
        for (jjj=iii+1;jjj<p;jjj++) {
            ret3(1+jjj,1+iii) = ret3(1+iii,1+jjj);
        }
    }
    return ret3;
}
//' @rdname centcosums 
//' @export
// [[Rcpp::export]]
NumericMatrix unjoin_cent_cosums(NumericMatrix ret3,NumericMatrix ret2) {
    // subtract ret2 from ret3
    if ((ret3.ncol() != ret3.nrow()) ||
        (ret2.ncol() != ret2.nrow())) {
        stop("malformed input"); // nocov
    }

    const int p=ret3.ncol() - 1;

    double n1,n2,ntot,n1rat,n2rat,muv;
    int iii,jjj;

    NumericVector della(p);
    NumericVector delnel(p);
    NumericMatrix ret1(p+1,p+1);
    
    ntot = ret3(0,0);
    n2 = ret2(0,0);
    n1 = ntot - n2;
    if (0 > n1) { stop("cannot subtract more observations than we have. Do you have the order of arguments right?"); }
    if (n1 == 0) { 
        // would be better to check they are all equal and throw if not...
        return ret1;
    }
    ret1(0,0) = n1;
    n1rat = n1 / ntot;
    n2rat = n2 / ntot;

    // start HERE

    for (iii=0;iii<p;iii++) {
        muv = ret3(iii+1,0);
        della(iii) = (ret2(iii+1,0) - muv) / n1rat;
        delnel(iii) = della(iii) * n2rat;
        ret1(iii+1,0) = muv - delnel(iii);
    }
    for (iii=0;iii<p;iii++) {
        for (jjj=iii;jjj<p;jjj++) {
            ret1(iii+1,jjj+1) = ret3(iii+1,jjj+1) - ret2(iii+1,jjj+1) - n1 * delnel(iii) * della(jjj);
        }
    }
    // fill out the upper triangle
    for (iii=0;iii<p;iii++) {
        ret1(0,1+iii) = ret1(1+iii,0);
        for (jjj=iii+1;jjj<p;jjj++) {
            ret1(1+jjj,1+iii) = ret1(1+iii,1+jjj);
        }
    }
    return ret1;
}
//UNFOLD

// running sums, moments, cumulants, approximate quantiles//FOLDUP

/*******************************************************************************
 * running moments
 */


// helper function; takes a double or integer windowsize and interprets w/out warning or vomit.
// if NULL, return NA_INTEGER;
// if integer, pass through as i ;
// if double, then 
//   if Inf or NA, return NA_INTEGER;
//   else convert to integer via as<int>( )
int get_wins(SEXP window) {
    if (!Rf_isNull(window)) {
        switch (TYPEOF(window)) {
            case  INTSXP: { return as<int>(window); }
            case REALSXP: { 
                              double wins = as<double>(window);
                              if ((NumericVector::is_na(wins)) || 
                                  (traits::is_infinite<REALSXP>(wins) && (wins > 0.0))) {
                                  return NA_INTEGER;
                              }
                              return (int)wins;
                          }
            default: stop("Unsupported input type");
        }
    }
    return NA_INTEGER; 
}

//template <T>
//Rcpp::Nullable<T> get_wts(SEXP wts) {
    //if (!Rf_isNull(window)) {
        //return Rcpp::Nullable<T>(wts);
    //}
    //return Rcpp::Nullable<T>();
//}

enum ReturnWhat { ret_centmaxonly, // maxonly is a *centered* moment
    ret_centmoments, 
    ret_stdmoments, 
    ret_sd3, ret_skew4, ret_exkurt5,
    centered, scaled, zscore, sharpe, tstat, sharpese, 
    ret_stdev, ret_skew, ret_exkurt, 
    ret_sum, ret_mean };

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

    jjj = 0;
    // now run through iii
    for (iii=0;iii < numel;++iii) {
        if (!do_recompute || (subcount < recom_period)) {
            // add one
            if (has_wts) { //FOLDUP
                nextw = wts[iii];
                if (check_wts && (nextw < 0)) { stop("negative weight detected"); }
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
                    if (check_wts && (prevw < 0)) { stop("negative weight detected"); }
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
                    if (check_wts && (nextw < 0)) { stop("negative weight detected"); }
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

// to keep it DRY, this function has a bunch of variants,
// depending on the templated bools in the vanilla form,
// this function returns a NumericMatrix with as many rows
// as elements in the input, and ord+1 columns.
// unlike the quasiWeightedMoments code, the *last* column
// is the number of elements,
// the last minus one is the mean and so on;
// this simplifies the transformation to moments later.
//
//   the first column is the number of elements, 
//   the second is the mean,
//   the (ord + 1 - k)th is the k-1th centered sum, defined as
//   as sum_j (v[j] - mean)^i
//
// moreover we adapt a sliding window of size window.
// for computational efficiency, we add and subtract
// observations. this can lead to roundoff issues,
// especially in the subtraction of observations.
// the algorithm checks for negative second moment and
// starts afresh when encountered. Also, the computation
// is periodically restarted.
//
// in other forms, depending on templated bools, this
// computes the centered input, the rescaled input, the z-scored input
// the running sharpe or running t-score, as matrices with a single column.
//
// we have a lookahead option for the centered, scaled, and Z-scored
// variants. Positive lookahead means take info from the future.
//
// there is also a 'minimum df' parameter. this is the minimum count
// required to return data for the ret_cent, _scald, _z, _sr, _srmer, and _t forms.
// the reasoning is that you might not want the z-score on fweer than 10
// observations. these do the right thing wrt NA, BTW. some moments
// come out as zero when computed on too few observations, and we blindly
// return Inf or NaN in that case. set the min_df to correct for this.
// srmer is 'Sharpe ratio and Mertens standard error'
//
// in summary:
// ret_mat return a rows x (1+ord) matrix of the running centered sums

template <typename T,ReturnWhat retwhat,typename W,typename oneW,bool has_wts,bool ord_beyond>
NumericMatrix runningQMoments(T v,
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

    Welford<oneW,has_wts,ord_beyond> frets = Welford<oneW,has_wts,ord_beyond>(ord);

    // 2FIX:
    double nextv, prevv;
    double nextw;
    const bool renormalize = has_wts && normalize_wts;
    double renorm,mydf,denom;

    if (has_wts) {
        if (wts.size() < v.size()) { stop("size of wts does not match v"); }
    }

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(window);
    if ((window < 1) && (!infwin)) { stop("must give positive window"); }

    if (min_df < 0) { stop("require positive min_df"); }
    if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

    if ((((retwhat==scaled) || 
          (retwhat==zscore) || 
          (retwhat==sharpe) || 
          (retwhat==tstat) || 
          (retwhat==ret_stdev) || 
          (retwhat==ret_sd3)) && (ord < 2)) ||
        (((retwhat==ret_skew) ||
          (retwhat==ret_skew4)) && (ord < 3)) ||
        (((retwhat==sharpese) || 
          (retwhat==ret_exkurt) ||
          (retwhat==ret_exkurt5)) && (ord < 4))) { 
        stop("bad code: order too small to support this computation"); 
    }
    // only for retwhat==sharpese, but cannot define outside its scope.
    // no bigs.
    double sigma,skew,exkurt,sr;

    int iii,jjj,lll,mmm,ppp,qqq,tr_iii,tr_jjj;
    int numel = v.size();
    const bool non_aligned = (lookahead != 0);
    // refers to the number of *subtractions* performed
    int subcount = 0;

    // preallocated with zeros; should
    // probably be NA?
    int ncols;
    if ((retwhat==ret_centmoments) ||
        (retwhat==ret_stdmoments) ||
        (retwhat==ret_sd3) || 
        (retwhat==ret_skew4) ||
        (retwhat==ret_exkurt5)) {
        ncols = 1+ord; 
    } else if ((retwhat==sharpese)) {
        ncols = 2; 
    } else {
        ncols = 1; 
    }
        
    NumericMatrix xret(numel,ncols);
    NumericVector vret;

    // as an invariant, we will start the computation
    // with vret, which is initialized as the summed
    // means on [jjj,iii]
    tr_iii = lookahead - 1;
    tr_jjj = lookahead - window;
    // sneakily set subcount large so we just recompute
    // at head of loop. sneaky.
    subcount = recom_period;

    // now run through lll index
    for (lll=0;lll < numel;++lll) {
        tr_iii++;
        // check subcount first and just recompute if needed.
        if (subcount >= recom_period) {
            // fix this
            iii = MIN(numel-1,tr_iii);
            jjj = MAX(0,tr_jjj+1);
            if (jjj <= iii) {
                frets = quasiWeightedThing<T,W,oneW,has_wts,true>(v,wts,ord,
                                                                  jjj,       //bottom
                                                                  iii+1,     //top
                                                                  na_rm, check_wts);
            }
            subcount = 0;
        } else {
            if ((tr_iii < numel) && (tr_iii >= 0)) {
                // add on nextv:
                nextv = double(v[tr_iii]);
                if (has_wts) { 
                    nextw = double(wts[tr_iii]); 
                    if (check_wts && (nextw < 0)) { stop("negative weight detected"); }
                    if (! (na_rm && (ISNAN(nextv) || ISNAN(nextw) || (nextw <= 0)))) {
                        frets.add_one(nextv,nextw);
                    }
                } else { 
                    if (! (na_rm && (ISNAN(nextv)))) {
                        frets.add_one(nextv,1);
                    }
                }
            }
            // remove prevv:
            if ((tr_jjj < numel) && (tr_jjj >= 0)) {
                prevv = double(v[tr_jjj]);
                if (has_wts) { 
                    nextw = double(wts[tr_jjj]); 
                    // unnecessary as we would have caught it prior?
                    //if (check_wts && (nextw < 0)) { stop("negative weight detected"); }
                    if (! (na_rm && (ISNAN(nextv) || ISNAN(nextw) || (nextw <= 0)))) {
                        frets.rem_one(prevv,nextw);
                        subcount++;
                    }
                } else { 
                    if (! (na_rm && (ISNAN(nextv)))) {
                        frets.rem_one(prevv,1);
                        subcount++;
                    }
                }
            }
        }
        tr_jjj++;

        // fill in the value in the output.//FOLDUP
        if ((retwhat==ret_centmoments) ||
            (retwhat==ret_stdmoments) ||
            (retwhat==ret_sd3) ||
            (retwhat==ret_skew4) ||
            (retwhat==ret_exkurt5)) {

            vret = frets.asvec();
            if (renormalize) { 
                mydf = double(frets.nel());
                renorm = mydf / vret[0];
                vret[0] = mydf;
                for (int ppp=2;ppp <= ord;ppp++) {
                    vret[ppp] *= renorm;
                }
            } else {
                mydf = vret[0];
            }
            if (retwhat==ret_centmoments) {//FOLDUP
                if (mydf >= min_df) {
                    denom = vret[0] - used_df;
                    xret(lll,ord) = vret[0];
                    xret(lll,ord-1) = vret[1];
                    // put them in backwards!
                    if (mydf >= ord) {
                        for (mmm=2;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm] / denom;
                        }
                    } else {
                        for (mmm=2;mmm <= mydf;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm] / denom;
                        }
                        for (mmm=int(ceil(mydf))+1;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = NAN;
                        }
                    }
                } else {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,mmm) = NAN;
                    }
                }
            } //UNFOLD
            else if (retwhat==ret_stdmoments) {//FOLDUP
                if (mydf >= min_df) {
                    denom = vret[0] - used_df;
                    sigma = COMP_SD_TWO(vret,used_df);
                    xret(lll,ord) = vret[0];
                    xret(lll,ord-1) = vret[1];
                    xret(lll,ord-2) = sigma;
                    // put them in backwards!
                    if (mydf >= ord) {
                        for (mmm=3;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm] / (denom * pow(sigma,mmm));
                        }
                    } else {
                        for (mmm=3;mmm <= mydf;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm] / (denom * pow(sigma,mmm));
                        }
                        for (mmm=int(ceil(mydf))+1;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = NAN;
                        }
                    }
                } else {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,mmm) = NAN;
                    }
                }
            } //UNFOLD
            else if (retwhat==ret_sd3) {//FOLDUP
                if (mydf >= min_df) {
                    // put them in backwards!
                    if (mydf >= ord) {
                        for (mmm=0;mmm < ord;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        xret(lll,0) = COMP_SD_TWO(vret,used_df);
                    } else {
                        for (mmm=0;mmm <= mydf;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        for (mmm=int(ceil(mydf))+1;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = NAN;
                        }
                    }
                } else {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,mmm) = NAN;
                    }
                }
            }//UNFOLD
            else if (retwhat==ret_skew4) {//FOLDUP
                if (mydf >= min_df) {
                    // put them in backwards!
                    if (mydf >= ord) {
                        for (mmm=0;mmm < (ord-1);++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        xret(lll,0) = COMP_SKEW(vret);
                        xret(lll,1) = COMP_SD_TWO(vret,used_df);
                    } else {
                        for (mmm=0;mmm <= mydf;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        for (mmm=int(ceil(mydf))+1;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = NAN;
                        }
                        if (mydf >= (ord-1)) {
                            xret(lll,1) = COMP_SD_TWO(vret,used_df);
                        }
                    }
                } else {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,mmm) = NAN;
                    }
                }
            }//UNFOLD
            else if (retwhat==ret_exkurt5) {//FOLDUP
                if (mydf >= min_df) {
                    // put them in backwards!
                    if (mydf >= ord) {
                        for (mmm=0;mmm < (ord-1);++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        xret(lll,0) = COMP_EXKURT(vret);
                        xret(lll,1) = COMP_SKEW(vret);
                        xret(lll,2) = COMP_SD_TWO(vret,used_df);
                    } else {
                        for (mmm=0;mmm <= mydf;++mmm) {
                            xret(lll,ord-mmm) = vret[mmm];
                        }
                        for (mmm=int(ceil(mydf))+1;mmm <= ord;++mmm) {
                            xret(lll,ord-mmm) = NAN;
                        }
                        if (mydf >= (ord-2)) {
                            xret(lll,2) = COMP_SD_TWO(vret,used_df);
                            if (mydf >= (ord-1)) {
                                xret(lll,1) = COMP_SKEW(vret);
                            }
                        }
                    }
                } else {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,mmm) = NAN;
                    }
                }
            }//UNFOLD
        } else if (retwhat==ret_centmaxonly) {
            vret = frets.asvec();
            if (renormalize) { 
                mydf = double(frets.nel());
                renorm = mydf / vret[0];
                vret[0] = mydf;
                vret[ord] *= renorm;
            } else {
                mydf = vret[0];
            }
            if ((mydf >= min_df) && (mydf >= ord)) {
                xret(lll,0) = vret[ord] / (vret[0] - used_df);
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
                if (retwhat==centered) { xret(lll,0) = frets.centered(double(v[lll])); }
                if (retwhat==scaled) { xret(lll,0) = frets.scaled(double(v[lll]),renormalize,used_df); }
                if (retwhat==zscore) { xret(lll,0) = frets.zscored(double(v[lll]),renormalize,used_df); }
                if (retwhat==tstat) {
                    if (renormalize) {
                        xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.nel()));
                    } else {
                        xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.wsum()));
                    }
                }
                if (retwhat==sharpe) { xret(lll,0) = frets.sharpe(renormalize,used_df); }
                if (retwhat==sharpese) {
                    skew = frets.skew();
                    exkurt = frets.exkurt();
                    sr = frets.sharpe(renormalize,used_df);
                    xret(lll,0) = sr;
                    if (renormalize) {
                        xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / double(frets.nel()));
                    } else {
                        xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / double(frets.wsum()));
                    }
                }
                if (retwhat==ret_stdev) { xret(lll,0) = frets.sd(renormalize,used_df); }
                if (retwhat==ret_skew) { xret(lll,0) = frets.skew(); }
                if (retwhat==ret_exkurt) { xret(lll,0) = frets.exkurt(); }
            } else {
                xret(lll,0) = NAN;
            }
        }
        //UNFOLD
    }

    return xret;
}

template <typename T,ReturnWhat retwhat>
NumericMatrix runningQMomentsCurryOne(T v, 
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

    if (wts.isNotNull()) {
        return runningQMoments<T,retwhat,NumericVector,double,true,true>(v, wts.get(), ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
    }
    NumericVector dummy_wts;
    return runningQMoments<T,retwhat,NumericVector,double,false,true>(v, dummy_wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); 
}

template <ReturnWhat retwhat>
NumericMatrix runningQMomentsCurryTwo(SEXP v, 
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
        case  INTSXP: { return runningQMomentsCurryOne<IntegerVector,retwhat>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); } 
        case REALSXP: { return runningQMomentsCurryOne<NumericVector,retwhat>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); } 
        case  LGLSXP: { return runningQMomentsCurryOne<LogicalVector,retwhat>(v, wts, ord, window, recom_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts); } 
        default: stop("Unsupported weight type"); // nocov
    }
    // have to have fallthrough for CRAN check.
    NumericMatrix retv;
    return retv;
}

//' @title
//' Compute first K moments over a sliding window
//' @description
//' Compute the (standardized) 2nd through kth moments, the mean, and the number of elements over
//' an infinite or finite sliding window, returning a matrix.
//' 
//' @param v a vector
//' @param window the window size. if given as finite integer or double, passed through.
//' If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, equivalent
//' to an infinite window size. If negative, an error will be thrown.
//' @param restart_period the recompute period. because subtraction of elements can cause
//' loss of precision, the computation of moments is restarted periodically based on 
//' this parameter. Larger values mean fewer restarts and faster, though less accurate
//' results. Note that the code checks for negative second and fourth moments and
//' recomputes when needed.
//' @param na_rm whether to remove NA, false by default.
//' @param max_order the maximum order of the centered moment to be computed.
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent moments from being computed on too few observations.
//' Defaults to zero, meaning no restriction.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the centered moments computation. These are subtracted from the number of
//' observations.
//'
//' @details
//'
//' Computes the number of elements, the mean, and the 2nd through kth
//' centered (and typically standardized) moments, for \eqn{k=2,3,4}{k=2,3,4}. These
//' are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//'
//' Given the length \eqn{n} vector \eqn{x}, we output matrix \eqn{M} where
//' \eqn{M_{i,j}}{M_i,j} is the \eqn{order - j + 1} moment (\emph{i.e.}
//' excess kurtosis, skewness, standard deviation, mean or number of elements)
//' of \eqn{x_{i-window+1},x_{i-window+2},...,x_{i}}{x_(i-window+1),x_(i-window+2),...,x_i}.
//' Barring \code{NA} or \code{NaN}, this is over a window of size \code{window}.
//' During the 'burn-in' phase, we take fewer elements.
//'
//' @return Typically a matrix, where the first columns are the kth, k-1th through 2nd standardized, 
//' centered moments, then a column of the mean, then a column of the number of (non-nan) elements in the input,
//' with the following exceptions:
//' \describe{
//' \item{running_cent_moments}{Computes arbitrary order centered moments. When \code{max_order_only} is set,
//' only a column of the maximum order centered moment is returned.}
//' \item{running_std_moments}{Computes arbitrary order standardized moments, then the standard deviation, the mean,
//' and the count. There is not yet an option for \code{max_order_only}, but probably should be.}
//' \item{running_cumulants}{Computes arbitrary order cumulants, and returns the kth, k-1th, through the second 
//' (which is the variance) cumulant, then the mean, and the count.}
//' }
//'
//' @note
//' the kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @examples
//' x <- rnorm(1e5)
//' xs3 <- running_sd3(x,10)
//' xs4 <- running_skew4(x,10)
//'
//' if (require(moments)) {
//'     set.seed(123)
//'     x <- rnorm(5e1)
//'     window <- 10L
//'     kt5 <- running_kurt5(x,window=window)
//'     rm1 <- t(sapply(seq_len(length(x)),function(iii) { 
//'                 xrang <- x[max(1,iii-window+1):iii]
//'                 c(moments::kurtosis(xrang)-3.0,moments::skewness(xrang),
//'                 sd(xrang),mean(xrang),length(xrang)) },
//'              simplify=TRUE))
//'     stopifnot(max(abs(kt5 - rm1),na.rm=TRUE) < 1e-12)
//' }
//'
//' xc6 <- running_cent_moments(x,window=100L,max_order=6L)
//'
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_sd3(SEXP v, SEXP window = R_NilValue, 
                          Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                          bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                          bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    NumericMatrix preval = runningQMomentsCurryTwo<ret_sd3>(v, wts, 2, wins, restart_period, 0, min_df, used_df, 
                                                            na_rm, check_wts, normalize_wts);
    return preval;
}
// return the skew, the standard deviation, the mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_skew4(SEXP v, SEXP window = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                            bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                            bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    NumericMatrix preval = runningQMomentsCurryTwo<ret_skew4>(v, wts, 3, wins, restart_period, 0, min_df, used_df, 
                                                              na_rm, check_wts, normalize_wts);
    return preval;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_kurt5(SEXP v, SEXP window = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                            bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                            bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    NumericMatrix preval = runningQMomentsCurryTwo<ret_exkurt5>(v, wts, 4, wins, restart_period, 0, min_df, used_df, 
                                                                na_rm, check_wts, normalize_wts);
    return preval;
}
// just the sd nothing else.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_sd(SEXP v, SEXP window = R_NilValue, 
                         Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                         bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                         bool check_wts=false, bool normalize_wts=true) {
//2FIX: introduce used_df ... 
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<ret_stdev>(v, wts, 2, wins, restart_period, 0, min_df, used_df,
                                              na_rm, check_wts, normalize_wts);
}
// just the skew nothing else.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_skew(SEXP v, SEXP window = R_NilValue, 
                           Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                           bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                           bool check_wts=false, bool normalize_wts=true) {
//2FIX: introduce used_df ... 
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<ret_skew>(v, wts, 3, wins, restart_period, 0, min_df, used_df, 
                                             na_rm, check_wts, normalize_wts);
}
// just the kurtosis nothing else.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_kurt(SEXP v, SEXP window = R_NilValue, 
                           Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                           bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                           bool check_wts=false, bool normalize_wts=true) {
//2FIX: introduce used_df ... 
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<ret_exkurt>(v, wts, 4, wins, restart_period, 0, min_df, used_df,
                                               na_rm, check_wts, normalize_wts);
}

// return the centered moments down to the 2nd, then the mean, and the dof.
//' @param max_order_only for \code{running_cent_moments}, if this flag is set, only compute
//' the maximum order centered moment, and return in a vector.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_cent_moments(SEXP v, SEXP window = R_NilValue, 
                                   Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                   int max_order=5, bool na_rm=false, bool max_order_only=false, 
                                   int min_df=0, double used_df=0.0, int restart_period=100, 
                                   bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    if (max_order_only) {
        return runningQMomentsCurryTwo<ret_centmaxonly>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
                                                    na_rm, check_wts, normalize_wts);
    } 
    return runningQMomentsCurryTwo<ret_centmoments>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
                                                    na_rm, check_wts, normalize_wts);
}

// return the standardized moments down to the 3rd, then the standard deviation, the mean, and the dof.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_std_moments(SEXP v, SEXP window = R_NilValue, 
                                  Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                  int max_order=5, bool na_rm=false, 
                                  int min_df=0, double used_df=0, int restart_period=100, 
                                  bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<ret_stdmoments>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
                                                   na_rm, check_wts, normalize_wts);
}

// return the cumulants down to the 2nd, then the standard deviation, the mean, and the dof.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_cumulants(SEXP v, SEXP window = R_NilValue, 
                                Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                int max_order=5, bool na_rm=false, int min_df=0, double used_df=0.0, int restart_period=100,
                                bool check_wts=false, bool normalize_wts=true) {
    NumericMatrix cumulants = running_cent_moments(v, window, wts, max_order, na_rm, 
                                                   false, min_df, used_df, restart_period, check_wts, normalize_wts);

    NumericVector temp_moments(1+max_order);
    int iii,jjj,mmm,ppp;
    // moments to cumulants. it's a snap! (c.f. PDQutils)
    for (iii=0;iii < cumulants.nrow();++iii) {
        // copy the row to avoid writeover; bleah;
        for (jjj=0;jjj <= max_order;++jjj) {
            temp_moments(jjj) = cumulants(iii,jjj);
        }
        for (jjj=4;jjj <= max_order;++jjj) {
            // compute the jth order cumulant.
            for (mmm=2;mmm <= jjj-2;mmm++) {
                cumulants(iii,max_order-jjj) -= bincoef[jjj-1][mmm-1] * cumulants(iii,max_order-mmm) * temp_moments(max_order-(jjj-mmm));
            }
        }
    }
    return cumulants;
}
//' @title
//' Compute approximate quantiles over a sliding window
//' @description
//' Computes cumulants up to some given order, then employs the Cornish-Fisher approximation
//' to compute approximate quantiles using a Gaussian basis.
//' 
//' @param p the probability points at which to compute the quantiles. Should be in the range (0,1).
//' @inheritParams running_cumulants
//'
//' @details
//'
//' Computes the cumulants, then approximates quantiles using AS269 of Lee & Lin.
//'
//' @references 
//'
//' Lee, Y-S., and Lin, T-K. "Algorithm AS269: High Order Cornish Fisher
//' Expansion." Appl. Stat. 41, no. 1 (1992): 233-240. 
//' \url{http://www.jstor.org/stable/2347649}
//'
//' Lee, Y-S., and Lin, T-K. "Correction to Algorithm AS269: High Order 
//' Cornish Fisher Expansion." Appl. Stat. 42, no. 1 (1993): 268-269. 
//' \url{http://www.jstor.org/stable/2347433}
//'
//' AS 269. \url{http://lib.stat.cmu.edu/apstat/269}
//'
//' Jaschke, Stefan R. "The Cornish-Fisher-expansion in the context of 
//' Delta-Gamma-normal approximations." No. 2001, 54. Discussion Papers, 
//' Interdisciplinary Research Project 373: Quantification and Simulation of 
//' Economic Processes, 2001. 
//' \url{http://www.jaschke-net.de/papers/CoFi.pdf}
//'
//' @return A matrix, with one row for each element of \code{x}, and one column for each element of \code{q}.
//'
//' @note
//' The current implementation is not as space-efficient as it could be, as it first computes
//' the cumulants for each row, then performs the Cornish-Fisher approximation on a row-by-row
//' basis. In the future, this computation may be moved earlier into the pipeline to be more
//' space efficient. File an issue if the memory footprint is an issue for you.
//'
//' @examples
//' x <- rnorm(1e5)
//' xq <- running_apx_quantiles(x,c(0.1,0.25,0.5,0.75,0.9))
//' xm <- running_apx_median(x)
//'
//' @seealso \code{\link{running_cumulants}}, \code{PDQutils::qapx_cf}, \code{PDQutils::AS269}.
//' @template etc
//' @template ref-romo
//' @template param-wts
//' @rdname runningquantiles
//' @export
// [[Rcpp::export]]
NumericMatrix running_apx_quantiles(SEXP v, NumericVector p, SEXP window = R_NilValue, 
                                    Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                    int max_order=5, bool na_rm=false, int min_df=0, double used_df=0.0, int restart_period=100,
                                    bool check_wts=false, bool normalize_wts=true) {
    NumericMatrix cumulants = running_cumulants(v, window, wts, max_order, na_rm, min_df, used_df, restart_period, check_wts, normalize_wts);

    int iii,jjj,mmm,nnn,qqq,lll;
    int ja,jb,jal,jbl;
    double fac,aa,bc;
    double mu,sigmasq;

    int nq=p.size();
    int nord=max_order-2;

    // yay for sugar! yay!
    NumericVector z = qnorm(p,0.0,1.0);

    // prealloc
    NumericMatrix retval(cumulants.nrow(),nq);
    // line 20
    NumericMatrix P(nq,3*bincoef[nord+1][2]);
    // line 30
    NumericMatrix D(nq,3*nord);
    NumericVector DEL(nq);

    NumericVector DD(nq);
    // standardized cumulants go here:
    NumericVector a(nord);

    // line 10
    // precompute the Hermite polynomials
    NumericMatrix H(nq,3*nord);
    for (qqq=0;qqq<nq;qqq++) {
        H(qqq,0) = -z(qqq);
        H(qqq,1) = z(qqq) * z(qqq) - 1.0;
        for (jjj=2;jjj < 3*nord;jjj++) {
            H(qqq,jjj) = - (z(qqq) * H(qqq,jjj-1) + jjj * H(qqq,jjj-2));
        }
    }

    // now begins the fun!
    for (iii=0;iii < cumulants.nrow();iii++) {
        // zero everything
        for (mmm=0;mmm < nq;mmm++) {
            DEL(mmm) = 0.0;
            for (nnn=0;nnn<P.ncol();nnn++) { P(mmm,nnn) = 0.0; }
            for (nnn=0;nnn<D.ncol();nnn++) { D(mmm,nnn) = 0.0; }
            // probably unecessary:
            DD(mmm) = 0.0;
        }
        mu = cumulants(iii,max_order-1);
        sigmasq = cumulants(iii,max_order-2);
        // change raw cumulants to standardized...
        for (jjj=0;jjj<nord;jjj++) {
            a(jjj) = pow(-1.0,(jjj+1)) * cumulants(iii,max_order-3-jjj) / 
                (pow(sigmasq,(jjj+3.0)/2.0) * (double)((jjj+2) * (jjj+3)));
        }
        for (qqq=0;qqq<nq;qqq++) {
            D(qqq,0) = -a(0) * H(qqq,1);
            DEL(qqq) = D(qqq,0);
            P(qqq,0) = D(qqq,0);
            P(qqq,2) = a(0);
        }
        ja = 0;
        fac = 1.0;

        for (jjj=2;jjj<=nord;++jjj) {//FOLDUP
            fac = fac * jjj;
            ja = ja + 3 * (jjj-1);
            jb = ja;
            bc = 1.0;
            for (mmm=1;mmm<jjj;mmm++) {
                for (qqq=0;qqq<nq;qqq++) {
                    DD(qqq) = bc * D(qqq,mmm-1);
                }
                aa = bc * a(mmm-1);
                jb -= 3 * (jjj - mmm);
                for (lll=1;lll<=3*(jjj-mmm);lll++) {
                    jbl = jb + lll;
                    jal = ja + lll;
                    for (qqq=0;qqq<nq;qqq++) {
                        P(qqq,jal) += DD(qqq) * P(qqq,jbl-1);
                        P(qqq,jal+mmm+1) += aa * P(qqq,jbl-1);
                    }
                }  // line 40
                bc *= (jjj - mmm) / mmm;
            } // line 50
            for (qqq=0;qqq<nq;qqq++) {
                P(qqq,ja + jjj + 1) += a(jjj-1);
                // calculate the adjustments
                D(qqq,jjj-1) = 0.0;
            }
            for (lll=2;lll <= 3*jjj;lll++) {
                for (qqq=0;qqq<nq;qqq++) {
                    D(qqq,jjj-1) -= P(qqq,ja+lll-1) * H(qqq,lll-2);
                }
            }
            // line 60
            for (qqq=0;qqq<nq;qqq++) {
                P(qqq,ja) = D(qqq,jjj-1);
                DEL(qqq) += (D(qqq,jjj-1) / fac);
            }
        } // line 70//UNFOLD

        // interpret as mean plus some sdevs://FOLDUP
        for (qqq=0;qqq<nq;qqq++) {
            retval(iii,qqq) = mu + sqrt(sigmasq) * (DEL(qqq) + z(qqq));
        }//UNFOLD
    }

    return retval;
}
//' @rdname runningquantiles
//' @export
// [[Rcpp::export]]
NumericMatrix running_apx_median(SEXP v, SEXP window = R_NilValue, 
                                 Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                                 int max_order=5, bool na_rm=false, int min_df=0, double used_df=0.0, int restart_period=100,
                                 bool check_wts=false, bool normalize_wts=true) {
    NumericVector p(1);
    p(0) = 0.5;
    NumericMatrix vret = running_apx_quantiles(v,p,window,wts,max_order,na_rm,min_df,used_df,restart_period,check_wts,normalize_wts);
    return vret;
}

//' @title
//' Compare data to moments computed over a sliding window.
//' @description
//' Computes moments over a sliding window, then adjusts the data accordingly, centering, or scaling,
//' or z-scoring, and so on.
//' 
//' @inheritParams running_cent_moments
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent \emph{e.g.} Z-scores from being computed on only 3
//' observations. Defaults to zero, meaning no restriction, which can result in 
//' infinite Z-scores during the burn-in period.
//' @param lookahead for some of the operations, the value is compared to 
//' mean and standard deviation possibly using 'future' or 'past' information
//' by means of a non-zero lookahead. Positive values mean data are taken from
//' the future.
//' @param compute_se for \code{running_sharpe}, return an extra column of the
//' standard error, as computed by Mertens' correction.
//'
//' @details
//'
//' Given the length \eqn{n} vector \eqn{x}, for
//' a given index \eqn{i}, define \eqn{x^{(i)}}{x^(i)}
//' as the vector of 
//' \eqn{x_{i-window+1},x_{i-window+2},...,x_{i}}{x_(i-window+1),x_(i-window+2),...,x_i},
//' where we do not run over the 'edge' of the vector. In code, this is essentially
//' \code{x[(max(1,i-window+1)):i]}. Then define \eqn{\mu_i}{mu_i}, \eqn{\sigma_i}{sigma_i}
//' and \eqn{n_i}{n_i} as, respectively, the sample mean, standard deviation and number of
//' non-NA elements in \eqn{x^{(i)}}{x^(i)}. 
//'
//' We compute output vector \eqn{m} the same size as \eqn{x}. 
//' For the 'centered' version of \eqn{x}, we have \eqn{m_i = x_i - \mu_i}{m_i = x_i - mu_i}.
//' For the 'scaled' version of \eqn{x}, we have \eqn{m_i = x_i / \sigma_i}{m_i = x_i / sigma_i}.
//' For the 'z-scored' version of \eqn{x}, we have \eqn{m_i = (x_i - \mu_i) / \sigma_i}{m_i = (x_i - mu_i) / sigma_i}.
//' For the 't-scored' version of \eqn{x}, we have \eqn{m_i = \sqrt{n_i} \mu_i / \sigma_i}{m_i = sqrt(n_i) mu_i / sigma_i}.
//'
//' We also allow a 'lookahead' for some of these operations.
//' If positive, the moments are computed using data from larger indices;
//' if negative, from smaller indices. Letting \eqn{j = i + lookahead}{j = i + lookahead}:
//' For the 'centered' version of \eqn{x}, we have \eqn{m_i = x_i - \mu_j}{m_i = x_i - mu_j}.
//' For the 'scaled' version of \eqn{x}, we have \eqn{m_i = x_i / \sigma_j}{m_i = x_i / sigma_j}.
//' For the 'z-scored' version of \eqn{x}, we have \eqn{m_i = (x_i - \mu_j) / \sigma_j}{m_i = (x_i - mu_j) / sigma_j}.
//'
//' @return a vector the same size as the input consisting of the adjusted version of the input.
//' When there are not sufficient (non-nan) elements for the computation, \code{NaN} are returned.
//'
//' @examples
//'
//' if (require(moments)) {
//'     set.seed(123)
//'     x <- rnorm(5e1)
//'     window <- 10L
//'     rm1 <- t(sapply(seq_len(length(x)),function(iii) { 
//'                   xrang <- x[max(1,iii-window+1):iii]
//'                   c(sd(xrang),mean(xrang),length(xrang)) },
//'                   simplify=TRUE))
//'     rcent <- running_centered(x,window=window)
//'     rscal <- running_scaled(x,window=window)
//'     rzsco <- running_zscored(x,window=window)
//'     rshrp <- running_sharpe(x,window=window)
//'     rtsco <- running_tstat(x,window=window)
//'     rsrse <- running_sharpe(x,window=window,compute_se=TRUE)
//'     stopifnot(max(abs(rcent - (x - rm1[,2])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rscal - (x / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rzsco - ((x - rm1[,2]) / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rshrp - (rm1[,2] / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rtsco - ((sqrt(rm1[,3]) * rm1[,2]) / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rsrse[,1] - rshrp),na.rm=TRUE) < 1e-12)
//'
//'     rm2 <- t(sapply(seq_len(length(x)),function(iii) { 
//'                   xrang <- x[max(1,iii-window+1):iii]
//'                   c(kurtosis(xrang)-3.0,skewness(xrang)) },
//'                   simplify=TRUE))
//'     mertens_se <- sqrt((1 + ((2 + rm2[,1])/4) * rshrp^2 - rm2[,2]*rshrp) / rm1[,3])
//'     stopifnot(max(abs(rsrse[,2] - mertens_se),na.rm=TRUE) < 1e-12)
//' }
//'
//' @seealso \code{\link{scale}}
//' @template etc
//' @template ref-romo
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_centered(SEXP v, 
                               SEXP window = R_NilValue, 
                               Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                               bool na_rm=false, int min_df=0, double used_df=1.0, int lookahead=0, int restart_period=100,
                               bool check_wts=false, bool normalize_wts=false) {
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<centered>(v, wts, 1, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
}
// scale the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_scaled(SEXP v, SEXP window = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                             bool na_rm=false, int min_df=0, double used_df=1.0, int lookahead=0, int restart_period=100,
                             bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<scaled>(v, wts, 2, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
}
// zscore the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_zscored(SEXP v, SEXP window = R_NilValue, 
                              Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                              bool na_rm=false, int min_df=0, double used_df=1.0, int lookahead=0, int restart_period=100,
                              bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<zscore>(v, wts, 2, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
}
// sharpe on the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_sharpe(SEXP v, SEXP window = R_NilValue, 
                             Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                             bool na_rm=false, bool compute_se=false, int min_df=0, double used_df=1.0, int restart_period=100,
                             bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    if (compute_se) {
        return runningQMomentsCurryTwo<sharpese>(v, wts, 4, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
    } 
    return runningQMomentsCurryTwo<sharpe>(v, wts, 2, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
}
// t stat of the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_tstat(SEXP v, SEXP window = R_NilValue, 
                            Rcpp::Nullable< Rcpp::NumericVector > wts = R_NilValue, 
                            bool na_rm=false, int min_df=0, double used_df=1.0, int restart_period=100,
                            bool check_wts=false, bool normalize_wts=true) {
    int wins=get_wins(window);
    return runningQMomentsCurryTwo<tstat>(v, wts, 2, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
}


//' @title
//' Convert between different types of moments, raw, central, standardized.
//' @description
//' Given raw or central or standardized moments, convert to another type.
//' 
//' @param input a vector of the count, then the mean, then the \code{2} through \code{k}
//' raw or central moments.
//'
//' @template etc
//' @rdname moment_conversions
//' @export
// [[Rcpp::export]]
NumericVector cent2raw(NumericVector input) {
    int ord = input.length() - 1;
    NumericVector output(ord+1);
    int ppp,qqq;
    output[0] = input[0];
    if (ord > 0) { 
        output[1] = input[1];
        for (ppp=2;ppp<=ord;++ppp) {
            output[ppp] = pow(output[1],ppp);
            for (qqq=2;qqq<=ppp;++qqq) {
                output[ppp] += bincoef[ppp][qqq] * input[qqq] * pow(output[1],ppp-qqq);
            }
        }
    }
    return output;
}
//UNFOLD

// junkyard//FOLDUP

// code that *would* have been used in RcppModules,
// but the latter cannot deal properly with generics.


//typedef std::vector<double> dubvec;  // convenience typedef

//class Moments {
    //public:
        //// truly empty
        //Moments() : order(3), na_rm(false) {
            //moments = Rcpp::NumericVector(1+order);
        //}
        //// empty
        //Moments(int order) : order(order), na_rm(false) {
            //moments = Rcpp::NumericVector(1+order);
        //}
        //Moments(int order, bool na_rm) : order(order), na_rm(na_rm) {
            //moments = Rcpp::NumericVector(1+order);
        //}
        //// 'pure'
        //Moments(dubvec input, int order, bool na_rm) : order(order), na_rm(na_rm) {
            //moments = quasiMoments<dubvec>(input, order, 0, -1, na_rm);
        //}
        //// access the normalized moments
        //NumericVector cent_moments(int used_df=1) {
            //NumericVector vret = NumericVector(1+order);
            //vret[order] = moments[0];
            //vret[order-1] = moments[1];
            //for (int mmm=2;mmm <= order;mmm++) {
                //vret[order-mmm] = moments[mmm] / (moments[0] - used_df);
            //}
            //return vret;
        //}
        //// append operation; this is the guy that shits the bed.
        //// much fun. do break. amaze.
        //void join(const Moments& rhs) {
            //moments = joinMoments(moments,rhs.moments);
        //}
        
        //const int order;
        //const bool na_rm;
    //private:
        //NumericVector moments;
//};

//RCPP_MODULE(moment_module) {
    //class_<Moments>( "Moments" )

    //.constructor()
    //.constructor<int>()
    //.constructor<int,bool>()
    //.constructor<dubvec,int,bool>()

    //.field_readonly("order", &Moments::order)
    //.field_readonly("na_rm", &Moments::na_rm)

    //.method("cent_moments", &Moments::cent_moments)
    //.method("%:%", &Moments::join)
    //;
//}
//UNFOLD
//


//
// 2FIX:
//
// running means need weights !!! 
// compensated summation for weights where necessary
// no Kahans for running sum of integers or logicals
//
// 2FIX: make a compensated summation class/object ?
// 2FIX: inline code for adding a (weighted) observation ?
// 2FIX: default to infinity for windows, not NULL...

// notes on Rcpp
// http://thecoatlessprofessor.com/programming/rcpp/unofficial-rcpp-api-docs/#nan-constants
//
// types in R!
// INTSXP REALSXP CPLXSXP
//
// https://stackoverflow.com/questions/25172419/how-can-i-get-the-sexptype-of-an-sexp-value
// https://teuder.gitbooks.io/introduction-to-rcpp/en/07_data_types.html
// http://adv-r.had.co.nz/C-interface.html
//
// overload += in c++
// https://stackoverflow.com/questions/34663785/c-operator-overload
// https://stackoverflow.com/questions/4581961/c-how-to-overload-operator
// 

// for help on dispatch, see:
// http://stackoverflow.com/a/25254680/164611
// http://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/


//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
