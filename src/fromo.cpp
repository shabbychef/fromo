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
#define MAX(a,b) ((a>b)? a:b)
#define MIN(a,b) ((a<b)? a:b)

#endif /* __DEF_FROMO__ */


#include <Rcpp.h>
using namespace Rcpp;

// this function returns a NumericVector of:
//   the number of elements, 
//   the mean, and 
//   an (ord - 1)-vector consisting of the
//   2nd through ord'th centered sum, defined
//   as sum_j (v[j] - mean)^i
// if top < 0, take the length of v.
template <typename T>
NumericVector quasiMoments(T v,
                           int ord = 3,
                           int bottom = 0,
                           int top = -1,
                           bool na_rm = false) {
    double nextv, nel, nelm, della, delnel, drat, ac_dn, ac_on, ac_de;

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // preallocated with zeros:
    NumericVector xret(1+ord);
    
    if ((top < 0) || (top > v.size())) {
        top = v.size(); 
    }

    for (int iii=bottom;iii < top;++iii) {
        nextv = v[iii];
        if (! (na_rm && ISNAN(nextv))) {
            della = nextv - xret[1];
            nelm = xret[0];
            nel = ++xret[0];
            delnel = della / nel;
            xret[1] += delnel;

            if (nelm > 0) {
                drat = delnel * nelm;
                ac_dn = pow(drat,ord);
                ac_on = pow(-1.0 / nelm,ord-1);

                for (int ppp=ord;ppp >= 2;ppp--) {
                    xret[ppp] += ac_dn * (1.0 - ac_on);
                    if (ppp > 2) {
                        if (drat != 0) { ac_dn /= drat; }
                        ac_on *= -nelm;
                    }
                    ac_de = -delnel;
                    for (int qqq=1;qqq <= ppp-2;qqq++) {
                        xret[ppp] += bincoef[ppp][qqq] * ac_de * xret[ppp-qqq];
                        if (qqq < ppp - 2) {
                            ac_de *= -delnel;
                        }
                    }
                }
            }
        }
    }

    return xret;
}

// wrap the call:
NumericVector wrapMoments(SEXP v, int ord, bool na_rm) {
    // FML, I cannot figure out how to get v.length()
    switch (TYPEOF(v)) {
        case  INTSXP: { return quasiMoments<IntegerVector>(v, ord, 0, -1, na_rm); }
        case REALSXP: { return quasiMoments<NumericVector>(v, ord, 0, -1, na_rm); }
        case  LGLSXP: { return quasiMoments<LogicalVector>(v, ord, 0, -1, na_rm); }
        default: stop("Unsupported input type");
    }
}

// specialization of the above for the case of ord=2
// this function returns a NumericVector of:
//   the number of elements, 
//   the mean, and 
//   an (ord - 1)-vector consisting of the
//   2nd through ord'th centered sum, defined
//   as sum_j (v[j] - mean)^i
template <typename T,typename iterT>
NumericVector quasiWelford(T v,
                           bool na_rm = false) {
    double nextv, nel, nelm, della, delnel;
    const int ord=2;

    // preallocated with zeros:
    NumericVector xret(1+ord);

    for (iterT it = v.begin(); it != v.end(); ++it) {
        nextv = double(*it);
        if (! (na_rm && ISNAN(nextv))) {
            della = nextv - xret[1];
            nelm = xret[0];
            nel = ++xret[0];
            delnel = della / nel;
            xret[1] += delnel;
            xret[2] += delnel * delnel * nelm * (nelm + 1.0);
        }
    }

    return xret;
}

// wrap the call:
NumericVector wrapWelford(SEXP v, bool na_rm) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return quasiWelford<IntegerVector,IntegerVector::iterator>(v, na_rm); }
        case REALSXP: { return quasiWelford<NumericVector,NumericVector::iterator>(v, na_rm); }
        case  LGLSXP: { return quasiWelford<LogicalVector,LogicalVector::iterator>(v, na_rm); }
        default: stop("Unsupported input type");
    }
}

#define COMP_SD(preval) (sqrt(preval[2]/(preval[0]-1.0)))
#define COMP_SKEW(preval) (sqrt(preval[0]) * preval[3] / pow(preval[2],1.5))
#define COMP_EXKURT(preval) ((preval[0] * preval[4] / (pow(preval[2],2.0))) - 3.0)

// for help on dispatch, see:
// http://stackoverflow.com/a/25254680/164611
// http://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
// compute the standard deviation and mean and # df
//' @title
//' Compute first K moments
//' @description
//' Compute the (standardized) 2nd through kth moments, the mean, and the number of elements.
//' 
//' @param v a vector
//' @param na_rm whether to remove NA, false by default.
//' @param used_df the number of degrees of freedom consumed, used in the denominator
//' of the centered moments computation. These are subtracted from the number of
//' observations.
//' @param max_order the maximum order of the centered moment to be computed.
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
//' @return a vector; the first elements are the kth, k-1th through 2nd standardized, centered moment,
//' then the mean, then the number of (non-nan) elements in the input.
//'
//' @note
//' the kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @examples
//' x <- rnorm(1e5)
//' sd3(x)[1] - sd(x)
//'   skew4(x)[4] - length(x)
//'   skew4(x)[3] - mean(x)
//'   skew4(x)[2] - sd(x)
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
//'
//' @template etc
//' @template ref-romo
//' @rdname firstmoments
//' @seealso runningmoments
//' @export
// [[Rcpp::export]]
NumericVector sd3(SEXP v, bool na_rm=false) {
    //NumericVector preval = wrapMoments(v, 2, na_rm);
    NumericVector preval = wrapWelford(v, na_rm);
    NumericVector vret = NumericVector::create(COMP_SD(preval),
                                               preval[1],
                                               preval[0]);

    return vret;
}

// return the skew, the standard deviation, the mean, and the dof
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector skew4(SEXP v, bool na_rm=false) {
    NumericVector preval = wrapMoments(v, 3, na_rm);
    NumericVector vret = NumericVector::create(COMP_SKEW(preval),
                                               COMP_SD(preval),
                                               preval[1],
                                               preval[0]);
    return vret;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector kurt5(SEXP v, bool na_rm=false) {
    NumericVector preval = wrapMoments(v, 4, na_rm);
    NumericVector vret = NumericVector::create(COMP_EXKURT(preval),
                                               COMP_SKEW(preval),
                                               COMP_SD(preval),
                                               preval[1],
                                               preval[0]);
    return vret;
}

// return the centered moments up to order max_order
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector cent_moments(SEXP v, int max_order=5, int used_df=1, bool na_rm=false) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = wrapMoments(v, max_order, na_rm);
    NumericVector vret = NumericVector(1+max_order);
    vret[max_order] = preval[0];
    vret[max_order-1] = preval[1];
    for (int mmm=2;mmm <= max_order;mmm++) {
        vret[max_order-mmm] = preval[mmm] / (preval[0] - used_df);
    }
    return vret;
}

/*******************************************************************************
 * running moments
 */

// to keep it DRY, this function has a bunch of variants,
// depending on the templated bools in the vanilla form,
// this function returns a NumericMatrix with as many rows
// as elements in the input, and ord+1 columns.
// unlike the quasiMoments code, the *last* column
// is the number of elements,
// the last minus one is the mean and so on;
// this simplifies the transformation to moments later.
//
//   the first column is the number of elements, 
//   the second is the mean,
//   the (ord + 1 - k)th is the k-1th centered sum, defined as
//   as sum_j (v[j] - mean)^i
//
// moreover we adapt a sliding window of size winsize.
// for computational efficiency, we add and subtract
// observations. this can lead to roundoff issues,
// especially in the subtraction of observations.
// the algorithm checks for negative second moment and
// starts afresh when encountered. Also, the computation
// is periodically restarted.
//
// in other forms, depending on templated bools, this
// computes the centered input, the rescaled input, the z-scored input
// or a t-scored input, as matrices with a single column.
template <typename T,bool ret_mat,bool ret_cent,bool ret_scald,bool ret_z,bool ret_t>
NumericMatrix runningQMoments(T v,
                              int ord = 3,
                              int winsize = NA_INTEGER,
                              int recom_period = 100, 
                              bool na_rm = false) {
    double nextv, prevv, nel, nelm, della, delnel, drat, ac_dn, ac_on, ac_de;

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    int iii,jjj,mmm;
    int runsize = 0;
    int numel = v.size();
    const bool infwin = IntegerVector::is_na(winsize);

    // preallocated with zeros; should
    // probably be NA?
    int ncols;
    if (ret_mat) {
        ncols = 1+ord;
    } else {
        ncols = 1;
    }
    NumericMatrix xret(numel,ncols);
    // this is the current estimate, fill it in as we go.
    NumericVector vret(1+ord);

    jjj = 0;
    for (iii=0;iii < numel;++iii) {
        nextv = double(v[iii]);
        if (! (na_rm && ISNAN(nextv))) {//FOLDUP
            della = nextv - vret[1];
            nelm = vret[0];
            nel = ++vret[0];
            delnel = della / nel;
            vret[1] += delnel;

            if (nelm > 0) {//FOLDUP
                drat = delnel * nelm;
                ac_dn = pow(drat,ord);
                ac_on = pow(-1.0 / nelm,ord-1);

                for (int ppp=ord;ppp >= 2;ppp--) {
                    vret[ppp] += ac_dn * (1.0 - ac_on);
                    if (ppp > 2) {
                        if (drat != 0) { ac_dn /= drat; }
                        ac_on *= -nelm;
                    }
                    ac_de = -delnel;
                    for (int qqq=1;qqq <= ppp-2;qqq++) {
                        vret[ppp] += bincoef[ppp][qqq] * ac_de * vret[ppp-qqq];
                        if (qqq < ppp - 2) {
                            ac_de *= -delnel;
                        }
                    }
                }
            }//UNFOLD
        }//UNFOLD
        // possibly remove old observations. wheee
        if ((!infwin) && (iii >= winsize)) {
            if (runsize >= recom_period) {
                vret = quasiMoments<T>(v, ord, iii - winsize + 1, iii+1, na_rm); 
                runsize = 0;
            } else {
                prevv = double(v[jjj]);
                if (! (na_rm && ISNAN(prevv))) {//FOLDUP
                    della = prevv - vret[1];
                    nelm = vret[0];
                    nel = --vret[0];

                    if (nel > 0) {
                        delnel = della / nel;
                        vret[1] -= delnel;
                        // correct delta
                        della = delnel * nelm;

                        drat = delnel * nel;
                        ac_dn = pow(drat,2);
                        ac_on = - (nel);

                        for (int ppp=2;ppp <= ord;ppp++) {
                            vret[ppp] -= ac_dn * (1.0 - ac_on);
                            if (ppp < ord) {
                                ac_dn *= drat;
                                ac_on /= -nel;
                            }
                            ac_de = -delnel;
                            for (int qqq=1;qqq <= ppp-2;qqq++) {
                                vret[ppp] -= bincoef[ppp][qqq] * ac_de * vret[ppp-qqq];
                                if (qqq < ppp - 2) {
                                    ac_de *= -delnel;
                                }
                            }
                        }
                    } else {
                        // else nel <= 0, meaning there are no observations.  reset vret to all zero
                        for (mmm=0;mmm <= ord;++mmm) { vret[mmm] = 0.0; }
                    }

                    runsize++;
                }//UNFOLD
            }
            // check for Heywood cases and recompute.//FOLDUP
            if (((ord > 1) && (vret[2] <= 0.0)) || ((ord > 3) && (vret[4] <= 0.0))) {
                vret = quasiMoments<T>(v, ord, iii - winsize + 1, iii+1, na_rm); 
                runsize = 0;
            }//UNFOLD
            jjj++;
        }
        // fill in the value in the output.//FOLDUP
        // backwards!
        if (ret_mat) {
            if (vret[0] >= ord) {
                for (mmm=0;mmm <= ord;++mmm) {
                    xret(iii,ord-mmm) = vret[mmm];
                }
            } else {
                for (mmm=0;mmm <= vret[0];++mmm) {
                    xret(iii,ord-mmm) = vret[mmm];
                }
                for (mmm=vret[0]+1;mmm <= ord;++mmm) {
                    xret(iii,ord-mmm) = NAN;
                }
            }
        } 
        if (ret_cent) {
            xret(iii,0) = nextv - vret[ord-1];
        }
        if (ret_scald) {
            xret(iii,0) = (nextv) / (sqrt(vret[ord-2] / (vret[ord]-1.0)));
        }
        if (ret_z) {
            xret(iii,0) = (nextv - vret[ord-1]) / (sqrt(vret[ord-2] / (vret[ord]-1.0)));
        }
        if (ret_t) {
            xret(iii,0) = (vret[ord-1]) / (sqrt(vret[ord-2] / (vret[ord] * (vret[ord]-1.0))));
        }
        //UNFOLD
    }

    return xret;
}

// wrap the call:
NumericMatrix wrapRunningQMoments(SEXP v, int ord, int winsize, int recom_period, bool na_rm) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return runningQMoments<IntegerVector, true, false, false, false, false>(v, ord, winsize, recom_period, na_rm); }
        case REALSXP: { return runningQMoments<NumericVector, true, false, false, false, false>(v, ord, winsize, recom_period, na_rm); }
        case  LGLSXP: { return runningQMoments<LogicalVector, true, false, false, false, false>(v, ord, winsize, recom_period, na_rm); }
        default: stop("Unsupported input type");
    }
}

//' @title
//' Compute first K moments over a sliding window
//' @description
//' Compute the (standardized) 2nd through kth moments, the mean, and the number of elements over
//' an infinite or finite sliding window, returning a matrix.
//' 
//' @param v a vector
//' @param winsize the window size. if NA, equivalent to infinity.
//' @param recoper the recompute period. because subtraction of elements can cause
//' loss of precision, the computation of moments is restarted periodically based on 
//' this parameter. Larger values mean fewer restarts and faster, though less accurate
//' results. Note that the code checks for negative second and fourth moments and
//' recomputes when needed.
//' @param na_rm whether to remove NA, false by default.
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
//' @return a matrix; the first columns are the kth, k-1th through 2nd standardized, centered moment,
//' then a column of the mean, then a column of the number of (non-nan) elements in the input.
//' When there are not sufficient (non-nan) elements for the computation, NaN are returned.
//'
//' @note
//' the kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @examples
//' x <- rnorm(1e5)
//' run_sd3(x,10)
//' run_skew4(x,10)
//' run_kurt5(x,500)
//'
//' @template etc
//' @template ref-romo
//' @seealso firstmoments
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_sd3(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval = wrapRunningQMoments(v, 2, winsize, recoper, na_rm);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = sqrt(preval(iii,0)/(preval(iii,2)-1.0));
    }
    return preval;
}

// return the skew, the standard deviation, the mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_skew4(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval = wrapRunningQMoments(v, 3, winsize, recoper, na_rm);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = sqrt(preval(iii,3)) * preval(iii,0) / pow(preval(iii,1),1.5);
        preval(iii,1) = sqrt(preval(iii,1)/(preval(iii,3)-1.0));
    }
    return preval;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_kurt5(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval = wrapRunningQMoments(v, 4, winsize, recoper, na_rm);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = (preval(iii,4) * preval(iii,0) / pow(preval(iii,2),2.0)) - 3.0;
        preval(iii,1) = sqrt(preval(iii,4)) * preval(iii,1) / pow(preval(iii,2),1.5);
        preval(iii,2) = sqrt(preval(iii,2)/(preval(iii,4)-1.0));
    }
    return preval;
}

// center the input
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_centered(SEXP v, int winsize=NA_INTEGER, int recoper=1000, bool na_rm=false) {
    NumericMatrix preval;
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, true, false, false, false>(v, 1, winsize, recoper, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, true, false, false, false>(v, 1, winsize, recoper, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, true, false, false, false>(v, 1, winsize, recoper, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// scale the input
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_scaled(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval;
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, true, false, false>(v, 2, winsize, recoper, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, true, false, false>(v, 2, winsize, recoper, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, true, false, false>(v, 2, winsize, recoper, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// zscore the input
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_zscored(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval;
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, true, false>(v, 2, winsize, recoper, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, true, false>(v, 2, winsize, recoper, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, true, false>(v, 2, winsize, recoper, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// tscore the input
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_tscored(SEXP v, int winsize=NA_INTEGER, int recoper=100, bool na_rm=false) {
    NumericMatrix preval;
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, true>(v, 2, winsize, recoper, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, true>(v, 2, winsize, recoper, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, true>(v, 2, winsize, recoper, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
