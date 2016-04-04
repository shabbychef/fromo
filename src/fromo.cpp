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

// univariate sums, moments, cumulants//FOLDUP

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
    NumericVector retv;
    // FML, I cannot figure out how to get v.length()
    switch (TYPEOF(v)) {
        case  INTSXP: { retv = quasiMoments<IntegerVector>(v, ord, 0, -1, na_rm); break; }
        case REALSXP: { retv = quasiMoments<NumericVector>(v, ord, 0, -1, na_rm); break; }
        case  LGLSXP: { retv = quasiMoments<LogicalVector>(v, ord, 0, -1, na_rm); break; }
        default: stop("Unsupported input type"); // nocov
    }
    return retv;
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
    NumericVector retv;
    switch (TYPEOF(v)) {
        case  INTSXP: { retv = quasiWelford<IntegerVector,IntegerVector::iterator>(v, na_rm); break; }
        case REALSXP: { retv = quasiWelford<NumericVector,NumericVector::iterator>(v, na_rm); break; }
        case  LGLSXP: { retv = quasiWelford<LogicalVector,LogicalVector::iterator>(v, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return retv;
}

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
//' @return a vector, filled out as follows:
//' \describe{
//' \item{sd3}{A vector of the (sample) standard devation, mean, and number of elements.}
//' \item{skew4}{A vector of the (sample) skewness, standard devation, mean, and number of elements.}
//' \item{kurt5}{A vector of the (sample) excess kurtosis, skewness, standard devation, mean, and number of elements.}
//' \item{cent_moments}{A vector of the (sample) \eqn{k}th centered moment, then \eqn{k-1}th centered moment, ..., 
//'  then the \emph{variance}, the mean, and number of elements.}
//' \item{std_moments}{A vector of the (sample) \eqn{k}th standardized (and centered) moment, then 
//'  \eqn{k-1}th, ..., then standard devation, mean, and number of elements.}
//' \item{cent_cumulants}{A vector of the (sample) \eqn{k}th (centered, but this is redundant) cumulant, then the \eqn{k-1}th, ...,
//'  then the \emph{variance} (which is the second cumulant), the mean, and number of elements.}
//' }
//'
//' @note
//' The first centered (and standardized) moment is often defined to be identically 0. Instead \code{cent_moments}
//' and \code{std_moments} returns the mean. 
//' Similarly, the second standardized moments defined to be identically 1; \code{std_moments} instead returns the standard
//' deviation. The reason is that a user can always decide to ignore the results and fill in a 0 or 1 as they need, but 
//' could not efficiently compute the mean and standard deviation from scratch if we discard it.
//' 
//' @note
//' The last minus two element of the output of \code{cent_moments} and \code{cent_cumulants} is the \emph{variance},
//' not the standard deviation. All other code return the standard deviation in that place.
//'
//' @note
//' The kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @note
//' 'centered cumulants' is redundant. The intent was to avoid collision with existing code named 'cumulants'.
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
//'
//' @template etc
//' @template ref-romo
//' @rdname firstmoments
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
NumericVector cent_moments(SEXP v, int max_order=5, int used_df=0, bool na_rm=false) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = wrapMoments(v, max_order, na_rm);
    NumericVector vret = sums2revm(preval,(double)used_df);
    return vret;
}
// return the standardized moments up to order max_order
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector std_moments(SEXP v, int max_order=5, int used_df=0, bool na_rm=false) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = wrapMoments(v, max_order, na_rm);
    double sigma = COMP_SD_TWO(preval,(double)used_df);
    int mmm;
    NumericVector vret(max_order+1);
    for (mmm=0;mmm <= MIN(1,max_order);++mmm) {
        vret[max_order-mmm] = preval[mmm];
    }
    if (max_order > 1) {
        vret[max_order-2] = sigma;
        for (mmm=3;mmm <= max_order;++mmm) {
            vret[max_order-mmm] = preval[mmm] / pow(sigma,((double)mmm)/2.0);
        }
    }
    return vret;
}
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector cent_cumulants(SEXP v, int max_order=5, int used_df=0, bool na_rm=false) {
    NumericVector cmoms = cent_moments(v,max_order,used_df,na_rm);
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
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector cent_sums(SEXP v, int max_order=5, bool na_rm=false) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = wrapMoments(v, max_order, na_rm);
    return preval;
}
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector join_cent_sums(NumericVector ret1,NumericVector ret2) {
    double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
    double ac_nfoo,ac_n2,ac_mn1;
    double ac_del,ac_mn2,ac_n1;

    if (ret1.size() != ret2.size()) { stop("mismatch in sizes."); }

    int ord = ret1.size() - 1;
    int ppp,qqq;

    n1 = ret1[0];
    if (n1 <= 0) { return ret2; }
    n2 = ret2[0];
    if (n2 <= 0) { return ret1; }
    
    // allocate output and copy;
    NumericVector ret3(ord+1);
    for (ppp=0;ppp<=ord;++ppp) { ret3[ppp] = ret1[ppp]; }

    ret3[0] += n2;
    ntot = ret3[0];
    n1rat = n1 / ntot;
    n2rat = n2 / ntot;
    del21 = ret2[1] - ret3[1];
    mupart = del21 * n2rat;

    ret3[1] += mupart;
    nfoo = n1 * mupart;
    ac_nfoo = pow(nfoo,ord);
    ac_n2 = pow(n2,1-ord);
    ac_mn1 = pow(-n1,1-ord);
    for (ppp=ord;ppp >= 2;ppp--) {
        ret3[ppp] += ret2[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
        if (ppp > 2) {
            if (nfoo != 0) { ac_nfoo /= nfoo; }
            ac_n2 *= n2;
            ac_mn1 *= (-n1);
        }
        ac_del = del21;
        ac_mn2 = -n2rat;
        ac_n1 = n1rat;
        for (int qqq=1;qqq <= (ppp-2); qqq++) {
            ret3[ppp] += bincoef[ppp][qqq] * ac_del * (ac_mn2 * ret3[ppp-qqq] + ac_n1 * ret2[ppp-qqq]);
            if (qqq < (ppp-2)) {
                ac_del *= del21;
                ac_mn2 *= (-n2rat);
                ac_n1  *= (n1rat);
            }
        }
    }
    return ret3;
}
//' @rdname centsums 
//' @export
// [[Rcpp::export]]
NumericVector unjoin_cent_sums(NumericVector ret3,NumericVector ret2) {
    // subtract ret2 from ret3
    double n1, n2, ntot, del21, mupart, nfoo, n1rat, n2rat;
    double ac_nfoo,ac_n2,ac_mn1;
    double ac_del,ac_mn2,ac_n1;

    if (ret3.size() != ret2.size()) { stop("mismatch in sizes."); }

    int ord = ret3.size() - 1;
    int ppp,qqq;

    n1 = ret3[0];
    n2 = ret2[0];
    if (n2 <= 0) { return ret3; }
    if (n2 > n1) { stop("cannot subtract more observations than were seen."); }

    // allocate output 
    NumericVector ret1(ord+1);
    // would be better to check they are equal, but whatever;
    if (n2 == n1) { 
        // ret1 is all zero, and should be, so return it
        return ret1;
    } else {
        // copy
        for (ppp=0;ppp<=ord;++ppp) { ret1[ppp] = ret3[ppp]; }
    }

    mupart = ret2[1] - ret1[1];

    ntot = ret1[0];
    ret1[0] -= n2;
    n1 = ret1[0];

    n1rat = n1 / ntot;
    n2rat = n2 / ntot;

    ret1[1] -= (n2/n1) * mupart;

    del21 = mupart / n1rat;
    nfoo = mupart * n2;

    ac_nfoo = nfoo * nfoo;
    ac_n2 = 1.0 / n2;
    ac_mn1 = -1.0 / n1;
    for (ppp=2;ppp <= ord;ppp++) {
        ret1[ppp] -= ret2[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
        if (ppp < ord) {
            ac_nfoo *= nfoo; 
            ac_n2 /= n2;
            ac_mn1 /= (-n1);
        }
        ac_del = del21;
        ac_mn2 = -n2rat;
        ac_n1 = n1rat;
        for (int qqq=1;qqq <= (ppp-2); qqq++) {
            ret1[ppp] -= bincoef[ppp][qqq] * ac_del * (ac_mn2 * ret1[ppp-qqq] + ac_n1 * ret2[ppp-qqq]);
            if (qqq < (ppp-2)) {
                ac_del *= del21;
                ac_mn2 *= (-n2rat);
                ac_n1  *= (n1rat);
            }
        }
    }
    return ret1;
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
        stop("malformed input");
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
        stop("malformed input");
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
// ret_extreme return a rows x 2 matrix of the count and the maximum centered sum
// ret_extreme return a rows x 2 matrix of the count and the maximum centered sum
template <typename T,bool ret_mat,bool ret_extreme,bool ret_cent,bool ret_scald,bool ret_z,bool ret_sr,bool ret_t,bool ret_srmer>
NumericMatrix runningQMoments(T v,
                              int ord = 3,
                              int window = NA_INTEGER,
                              int recom_period = 100, 
                              int lookahead = 0,
                              const int min_df = 0,
                              bool na_rm = false) {
    double nextv, prevv, compv, nel, nelm, della, delnel, drat, ac_dn, ac_on, ac_de;

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(window);
    if ((window < 1) && (!infwin)) { stop("must give positive window"); }

    if (min_df < 0) { stop("require positive min_df"); }
    if (!infwin && (min_df > window)) { stop("must have min_df <= window"); }

    if (((ret_scald || ret_z || ret_sr || ret_t) && (ord < 2)) ||
        ((ret_srmer) && (ord < 4))) { 
        stop("bad code: order too small to support this computation"); 
    }
    // only for ret_srmer, but cannot define outside its scope.
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
    if (ret_mat) { 
        ncols = 1+ord; 
    } else if (ret_extreme || ret_srmer) {
        ncols = 2; 
    } else {
        ncols = 1; 
    }
        
    NumericMatrix xret(numel,ncols);
    // this is the current estimate, fill it in as we go.
    NumericVector vret(1+ord);

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
                vret = quasiMoments<T>(v, ord, jjj, iii + 1, na_rm);
            }
            subcount = 0;
        } else {
            if ((tr_iii < numel) && (tr_iii >= 0)) {
                // add on nextv:
                nextv = double(v[tr_iii]);
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
            }
            // remove prevv:
            if ((tr_jjj < numel) && (tr_jjj >= 0)) {
                prevv = double(v[tr_jjj]);
                if (! (na_rm && ISNAN(prevv))) {//FOLDUP
                    della = prevv - vret[1];
                    nel = --vret[0];
                    nelm = nel + 1;

                    if (nel > 0) {
                        delnel = della / nel;
                        vret[1] -= delnel;
                        // correct delta
                        della = delnel * nelm;

                        drat = delnel * nel;
                        ac_dn = drat * drat;
                        ac_on = - 1.0 / (nel);

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
                    // increment the subcount counter
                    subcount++;
                    // check for Heywood cases and recompute if hit.//FOLDUP
                    for (ppp=2;ppp <= ord;ppp += 2) {
                        if (vret[ppp] <= 0.0) {
                            iii = MIN(numel-1,tr_iii);
                            jjj = MAX(0,tr_jjj+1);
                            if (jjj <= iii) {
                                vret = quasiMoments<T>(v, ord, jjj, iii + 1, na_rm);
                            }
                            subcount = 0;
                            break;
                        }
                    }//UNFOLD
                }//UNFOLD
            }
        }
        tr_jjj++;

        // fill in the value in the output.//FOLDUP
        if (ret_mat) {
            if (vret[0] >= min_df) {
                // put them in backwards!
                if (vret[0] >= ord) {
                    for (mmm=0;mmm <= ord;++mmm) {
                        xret(lll,ord-mmm) = vret[mmm];
                    }
                } else {
                    for (mmm=0;mmm <= vret[0];++mmm) {
                        xret(lll,ord-mmm) = vret[mmm];
                    }
                    for (mmm=vret[0]+1;mmm <= ord;++mmm) {
                        xret(lll,ord-mmm) = NAN;
                    }
                }
            } else {
                for (mmm=0;mmm <= ord;++mmm) {
                    xret(lll,mmm) = NAN;
                }
            }
        } else if (ret_extreme) {
            if (vret[0] >= min_df) {
                // put them in backwards!
                if (vret[0] >= ord) {
                    xret(lll,0) = vret[ord];
                    xret(lll,1) = vret[0];
                } else {
                    xret(lll,0) = vret[ord];
                    xret(lll,1) = NAN;
                }
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
            }
        } else {
            if (vret[0] >= min_df) {
                if (ret_cent) {
                    compv = double(v[lll]);
                    xret(lll,0) = COMP_CENTERED(compv,vret);
                }
                if (ret_scald) {
                    compv = double(v[lll]);
                    xret(lll,0) = (compv) / COMP_SD(vret);
                }
                if (ret_z) {
                    compv = double(v[lll]);
                    xret(lll,0) = COMP_CENTERED(compv,vret) / COMP_SD(vret);
                }
                if (ret_sr) {
                    xret(lll,0) = COMP_SHARPE(vret);
                }
                if (ret_t) {
                    xret(lll,0) = (vret[1]) / (sqrt(vret[2] / (vret[0] * (vret[0]-1.0))));
                }
                if (ret_srmer) {
                    skew = COMP_SKEW(vret);
                    exkurt = COMP_EXKURT(vret);
                    sr = COMP_SHARPE(vret);
                    xret(lll,0) = sr;
                    xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / vret[0]);
                }
            } else {
                xret(lll,0) = NAN;
            }
        }
        //UNFOLD
    }

    return xret;
}

// wrap the call:
NumericMatrix wrapRunningQMoments(SEXP v, int ord, int window, int recom_period, const int min_df, bool na_rm, bool max_order_only) {
    if (max_order_only) {
        switch (TYPEOF(v)) {
            case  INTSXP: { return runningQMoments<IntegerVector, false, true, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            case REALSXP: { return runningQMoments<NumericVector, false, true, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            case  LGLSXP: { return runningQMoments<LogicalVector, false, true, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            default: stop("Unsupported input type");
        }
    } else {
        switch (TYPEOF(v)) {
            case  INTSXP: { return runningQMoments<IntegerVector, true, false, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            case REALSXP: { return runningQMoments<NumericVector, true, false, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            case  LGLSXP: { return runningQMoments<LogicalVector, true, false, false, false, false, false, false, false>(v, ord, window, recom_period, 0, min_df, na_rm); }
            default: stop("Unsupported input type");
        }
    }
    // CRAN checks are broken: 'warning: control reaches end of non-void function'
    // ... only for a crappy automated warning.
    return NumericMatrix(1,1);
}

// helper function; takes a double or integer windowsize and interprets w/out warning or vomit.
// if NULL, return NA_INTEGER;
// if integer, pass through as i ;
// if double, then 
//   if Inf or NA, return NA_INTEGER;
//   else convert to integer via as<int>( )
int get_wins(SEXP window) {
    if (Rf_isNull(window)) { return NA_INTEGER; }
    switch (TYPEOF(window)) {
        case  INTSXP: { return as<int>(window); break; }
        case REALSXP: { 
                          double wins = as<double>(window);
                          if ((NumericVector::is_na(wins)) || 
                              (traits::is_infinite<REALSXP>(wins) && (wins > 0.0))) {
                              return NA_INTEGER;
                          }
                          return (int)wins;
                          break;
                      }
        default: stop("Unsupported input type");
    }
    // CRAN checks are broken: 'warning: control reaches end of non-void function'
    // ... only for a crappy automated warning.
    return -1;
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
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_sd3(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int used_df=1, int restart_period=100) {
    int wins=get_wins(window);
    double udf=(double)used_df;
    NumericMatrix preval = wrapRunningQMoments(v, 2, wins, restart_period, min_df, na_rm, false);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = sqrt(preval(iii,0)/(preval(iii,2)-udf));
    }
    return preval;
}

// return the skew, the standard deviation, the mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_skew4(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int used_df=1, int restart_period=100) {
    int wins=get_wins(window);
    double udf=(double)used_df;
    NumericMatrix preval = wrapRunningQMoments(v, 3, wins, restart_period, min_df, na_rm, false);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = sqrt(preval(iii,3)) * preval(iii,0) / pow(preval(iii,1),1.5);
        preval(iii,1) = sqrt(preval(iii,1)/(preval(iii,3)-udf));
    }
    return preval;
}

// return the //excess// kurtosis, skew, standard deviation, mean, and the dof
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_kurt5(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int used_df=1, int restart_period=100) {
    int wins=get_wins(window);
    double udf=(double)used_df;
    NumericMatrix preval = wrapRunningQMoments(v, 4, wins, restart_period, min_df, na_rm, false);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = (preval(iii,4) * preval(iii,0) / pow(preval(iii,2),2.0)) - 3.0;
        preval(iii,1) = sqrt(preval(iii,4)) * preval(iii,1) / pow(preval(iii,2),1.5);
        preval(iii,2) = sqrt(preval(iii,2)/(preval(iii,4)-udf));
    }
    return preval;
}
// return the centered moments down to the 2nd, then the mean, and the dof.
//' @param max_order_only for \code{running_cent_moments}, if this flag is set, only compute
//' the maximum order centered moment, and return in a vector.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_cent_moments(SEXP v, SEXP window = R_NilValue, int max_order=5, bool na_rm=false, bool max_order_only=false, 
                                   int min_df=0, int used_df=0, int restart_period=100) {
    int wins=get_wins(window);
    int ncolm1,ncolm2;
    double denom;
    double udf=(double)used_df;
    NumericMatrix preval = wrapRunningQMoments(v, max_order, wins, restart_period, min_df, na_rm, max_order_only);
    ncolm1 = preval.ncol() - 1;
    ncolm2 = ncolm1 - 1;
    if (max_order > 1) {
        if (max_order_only) {
            // 2FIX: this is not as space efficient as it could be,
            // and a bit of a mess!
            NumericMatrix retval(Dimension(preval.nrow(),1));
            for (int iii=0;iii < preval.nrow();++iii) {
                denom = preval(iii,1) - udf;
                retval(iii,0) = preval(iii,0) / denom;
            }
            return retval;
        } else {
            // fix the higher than mean columns;
            for (int iii=0;iii < preval.nrow();++iii) {
                denom = preval(iii,ncolm1) - udf;
                for (int mmm=0;mmm < ncolm2;++mmm) {
                    preval(iii,mmm) = preval(iii,mmm) / denom;
                }
            }
        }
    }
    return preval;
}

// return the standardized moments down to the 3rd, then the standard deviation, the mean, and the dof.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_std_moments(SEXP v, SEXP window = R_NilValue, int max_order=5, bool na_rm=false, int min_df=0, int used_df=0, int restart_period=100) {
    int wins=get_wins(window);
    double denom;
    double udf = (double)used_df;
    double sigma;
    NumericMatrix preval = wrapRunningQMoments(v, max_order, wins, restart_period, min_df, na_rm, false);
    if (max_order > 1) {
        // fix the higher than mean columns;
        for (int iii=0;iii < preval.nrow();++iii) {
            denom = preval(iii,max_order) - udf;
            sigma = sqrt(preval(iii,max_order-2) / denom);
            preval(iii,max_order-2) = sigma;

            for (int mmm=0;mmm < (max_order-2);++mmm) {
                preval(iii,mmm) = preval(iii,mmm) / (denom * pow(sigma,max_order-mmm));
            }
        }
    }
    return preval;
}

// return the cumulants down to the 2nd, then the standard deviation, the mean, and the dof.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix running_cumulants(SEXP v, SEXP window = R_NilValue, int max_order=5, bool na_rm=false, int min_df=0, int used_df=0, int restart_period=100) {
    NumericMatrix cumulants = running_cent_moments(v, window, max_order, na_rm, false, min_df, used_df, restart_period);
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
//' AS 269. \url{https://web.archive.org/web/20110103030111/http://lib.stat.cmu.edu/apstat/269}
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
//' @rdname runningquantiles
//' @export
// [[Rcpp::export]]
NumericMatrix running_apx_quantiles(SEXP v, NumericVector p, SEXP window = R_NilValue, int max_order=5, bool na_rm=false, int min_df=0, int used_df=0, int restart_period=100) {
    NumericMatrix cumulants = running_cumulants(v, window, max_order, na_rm, min_df, used_df, restart_period);
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
NumericMatrix running_apx_median(SEXP v, SEXP window = R_NilValue, int max_order=5, bool na_rm=false, int min_df=0, int used_df=0, int restart_period=100) {
    NumericVector p(1);
    p(0) = 0.5;
    NumericMatrix vret = running_apx_quantiles(v,p,window,max_order,na_rm,min_df,used_df,restart_period);
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
NumericMatrix running_centered(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int lookahead=0, int restart_period=100) {
    NumericMatrix preval;
    int wins=get_wins(window);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, true, false, false, false, false, false>(v, 1, wins, restart_period, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, true, false, false, false, false, false>(v, 1, wins, restart_period, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, true, false, false, false, false, false>(v, 1, wins, restart_period, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// scale the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_scaled(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int lookahead=0, int restart_period=100) {
    NumericMatrix preval;
    int wins=get_wins(window);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, true, false, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, true, false, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, true, false, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// zscore the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_zscored(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int lookahead=0, int restart_period=100) {
    NumericMatrix preval;
    int wins=get_wins(window);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, true, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, true, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, true, false, false, false>(v, 2, wins, restart_period, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// sharpe on the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_sharpe(SEXP v, SEXP window = R_NilValue, bool na_rm=false, bool compute_se=false, int min_df=0, int restart_period=100) {
    NumericMatrix preval;
    int wins=get_wins(window);
    if (compute_se) {
        switch (TYPEOF(v)) {
            case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, false, false, false, true>(v, 4, wins, restart_period, 0, min_df, na_rm); break; }
            case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, false, false, false, true>(v, 4, wins, restart_period, 0, min_df, na_rm); break; }
            case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, false, false, false, true>(v, 4, wins, restart_period, 0, min_df, na_rm); break; }
            default: stop("Unsupported input type");
        }
    } else {
        switch (TYPEOF(v)) {
            case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, false, true, false, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
            case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, false, true, false, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
            case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, false, true, false, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
            default: stop("Unsupported input type");
        }
    }
    return preval;
}
// t stat of the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix running_tstat(SEXP v, SEXP window = R_NilValue, bool na_rm=false, int min_df=0, int restart_period=100) {
    NumericMatrix preval;
    int wins=get_wins(window);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, false, false, true, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, false, false, true, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, false, false, true, false>(v, 2, wins, restart_period, 0, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
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

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
