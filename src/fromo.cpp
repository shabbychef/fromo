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

#define COMP_CENTERED(x,preval) (x - preval[1])

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
NumericVector cent_moments(SEXP v, int max_order=5, int used_df=0, bool na_rm=false) {
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

// return the centered moments up to order max_order
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
NumericVector raw_sums(SEXP v, int max_order=5, bool na_rm=false) {
    if (max_order < 1) { stop("must give largeish max_order"); }
    NumericVector preval = wrapMoments(v, max_order, na_rm);
    return preval;
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
//
// we have a lookahead option for the centered, scaled, and Z-scored
// variants. Positive lookahead means take info from the future.
//
// there is also a 'minimum df' parameter. this is the minimum count
// required to return data for the ret_cent, _scald, _z, and _t forms.
// the reasoning is that you might not want the z-score on fweer than 10
// observations. these do the right thing wrt NA, BTW. some moments
// come out as zero when computed on too few observations, and we blindly
// return Inf or NaN in that case. set the min_df to correct for this.
template <typename T,bool ret_mat,bool ret_cent,bool ret_scald,bool ret_z,bool ret_t>
NumericMatrix runningQMoments(T v,
                              int ord = 3,
                              int winsize = NA_INTEGER,
                              int recom_period = 100, 
                              int lookahead = 0,
                              const int min_df = 0,
                              bool na_rm = false) {
    double nextv, prevv, compv, nel, nelm, della, delnel, drat, ac_dn, ac_on, ac_de;

    if (ord < 1) { stop("require positive order"); }
    if (ord > MAX_ORD) { stop("too many moments requested, weirdo"); }

    // 2FIX: later you should use the infwin to prevent some computations
    // from happening. like subtracting old observations, say.
    const bool infwin = IntegerVector::is_na(winsize);
    if ((winsize < 1) && (!infwin)) { stop("must give positive winsize"); }

    if (min_df < 0) { stop("require positive min_df"); }
    if (!infwin && (min_df > winsize)) { stop("must have min_df <= winsize"); }

    int iii,jjj,lll,mmm,ppp,qqq,tr_iii,tr_jjj;
    int numel = v.size();
    const bool non_aligned = (lookahead != 0);
    // refers to the number of *subtractions* performed
    int subcount = 0;

    // preallocated with zeros; should
    // probably be NA?
    int ncols;
    if (ret_mat) { ncols = 1+ord; } else { ncols = 1; }
    NumericMatrix xret(numel,ncols);
    // this is the current estimate, fill it in as we go.
    NumericVector vret(1+ord);

    // as an invariant, we will start the computation
    // with vret, which is initialized as the summed
    // means on [jjj,iii]
    tr_iii = lookahead - 1;
    tr_jjj = lookahead - winsize;
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
                    // check for Heywood cases and recompute.//FOLDUP
                    if (((ord > 1) && (vret[2] <= 0.0)) || ((ord > 3) && (vret[4] <= 0.0))) {
                        iii = MIN(numel-1,tr_iii);
                        jjj = MAX(0,tr_jjj+1);
                        if (jjj <= iii) {
                            vret = quasiMoments<T>(v, ord, jjj, iii + 1, na_rm);
                        }
                        subcount = 0;
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
                if (ret_t) {
                    xret(lll,0) = (vret[1]) / (sqrt(vret[2] / (vret[0] * (vret[0]-1.0))));
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
NumericMatrix wrapRunningQMoments(SEXP v, int ord, int winsize, int recom_period, const int min_df, bool na_rm) {
    switch (TYPEOF(v)) {
        case  INTSXP: { return runningQMoments<IntegerVector, true, false, false, false, false>(v, ord, winsize, recom_period, 0, min_df, na_rm); }
        case REALSXP: { return runningQMoments<NumericVector, true, false, false, false, false>(v, ord, winsize, recom_period, 0, min_df, na_rm); }
        case  LGLSXP: { return runningQMoments<LogicalVector, true, false, false, false, false>(v, ord, winsize, recom_period, 0, min_df, na_rm); }
        default: stop("Unsupported input type");
    }
}

// helper function; takes a double or integer windowsize and interprets w/out warning or vomit.
// if NULL, return NA_INTEGER;
// if integer, pass through as i ;
// if double, then 
//   if Inf or NA, return NA_INTEGER;
//   else convert to integer via as<int>( )
int get_wins(SEXP winsize) {
    if (Rf_isNull(winsize)) { return NA_INTEGER; }
    switch (TYPEOF(winsize)) {
        case  INTSXP: { return as<int>(winsize); break; }
        case REALSXP: { 
                          double wins = as<double>(winsize);
                          if ((NumericVector::is_na(wins)) || 
                              (traits::is_infinite<REALSXP>(wins) && (wins > 0.0))) {
                              return NA_INTEGER;
                          }
                          return (int)wins;
                          break;
                      }
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
//' @param winsize the window size. if given as finite integer or double, passed through.
//' If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, equivalent
//' to an infinite window size. If negative, an error will be thrown.
//' @param recoper the recompute period. because subtraction of elements can cause
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
//' centered standardized moment, for \eqn{k=2,3,4}{k=2,3,4}. These
//' are computed via the numerically robust one-pass method of Bennett \emph{et. al.}
//'
//' Given the length \eqn{n} vector \eqn{x}, we output matrix \eqn{M} where
//' \eqn{M_{i,j}}{M_i,j} is the \eqn{order - j + 1} moment (\emph{i.e.}
//' excess kurtosis, skewness, standard deviation, mean or number of elements)
//' of \eqn{x_{i-winsize+1},x_{i-winsize+2},...,x_{i}}{x_(i-winsize+1),x_(i-winsize+2),...,x_i}.
//' Barring \code{NA} or \code{NaN}, this is over a window of size \code{winsize}.
//' During the 'burn-in' phase, we take fewer elements.
//'
//' @return a matrix; the first columns are the kth, k-1th through 2nd standardized, centered moment,
//' then a column of the mean, then a column of the number of (non-nan) elements in the input.
//' When there are not sufficient (non-nan) elements for the computation, \code{NaN} are returned.
//'
//' @note
//' the kurtosis is \emph{excess kurtosis}, with a 3 subtracted, and should be nearly zero
//' for Gaussian input.
//'
//' @examples
//' x <- rnorm(1e5)
//' run_sd3(x,10)
//' run_skew4(x,10)
//'
//' if (require(moments)) {
//'     set.seed(123)
//'     x <- rnorm(5e1)
//'     winsize <- 10L
//'     kt5 <- run_kurt5(x,winsize=winsize)
//'     rm1 <- t(sapply(seq_len(length(x)),function(iii) { 
//'                 xrang <- x[max(1,iii-winsize+1):iii]
//'                 c(moments::kurtosis(xrang)-3.0,moments::skewness(xrang),
//'                 sd(xrang),mean(xrang),length(xrang)) },
//'              simplify=TRUE))
//'     stopifnot(max(abs(kt5 - rm1),na.rm=TRUE) < 1e-12)
//' }
//'
//' run_cent_moments(x,winsize=100L,max_order=6L)
//'
//' @template etc
//' @template ref-romo
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_sd3(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int min_df=0, bool na_rm=false) {
    int wins=get_wins(winsize);
    NumericMatrix preval = wrapRunningQMoments(v, 2, wins, recoper, min_df, na_rm);
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
NumericMatrix run_skew4(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int min_df=0, bool na_rm=false) {
    int wins=get_wins(winsize);
    NumericMatrix preval = wrapRunningQMoments(v, 3, wins, recoper, min_df, na_rm);
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
NumericMatrix run_kurt5(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int min_df=0, bool na_rm=false) {
    int wins=get_wins(winsize);
    NumericMatrix preval = wrapRunningQMoments(v, 4, wins, recoper, min_df, na_rm);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        preval(iii,0) = (preval(iii,4) * preval(iii,0) / pow(preval(iii,2),2.0)) - 3.0;
        preval(iii,1) = sqrt(preval(iii,4)) * preval(iii,1) / pow(preval(iii,2),1.5);
        preval(iii,2) = sqrt(preval(iii,2)/(preval(iii,4)-1.0));
    }
    return preval;
}
// return the centered moments down to the 2nd, then the mean, and the dof.
//' @rdname runningmoments
//' @export
// [[Rcpp::export]]
NumericMatrix run_cent_moments(SEXP v, SEXP winsize = R_NilValue, int max_order=5, int recoper=100, int min_df=0, int used_df=0, bool na_rm=false) {
    int wins=get_wins(winsize);
    double denom;
    double udf = (double)used_df;
    NumericMatrix preval = wrapRunningQMoments(v, max_order, wins, recoper, min_df, na_rm);
    // fix the higher than mean columns;
    for (int iii=0;iii < preval.nrow();++iii) {
        denom = preval(iii,0) - udf;
        for (int mmm=0;mmm < (max_order-1);++mmm) {
            preval(iii,mmm) = preval(iii,mmm) / denom;
        }
    }
    return preval;
}

//' @title
//' Compare data to moments computed over a sliding window.
//' @description
//' Computes moments over a sliding window, then adjusts the data accordingly, centering, or scaling,
//' or z-scoring, and so on.
//' 
//' @inheritParams run_cent_moments
//' @param min_df the minimum df to return a value, otherwise \code{NaN} is returned.
//' This can be used to prevent \emph{e.g.} Z-scores from being computed on only 3
//' observations. Defaults to zero, meaning no restriction, which can result in 
//' infinite Z-scores during the burn-in period.
//' @param lookahead for some of the operations, the value is compared to 
//' mean and standard deviation possibly using 'future' or 'past' information
//' by means of a non-zero lookahead. Positive values mean data are taken from
//' the future.
//'
//' @details
//'
//' Given the length \eqn{n} vector \eqn{x}, for
//' a given index \eqn{i}, define \eqn{x^{(i)}}{x^(i)}
//' as the vector of 
//' \eqn{x_{i-winsize+1},x_{i-winsize+2},...,x_{i}}{x_(i-winsize+1),x_(i-winsize+2),...,x_i},
//' where we do not run over the 'edge' of the vector. In code, this is essentially
//' \code{x[(max(1,i-winsize+1)):i]}. Then define \eqn{\mu_i}{mu_i}, \eqn{\sigma_i}{sigma_i}
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
//'     winsize <- 10L
//'     rm1 <- t(sapply(seq_len(length(x)),function(iii) { 
//'                   xrang <- x[max(1,iii-winsize+1):iii]
//'                   c(sd(xrang),mean(xrang),length(xrang)) },
//'                   simplify=TRUE))
//'     rcent <- run_centered(x,winsize=winsize)
//'     rscal <- run_scaled(x,winsize=winsize)
//'     rzsco <- run_zscored(x,winsize=winsize)
//'     rtsco <- run_tscored(x,winsize=winsize)
//'     stopifnot(max(abs(rcent - (x - rm1[,2])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rscal - (x / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rzsco - ((x - rm1[,2]) / rm1[,1])),na.rm=TRUE) < 1e-12)
//'     stopifnot(max(abs(rtsco - ((sqrt(rm1[,3]) * rm1[,2]) / rm1[,1])),na.rm=TRUE) < 1e-12)
//' }
//'
//' @template etc
//' @template ref-romo
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix run_centered(SEXP v, SEXP winsize = R_NilValue, int recoper=1000, int lookahead=0, int min_df=0, bool na_rm=false) {
    NumericMatrix preval;
    int wins=get_wins(winsize);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, true, false, false, false>(v, 1, wins, recoper, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, true, false, false, false>(v, 1, wins, recoper, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, true, false, false, false>(v, 1, wins, recoper, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// scale the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix run_scaled(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int lookahead=0, int min_df=0, bool na_rm=false) {
    NumericMatrix preval;
    int wins=get_wins(winsize);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, true, false, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, true, false, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, true, false, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// zscore the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix run_zscored(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int lookahead=0, int min_df=0, bool na_rm=false) {
    NumericMatrix preval;
    int wins=get_wins(winsize);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, true, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, true, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, true, false>(v, 2, wins, recoper, lookahead, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}
// tscore the input
//' @rdname runningadjustments
//' @export
// [[Rcpp::export]]
NumericMatrix run_tscored(SEXP v, SEXP winsize = R_NilValue, int recoper=100, int min_df=0, bool na_rm=false) {
    NumericMatrix preval;
    int wins=get_wins(winsize);
    switch (TYPEOF(v)) {
        case  INTSXP: { preval = runningQMoments<IntegerVector, false, false, false, false, true>(v, 2, wins, recoper, 0, min_df, na_rm); break; }
        case REALSXP: { preval = runningQMoments<NumericVector, false, false, false, false, true>(v, 2, wins, recoper, 0, min_df, na_rm); break; }
        case  LGLSXP: { preval = runningQMoments<LogicalVector, false, false, false, false, true>(v, 2, wins, recoper, 0, min_df, na_rm); break; }
        default: stop("Unsupported input type");
    }
    return preval;
}

//' @title
//' Join or unjoin moments computations.
//'
//' @param ret1 an \eqn{ord+1} vector as output by \code{\link{raw_sums}}? consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//' @param ret2 an \eqn{ord+1} vector as output by \code{\link{raw_sums}}? consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//' @param ret3 an \eqn{ord+1} vector as output by \code{\link{raw_sums}}? consisting of
//' the count, the mean, then the k through ordth centered sum of some observations.
//'
//' @details
//'
//' merge or unmerge sums of centered variables.
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
//'  rs1 <- raw_sums(x1,max_ord)
//'  rs2 <- raw_sums(x2,max_ord)
//'  rs3 <- raw_sums(c(x1,x2),max_ord)
//'  rs3alt <- join_moments(rs1,rs2)
//'  stopifnot(max(abs(rs3 - rs3alt)) < 1e-7)
//'  rs1alt <- unjoin_moments(rs3,rs2)
//'  rs2alt <- unjoin_moments(rs3,rs1)
//'  stopifnot(max(abs(rs1 - rs1alt)) < 1e-7)
//'  stopifnot(max(abs(rs2 - rs2alt)) < 1e-7)
//'
//' @template etc
//' @template ref-romo
//' @rdname joinmoments 
//' @export
// [[Rcpp::export]]
NumericVector join_moments(NumericVector ret1,NumericVector ret2) {
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

    NumericVector vret(ord+1);
    // copy
    for (ppp=0;ppp <= ord;++ppp) { vret[ppp] = ret1[ppp]; }

    vret[0] += n2;
    ntot = vret[0];
    n1rat = n1 / ntot;
    n2rat = n2 / ntot;
    del21 = ret2[1] - vret[1];
    mupart = del21 * n2rat;

    vret[1] += mupart;
    nfoo = n1 * mupart;
    ac_nfoo = pow(nfoo,ord);
    ac_n2 = pow(n2,1-ord);
    ac_mn1 = pow(-n1,1-ord);
    for (ppp=ord;ppp >= 2;ppp--) {
        //vret[ppp] = vret[ppp] + ret2[ppp] + (pow(nfoo,ppp) * (pow(n2,1-ppp) - pow(-n1,1-ppp)));
        vret[ppp] += ret2[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
        if (ppp > 2) {
            if (nfoo != 0) { ac_nfoo /= nfoo; }
            ac_n2 *= n2;
            ac_mn1 *= (-n1);
        }
        ac_del = del21;
        ac_mn2 = -n2rat;
        ac_n1 = n1rat;
        for (int qqq=1;qqq <= (ppp-2); qqq++) {
            //vret[ppp] += bincoef[ppp][qqq] * pow(del21,qqq) * (pow(-n2/ntot,qqq) * vret[ppp-qqq] + pow(n1/ntot,qqq) * ret2[ppp-qqq]);
            vret[ppp] += bincoef[ppp][qqq] * ac_del * (ac_mn2 * vret[ppp-qqq] + ac_n1 * ret2[ppp-qqq]);
            if (qqq < (ppp-2)) {
                ac_del *= del21;
                ac_mn2 *= (-n2rat);
                ac_n1  *= (n1rat);
            }
        }
    }
    return vret;
}
//' @rdname joinmoments 
//' @export
// [[Rcpp::export]]
NumericVector unjoin_moments(NumericVector ret3,NumericVector ret2) {
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

    NumericVector vret(ord+1);
    // would be better to check they are equal, but whatever;
    if (n2 == n1) { 
        // vret is all zero, just return it
        return vret;
    } else {
        // else copy
        for (ppp=0;ppp <= ord;++ppp) { vret[ppp] = ret3[ppp]; }
    }

    mupart = ret2[1] - vret[1];

    ntot = vret[0];
    vret[0] -= n2;
    n1 = vret[0];

    n1rat = n1 / ntot;
    n2rat = n2 / ntot;

    vret[1] -= (n2/n1) * mupart;

    del21 = mupart / n1rat;
    nfoo = mupart * n2;

    ac_nfoo = nfoo * nfoo;
    ac_n2 = 1.0 / n2;
    ac_mn1 = -1.0 / n1;
    for (ppp=2;ppp <= ord;ppp++) {
        vret[ppp] -= ret2[ppp] + (ac_nfoo * (ac_n2 - ac_mn1));
        if (ppp < ord) {
            ac_nfoo *= nfoo; 
            ac_n2 /= n2;
            ac_mn1 /= (-n1);
        }
        ac_del = del21;
        ac_mn2 = -n2rat;
        ac_n1 = n1rat;
        for (int qqq=1;qqq <= (ppp-2); qqq++) {
            vret[ppp] -= bincoef[ppp][qqq] * ac_del * (ac_mn2 * vret[ppp-qqq] + ac_n1 * ret2[ppp-qqq]);
            if (qqq < (ppp-2)) {
                ac_del *= del21;
                ac_mn2 *= (-n2rat);
                ac_n1  *= (n1rat);
            }
        }
    }
    return vret;
}

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


//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
