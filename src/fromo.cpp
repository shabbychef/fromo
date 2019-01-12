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

#include "common.h"
#include "kahan.h"
#include "welford.h"
#include "running.h"

#endif /* __DEF_FROMO__ */

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// univariate sums, moments, cumulants//FOLDUP

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
//' Moments are computed as follows, given some values \eqn{x_i} and optional weights \eqn{w_i},
//' defaulting to 1, the weighted mean is computed as
//' \deqn{\mu = \frac{\sum_i x_i w_i}{\sum w_i}.}
//' The weighted kth central sum is computed as
//' \deqn{\mu = \sum_i \left(x_i - \mu\right)^k w_i.}
//' Let \eqn{n = \sum_i w_i} be the sum of weights (or number of observations in the unweighted case).
//' Then the weighted kth central moment is computed as that weighted sum divided by the
//' adjusted sum weights:
//' \deqn{\mu_k = \frac{\sum_i \left(x_i - \mu\right)^k w_i}{n - \nu},}
//' where \eqn{\nu} is the \sQuote{used df}, provided by the user to adjust the denominator.
//' (Typical values are 0 or 1.)
//' The weighted kth standardized moment is the central moment divided by the second central moment
//' to the \eqn{k/2} power:
//' \deqn{\tilde{\mu}_k = \frac{\mu_k}{\mu_2^{k/2}}.}
//' The (centered) rth cumulant, for \eqn{r \ge 2} is then computed using the formula of Willink, namely
//' \deqn{\kappa_r = \mu_r - \sum_{j=0}^{r - 2} {r - 1 \choose j} \mu_j \kappa {r-j}.}
//' The standardized rth cumulant is the rth centered cumulant divided by \eqn{\mu_2^{r/2}}.
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
//' The term 'centered cumulants' is redundant. The intent was to avoid possible collision with existing code named 'cumulants'.
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
//'
//'
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
//' # check that skew matches moments::skewness
//' if (require(moments)) {
//'     set.seed(1234)
//'     x <- rnorm(1000)
//'     resu <- fromo::skew4(x)
//' 
//'     msku <- moments::skewness(x)
//'     stopifnot(abs(msku - resu[1]) < 1e-14)
//' }
//' 
//' # check skew vs e1071 skewness, which has a different denominator
//' if (require(e1071)) {
//'     set.seed(1234)
//'     x <- rnorm(1000)
//'     resu <- fromo::skew4(x)
//' 
//'     esku <- e1071::skewness(x,type=3)
//'     nobs <- resu[4]
//'     stopifnot(abs(esku - resu[1] * ((nobs-1)/nobs)^(3/2)) < 1e-14)
//'
//'     # similarly:
//'     resu <- fromo::std_moments(x,max_order=3,used_df=0)
//'     stopifnot(abs(esku - resu[1] * ((nobs-1)/nobs)^(3/2)) < 1e-14)
//' }
//'
//' @template etc
//' @template ref-romo
//' @template ref-cumulants
//' @template param-wts
//' @template note-wts
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
    // WTF, why does it not use the used_df? 
    
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
        adj = cmoms[max_order - 2];
        sigma = sqrt(adj);
        cmoms(max_order-2) = sigma;  // put the stdev in
        for (mmm=3;mmm <= max_order;++mmm) {
            adj *= sigma;
            cmoms[max_order-mmm] /= adj;
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
    Welford<double,true,true,true> frets1 = Welford<double,true,true,true>(ord,cret1);
    Welford<double,true,true,true> frets2 = Welford<double,true,true,true>(ord,cret2);
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
    Welford<double,true,true,true> frets3 = Welford<double,true,true,true>(ord,cret3);
    Welford<double,true,true,true> frets2 = Welford<double,true,true,true>(ord,cret2);
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


//
// 2FIX:
//
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
     
//UNFOLD

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
