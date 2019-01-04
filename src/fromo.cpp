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

// running sums, moments, cumulants, approximate quantiles//FOLDUP

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
    NumericMatrix preval = runQMCurryThree<ret_sd3>(v, wts, 2, wins, restart_period, 0, min_df, used_df, 
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
    NumericMatrix preval = runQMCurryThree<ret_skew4>(v, wts, 3, wins, restart_period, 0, min_df, used_df, 
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
    NumericMatrix preval = runQMCurryThree<ret_exkurt5>(v, wts, 4, wins, restart_period, 0, min_df, used_df, 
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
    return runQMCurryThree<ret_stdev>(v, wts, 2, wins, restart_period, 0, min_df, used_df,
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
    return runQMCurryThree<ret_skew>(v, wts, 3, wins, restart_period, 0, min_df, used_df, 
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
    return runQMCurryThree<ret_exkurt>(v, wts, 4, wins, restart_period, 0, min_df, used_df,
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
        return runQMCurryThree<ret_centmaxonly>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
                                                    na_rm, check_wts, normalize_wts);
    } 
    return runQMCurryThree<ret_centmoments>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
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
    return runQMCurryThree<ret_stdmoments>(v, wts, max_order, wins, restart_period, 0, min_df, used_df, 
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
    return runQMCurryThree<ret_centered>(v, wts, 1, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
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
    return runQMCurryThree<ret_scaled>(v, wts, 2, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
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
    return runQMCurryThree<ret_zscore>(v, wts, 2, wins, restart_period, lookahead, min_df, used_df, na_rm, check_wts, normalize_wts);
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
        return runQMCurryThree<ret_sharpese>(v, wts, 4, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
    } 
    return runQMCurryThree<ret_sharpe>(v, wts, 2, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
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
    return runQMCurryThree<ret_tstat>(v, wts, 2, wins, restart_period, 0, min_df, used_df, na_rm, check_wts, normalize_wts);
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
