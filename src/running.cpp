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

#ifndef __DEF_RUNCODE__
#define __DEF_RUNCODE__

#include "common.h"
#include "kahan.h"
#include "welford.h"
#include "running.h"

#endif /* __DEF_RUNCODE__ */

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

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
//UNFOLD

// for help on dispatch, see:
// http://stackoverflow.com/a/25254680/164611
// http://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
     
//UNFOLD

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
