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

  running sum and mean functions.

  Created: 2018.12.27
  Copyright: Steven E. Pav, 2016-2018
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

// 2FIX: rename to utils.cpp ?
#ifndef __DEF_COMMON_CODE__
#define __DEF_COMMON_CODE__

#include <cmath>
#include <string.h>
#include "common.h"

#include <Rcpp.h>

#endif /* __DEF_COMMON_CODE__ */

using namespace Rcpp;

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

double get_double_wins(SEXP window) {
	if (!Rf_isNull(window)) {
		switch (TYPEOF(window)) {
			case  INTSXP: { 
											return as<double>(window); }
			case REALSXP: { 
											double wins = as<double>(window);
											if ((NumericVector::is_na(wins)) || 
													(traits::is_infinite<REALSXP>(wins) && (wins > 0.0))) {
												return NA_REAL;
											}
											return wins;
										}
			default: stop("Unsupported input type");
		}
	}
	return NA_REAL;
}



// modifies in place.
void centmom2cumulants(NumericMatrix cumulants,int max_order) {
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
}
// cf:
//
// Lee, Y-S., and Lin, T-K. "Algorithm AS269: High Order Cornish Fisher
// Expansion." Appl. Stat. 41, no. 1 (1992): 233-240. 
// \url{http://www.jstor.org/stable/2347649}
//
// Lee, Y-S., and Lin, T-K. "Correction to Algorithm AS269: High Order 
// Cornish Fisher Expansion." Appl. Stat. 42, no. 1 (1993): 268-269. 
// \url{http://www.jstor.org/stable/2347433}
//
// AS 269. \url{http://lib.stat.cmu.edu/apstat/269}
//
// turns cumulants into approximate quantiles
NumericMatrix cumulants2quantiles(NumericMatrix cumulants,NumericVector p, int max_order) {
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


//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
