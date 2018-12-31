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

#ifndef __DEF_COMMON__
#define __DEF_COMMON__

#include <cmath>
#include <string.h>

#include <Rcpp.h>

#endif /* __DEF_COMMON__ */

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

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
