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

  this is horrible; 
  originally this was a templated function, but it was too slow:
  flattening inlines did not work in terms of speed.
  so include this as a file.
  sorry.

  Created: 2024.11.20
  Copyright: Steven E. Pav, 2016-2024
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

// to keep it DRY, this function has a bunch of variants,
// depending on the templated bools in the vanilla form,
// this function returns a NumericMatrix with as many rows
// as elements in the input, and some number of columns.
//
// in certain forms, depending on templated retwhat, this computes
// the correlation, covariance, regression coefficient and so on.
//
// there is also a 'minimum df' parameter. this is the minimum count
// required to return data for some of these forms.

// assumes that template parameters retwhat, T, and renormalize are set.
// note that renormalize could be a template bool or an actual variable.
// and that xret, lll (rownum), frets  and min_df are set.
// do_interp<retwhat, Welford<oneW,has_wts,ord_beyond> ,T,renormalize>(xret,lll,ord,frets,v,used_df,min_df);

// assumes these have been defined:
//double denom;

    if (retwhat==ret_correlation) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.correlation();
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.correlation();
            } else {
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_regression_slope) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.regression_slope();
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.regression_slope();
            } else {
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_regression_intercept) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.regression_intercept();
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.regression_intercept();
            } else {
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_regression_fit) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                frets.assign_regression_fit(xret, lll);
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                frets.assign_regression_fit(xret, lll);
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_regression_diagnostics) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                frets.assign_regression_diagnostics(xret, lll, renormalize, used_df);
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
                xret(lll,2) = NAN;
                xret(lll,3) = NAN;
                xret(lll,4) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                frets.assign_regression_diagnostics(xret, lll, renormalize, used_df);
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
                xret(lll,2) = NAN;
                xret(lll,3) = NAN;
                xret(lll,4) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_covariance) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.covariance(renormalize, used_df);
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.covariance(renormalize, used_df);
            } else {
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_covariance_matrix) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                denom = frets.var_denominator(renormalize, used_df);
                xret(lll,0) = frets.m_xx[3] / denom;
                xret(lll,1) = frets.m_xx[4] / denom;
                xret(lll,2) = frets.m_xx[5] / denom;
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
                xret(lll,2) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                denom = frets.var_denominator(renormalize, used_df);
                xret(lll,0) = frets.m_xx[3] / denom;
                xret(lll,1) = frets.m_xx[4] / denom;
                xret(lll,2) = frets.m_xx[5] / denom;
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
                xret(lll,2) = NAN;
            }
        }
    } //UNFOLD

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
