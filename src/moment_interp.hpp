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

  this is horrible; including this file because flattening inlines is
  not working...

  Created: 2019.01.02
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

// assumes that template parameters retwhat, T, and renormalize are set.
// and that xret, lll (rownum), ord, frets, v, used_df and min_df are set.
// do_interp<retwhat, Welford<oneW,has_wts,ord_beyond> ,T,renormalize>(xret,lll,ord,frets,v,used_df,min_df);

    if (retwhat==ret_centmoments) {//FOLDUP
        double sg_denom,renorm;
        double mydf,dwsum;
        int mmm;
        dwsum = double(frets.wsum());
        if (renormalize) { 
            mydf = double(frets.nel());
            renorm = mydf / dwsum;
        } else {
            mydf = dwsum;
        }

        if (mydf >= min_df) {
            sg_denom = mydf - used_df;
            if (renormalize) { sg_denom /= renorm; }
            xret(lll,ord) = mydf; 
            xret(lll,ord-1) = frets.m_xx[1];

            // put them in backwards!
            if (mydf >= ord) {
                if (ord >= 2) {
                    xret(lll,ord-2) = frets.m_xx[2] / sg_denom;
                    for (mmm=3;mmm <= ord;++mmm) {
                        xret(lll,ord-mmm) = frets.m_xx[mmm] / dwsum;
                    }
                }
            } else {
                if (ord >= 2) {
                    xret(lll,ord-2) = frets.m_xx[2] / sg_denom;
                    for (mmm=3;mmm <= mydf;++mmm) {
                        xret(lll,ord-mmm) = frets.m_xx[mmm] / dwsum;
                    }
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
    if (retwhat==ret_stdmoments) { //FOLDUP
        double sg_denom,renorm;
        double mydf,dwsum;
        double sigmasq,sigma,sigmapow;
        int mmm;
        dwsum = double(frets.wsum());
        if (renormalize) { 
            mydf = double(frets.nel());
            renorm = mydf / dwsum;
        } else {
            mydf = dwsum;
        }

        if (mydf >= min_df) {
            sg_denom = mydf - used_df;
            if (renormalize) { sg_denom /= renorm; }
            sigmasq = frets.m_xx[2] / sg_denom;
            sigma = sqrt(sigmasq);
            xret(lll,ord) = mydf; 
            xret(lll,ord-1) = frets.m_xx[1];
            xret(lll,ord-2) = sigma;

            // put them in backwards!
            if (mydf >= ord) {
                for (mmm=3;mmm <= ord;++mmm) {
                    sigmasq *= sigma;
                    xret(lll,ord-mmm) = frets.m_xx[mmm] / (dwsum * sigmasq);
                }
            } else {
                for (mmm=3;mmm <= mydf;++mmm) {
                    sigmasq *= sigma;
                    xret(lll,ord-mmm) = frets.m_xx[mmm] / (dwsum * sigmasq);
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
    }
    //UNFOLD
    if (retwhat==ret_sd3) { //FOLDUP
        if (renormalize) {
            double sg_denom,renorm;
            double mydf,dwsum,sigma;
            dwsum = double(frets.wsum());
            mydf = double(frets.nel());
            renorm = mydf / dwsum;

            if (mydf >= min_df) {
                xret(lll,2) = mydf; 

                // put them in backwards!
                if (mydf >= ord) {
                    sg_denom = (mydf - used_df) / renorm;
                    sigma = sqrt(frets.m_xx[2] / sg_denom);

                    xret(lll,1) = frets.m_xx[1];
                    xret(lll,0) = sigma;
                } else {
                    if (mydf >= 1) { 
                        xret(lll,1) = frets.m_xx[1]; 
                    } else  {
                        xret(lll,1) = NAN;
                    }
                    xret(lll,0) = NAN;
                }
            } else {
                xret(lll,2) = NAN;
                xret(lll,1) = NAN;
                xret(lll,0) = NAN;
            }
        }
        if (!renormalize) {
            double sg_denom,mydf,sigma;
            mydf = double(frets.wsum());

            if (mydf >= min_df) {
                // put them in backwards!
                xret(lll,2) = mydf;
                if (mydf >= ord) {
                    sigma = sqrt(frets.m_xx[2] / (mydf - used_df));
                    xret(lll,1) = frets.m_xx[1];
                    xret(lll,0) = sigma;
                } else {
                    if (mydf >= 1) { 
                        xret(lll,1) = frets.m_xx[1]; 
                    } else  {
                        xret(lll,1) = NAN;
                    }
                    xret(lll,0) = NAN;
                }
            } else {
                xret(lll,2) = NAN;
                xret(lll,1) = NAN;
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_skew4) { //FOLDUP
        double sg_denom,renorm;
        double mydf,dwsum;
        double sigmasq,sigma,sigmapow;
        int mmm;
        dwsum = double(frets.wsum());
        if (renormalize) { 
            mydf = double(frets.nel());
            renorm = mydf / dwsum;
        } else {
            mydf = dwsum;
        }

        if (mydf >= min_df) {
            sg_denom = mydf - used_df;
            if (renormalize) { sg_denom /= renorm; }
            sigmasq = frets.m_xx[2] / sg_denom;
            sigma = sqrt(sigmasq);

            // put them in backwards!
            if (mydf >= ord) {
                xret(lll,3) = mydf; 
                xret(lll,2) = frets.m_xx[1];
                xret(lll,1) = sigma;
                xret(lll,0) = COMP_SKEW_TWO(frets.m_xx,dwsum);
            } else {
                xret(lll,3) = mydf; 
                if (mydf >= 1) {
                    xret(lll,2) = frets.m_xx[1];
                    if (mydf >= 2) {
                        xret(lll,1) = sigma;
                    } else {
                        xret(lll,1) = NAN;
                    }
                } else {
                    xret(lll,2) = NAN;
                    xret(lll,1) = NAN;
                }
                xret(lll,0) = NAN;
            }
        } else {
            xret(lll,3) = NAN;
            xret(lll,2) = NAN;
            xret(lll,1) = NAN;
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_exkurt5) { //FOLDUP
        double sg_denom,renorm;
        double mydf,dwsum;
        double sigmasq,sigma,sigmapow;
        int mmm;
        dwsum = double(frets.wsum());
        if (renormalize) { 
            mydf = double(frets.nel());
            renorm = mydf / dwsum;
        } else {
            mydf = dwsum;
        }

        if (mydf >= min_df) {
            sg_denom = mydf - used_df;
            if (renormalize) { sg_denom /= renorm; }
            sigmasq = frets.m_xx[2] / sg_denom;
            sigma = sqrt(sigmasq);

            // put them in backwards!
            if (mydf >= ord) {
                xret(lll,4) = mydf; 
                xret(lll,3) = frets.m_xx[1];
                xret(lll,2) = sigma;
                // uhoh! renormalization!
                xret(lll,1) = COMP_SKEW_TWO(frets.m_xx,dwsum);
                // uhoh! renormalization!
                xret(lll,0) = (COMP_KURT_TWO(frets.m_xx,dwsum)) - 3.0;
            } else {
                xret(lll,4) = mydf; 
                if (mydf >= 1) {
                    xret(lll,3) = frets.m_xx[1];
                    if (mydf >= 2) {
                        xret(lll,2) = sigma;
                        if (mydf >= 3) {
                            // uhoh! renormalization!
                            xret(lll,1) = COMP_SKEW_TWO(frets.m_xx,dwsum);
                        } else {
                            xret(lll,1) = NAN;
                        }
                    } else {
                        xret(lll,2) = NAN;
                        xret(lll,1) = NAN;
                    }
                } else {
                    xret(lll,3) = NAN;
                    xret(lll,2) = NAN;
                    xret(lll,1) = NAN;
                }
                xret(lll,0) = NAN;
            }
        } else {
            xret(lll,4) = NAN;
            xret(lll,3) = NAN;
            xret(lll,2) = NAN;
            xret(lll,1) = NAN;
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_centmaxonly) { //FOLDUP
        double denom,renorm;
        double mydf,dwsum;
        int mmm;
        dwsum = double(frets.wsum());
        if (renormalize) { 
            mydf = double(frets.nel());
            renorm = mydf / dwsum;
        } else {
            mydf = dwsum;
        }

        if ((mydf >= min_df) && (mydf >= ord)) {
            if (ord==2) {
                denom = mydf - used_df;
            } else {
                denom = dwsum;
            }
            if (renormalize) { 
                xret(lll,0) = renorm * frets.m_xx[ord] / denom;
            } else {
                xret(lll,0) = frets.m_xx[ord] / denom;
            }
        } else {
            xret(lll,0) = NAN;
        }
    } //UNFOLD
    if (retwhat==ret_centered) {//FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.centered(double(v[lll]));
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_scaled) {//FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.scaled(double(v[lll]),renormalize,used_df);
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_zscore) { //FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.zscored(double(v[lll]),renormalize,used_df);
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    //if (retwhat==ret_tstat) {//FOLDUP
    //if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
    //if (renormalize) {
    //xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.nel()));
    //} else {
    //xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.wsum()));
    //}
    //} else {
    //xret(lll,0) = NAN;
    //}
    //}//UNFOLD
    if (retwhat==ret_sharpe) { //FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.sharpe(renormalize,used_df); 
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_sharpese) { //FOLDUP
        double skew,exkurt,sr;
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            skew = frets.skew();
            exkurt = frets.exkurt();
            sr = frets.sharpe(renormalize,used_df);
            xret(lll,0) = sr;
            if (renormalize) {
                xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / double(frets.nel()));
            } else {
                xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / double(frets.wsum()));
            }
        } else {
            xret(lll,0) = NAN;
            xret(lll,1) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_stdev) { //FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.sd(renormalize,used_df); 
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD
    if (retwhat==ret_skew) { //FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.skew();
        } else {
            xret(lll,0) = NAN;
        }
    } //UNFOLD
    if (retwhat==ret_exkurt) { //FOLDUP
        if ((!renormalize && (frets.wsum() >= min_df)) || (renormalize && (frets.nel() >= min_df))) {
            xret(lll,0) = frets.exkurt();
        } else {
            xret(lll,0) = NAN;
        }
    }//UNFOLD

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
