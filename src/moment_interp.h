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

  Created: 2019.01.02
  Copyright: Steven E. Pav, 2016-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

// to keep it DRY, this function has a bunch of variants,
// depending on the templated bools in the vanilla form,
// this function returns a NumericMatrix with as many rows
// as elements in the input, and ord+1 columns.
// unlike the quasiWeightedMoments code, the *last* column
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
//
// pipe in used_df

// assumes that template parameters retwhat, T, and renormalize are set.
// note that renormalize could be a template bool or an actual variable.
// and that xret, lll (rownum), ord, frets, v, used_df and min_df are set.
// do_interp<retwhat, Welford<oneW,has_wts,ord_beyond> ,T,renormalize>(xret,lll,ord,frets,v,used_df,min_df);

// assumes these have been defined:
//double sg_denom,renorm,denom,sigmasq,sigma,sigmapow,mydf,dwsum,skew,exkurt,sr;
//int mmm;

// 2FIX: the various checks here for mydf >= # seem to go against
// the min_df checks???

    if (retwhat==ret_centmoments) {//FOLDUP
        //double sg_denom,renorm;
        //double mydf,dwsum;
        //int mmm;
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
        //double sg_denom,renorm;
        //double mydf,dwsum;
        //double sigmasq,sigma,sigmapow;
        //int mmm;
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
        // could be made more compact, I think.
        if (renormalize) {
            mydf = double(frets.nel());
            if (mydf >= min_df) {
                xret(lll,2) = mydf;
                // put them in backwards!
                if (mydf >= 2) {
                    xret(lll,1) = frets.m_xx[1];
                    xret(lll,0) = frets.sd(renormalize,used_df); 
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
            mydf = double(frets.wsum());
            if (mydf >= min_df) {
                // put them in backwards!
                xret(lll,2) = mydf;
                if (mydf >= 2) {
                    xret(lll,1) = frets.m_xx[1];
                    xret(lll,0) = frets.sd(renormalize,used_df);
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
        if (renormalize) { 
            mydf = double(frets.nel());
        } else {
            mydf = double(frets.wsum());
        }
        if (mydf >= min_df) {
            // put them in backwards!
            if (mydf >= 3) {
                xret(lll,3) = mydf; 
                xret(lll,2) = frets.m_xx[1];
                xret(lll,1) = frets.sd(renormalize,used_df);
                xret(lll,0) = frets.skew();
            } else {
                xret(lll,3) = mydf; 
                if (mydf >= 1) {
                    xret(lll,2) = frets.m_xx[1];
                    if (mydf >= 2) {
                        xret(lll,1) = frets.sd(renormalize,used_df);
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
        if (renormalize) { 
            mydf = double(frets.nel());
        } else {
            mydf = double(frets.wsum());
        }

        if (mydf >= min_df) {
            // put them in backwards!
            if (mydf >= 4) {
                xret(lll,4) = mydf; 
                xret(lll,3) = frets.m_xx[1];
                xret(lll,2) = frets.sd(renormalize,used_df);
                // uhoh! renormalization!?
                xret(lll,1) = frets.skew();
                // uhoh! renormalization!?
                xret(lll,0) = frets.exkurt();
            } else {
                xret(lll,4) = mydf; 
                if (mydf >= 1) {
                    xret(lll,3) = frets.m_xx[1];
                    if (mydf >= 2) {
                        xret(lll,2) = frets.sd(renormalize,used_df);
                        if (mydf >= 3) {
                            // uhoh! renormalization!?
                            xret(lll,1) = frets.skew();
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
    //crap versions of skew4 and kurt5
    //if (retwhat==ret_skew4) { //FOLDUP


        //dwsum = double(frets.wsum());
        //if (renormalize) { 
            //mydf = double(frets.nel());
            //renorm = mydf / dwsum;
        //} else {
            //mydf = dwsum;
        //}

        //if (mydf >= min_df) {
            //sg_denom = mydf - used_df;
            //if (renormalize) { sg_denom /= renorm; }
            //sigmasq = frets.m_xx[2] / sg_denom;
            //sigma = sqrt(sigmasq);

            //// put them in backwards!
            //if (mydf >= 3) {
                //xret(lll,3) = mydf; 
                //xret(lll,2) = frets.m_xx[1];
                //xret(lll,1) = sigma;
                //xret(lll,0) = COMP_SKEW_TWO(frets.m_xx,dwsum);
            //} else {
                //xret(lll,3) = mydf; 
                //if (mydf >= 1) {
                    //xret(lll,2) = frets.m_xx[1];
                    //if (mydf >= 2) {
                        //xret(lll,1) = sigma;
                    //} else {
                        //xret(lll,1) = NAN;
                    //}
                //} else {
                    //xret(lll,2) = NAN;
                    //xret(lll,1) = NAN;
                //}
                //xret(lll,0) = NAN;
            //}
        //} else {
            //xret(lll,3) = NAN;
            //xret(lll,2) = NAN;
            //xret(lll,1) = NAN;
            //xret(lll,0) = NAN;
        //}
    //}//UNFOLD
    //if (retwhat==ret_exkurt5) { //FOLDUP
        ////double sg_denom,renorm;
        ////double mydf,dwsum;
        ////double sigmasq,sigma,sigmapow;
        ////int mmm;
        //dwsum = double(frets.wsum());
        //if (renormalize) { 
            //mydf = double(frets.nel());
            //renorm = mydf / dwsum;
        //} else {
            //mydf = dwsum;
        //}

        //if (mydf >= min_df) {
            //sg_denom = mydf - used_df;
            //if (renormalize) { sg_denom /= renorm; }
            //sigmasq = frets.m_xx[2] / sg_denom;
            //sigma = sqrt(sigmasq);

            //// put them in backwards!
            //if (mydf >= 4) {
                //xret(lll,4) = mydf; 
                //xret(lll,3) = frets.m_xx[1];
                //xret(lll,2) = sigma;
                //// uhoh! renormalization!
                //xret(lll,1) = COMP_SKEW_TWO(frets.m_xx,dwsum);
                //// uhoh! renormalization!
                //xret(lll,0) = (COMP_KURT_TWO(frets.m_xx,dwsum)) - 3.0;
            //} else {
                //xret(lll,4) = mydf; 
                //if (mydf >= 1) {
                    //xret(lll,3) = frets.m_xx[1];
                    //if (mydf >= 2) {
                        //xret(lll,2) = sigma;
                        //if (mydf >= 3) {
                            //// uhoh! renormalization!
                            //xret(lll,1) = COMP_SKEW_TWO(frets.m_xx,dwsum);
                        //} else {
                            //xret(lll,1) = NAN;
                        //}
                    //} else {
                        //xret(lll,2) = NAN;
                        //xret(lll,1) = NAN;
                    //}
                //} else {
                    //xret(lll,3) = NAN;
                    //xret(lll,2) = NAN;
                    //xret(lll,1) = NAN;
                //}
                //xret(lll,0) = NAN;
            //}
        //} else {
            //xret(lll,4) = NAN;
            //xret(lll,3) = NAN;
            //xret(lll,2) = NAN;
            //xret(lll,1) = NAN;
            //xret(lll,0) = NAN;
        //}
    //}//UNFOLD
    if (retwhat==ret_centmaxonly) { //FOLDUP
        mydf = double(frets.wsum());
        if ((mydf >= min_df) && (mydf >= ord)) {
            xret(lll,0) = frets.a_cent_mom(ord,renormalize,used_df);
        } else {
            xret(lll,0) = NAN;
        }
    } //UNFOLD
    if (retwhat==ret_centered) {//FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.centered(double(v[lll]));
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.centered(double(v[lll]));
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_scaled) {//FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.scaled(double(v[lll]),renormalize,used_df);
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.scaled(double(v[lll]),renormalize,used_df);
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_zscore) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.zscored(double(v[lll]),renormalize,used_df);
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.zscored(double(v[lll]),renormalize,used_df);
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_tstat) {//FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.nel()));
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = (frets.mean() / frets.sd(renormalize,used_df)) * sqrt(double(frets.wsum()));
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_sharpe) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.sharpe(renormalize,used_df); 
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.sharpe(renormalize,used_df); 
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_sharpese) { //FOLDUP
        //double skew,exkurt,sr;
        if (renormalize) {
            if (frets.nel() >= min_df) {
                // same
                skew = frets.skew();
                exkurt = frets.exkurt();
                sr = frets.sharpe(renormalize,used_df);
                xret(lll,0) = sr;
                // different
                mydf = double(frets.nel());

                xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / mydf);
            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                // same
                skew = frets.skew();
                exkurt = frets.exkurt();
                sr = frets.sharpe(renormalize,used_df);
                xret(lll,0) = sr;
                // different
                mydf = double(frets.wsum());

                xret(lll,1) = sqrt((1.0 + sr * (0.25 * (2.0 + exkurt) * sr - skew)) / mydf);

            } else {
                xret(lll,0) = NAN;
                xret(lll,1) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_stdev) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.sd(renormalize,used_df); 
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.sd(renormalize,used_df); 
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD
    if (retwhat==ret_skew) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.skew();
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.skew();
            } else {
                xret(lll,0) = NAN;
            }
        }
    } //UNFOLD
    if (retwhat==ret_exkurt) { //FOLDUP
        if (renormalize) {
            if (frets.nel() >= min_df) {
                xret(lll,0) = frets.exkurt();
            } else {
                xret(lll,0) = NAN;
            }
        } else {
            if (frets.wsum() >= min_df) {
                xret(lll,0) = frets.exkurt();
            } else {
                xret(lll,0) = NAN;
            }
        }
    }//UNFOLD

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
