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

#ifndef __DEF_REF_IMPL__
#define __DEF_REF_IMPL__

#include "common.h"
#include "welford.h"
#include "running.h"

#endif /* __DEF_REF_IMPL__ */

#include <Rcpp.h>
using namespace Rcpp;

// try to use c++11 ? 
// [[Rcpp::plugins(cpp11)]]

// reference implementations for speed;
// these should be the fastest possible versions of sd and running_sd

// return the sd
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
double ref_sd(NumericVector v) {
    double mu,sd,delta;
    double x;
    
    int numel=v.size();
    sd = 0.0;
    mu = v[0];
    for (int iii=2;iii <= numel;++iii) {
        x = v[iii-1];
        delta = x - mu;
        mu += delta / double(iii);
        sd += delta * (x - mu);
    }
    //NumericVector vret = NumericVector::create(sqrt(sd / (nel - 1)));
    //return vret;
    return sqrt(sd / (numel - 1));
}
// same, but use the Welford object. another comparison point.
//' @rdname firstmoments
//' @export
// [[Rcpp::export]]
double ref_sd_objecty(NumericVector v) {
    Welford<double,false,false,false> frets = Welford<double,false,false,false>(2);
    int numel=v.size();
    for (int iii=0;iii < numel;++iii) { frets.add_one(v[iii],1.0); }
    return frets.sd(false,1.0);
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd(NumericVector v,int window=1000) {
    double nel,mu,sd,delta;
    double diffmu,sumx,summu,prevmu;
    double prevx,nextx;
    double x;
    int jjj;
    
    int numel=v.size();
    NumericVector vret = NumericVector(numel);
    const int firstpart = MIN(numel,window);

    nel = 0.0;
    sd = 0.0;
    mu = 0.0;
    for (int iii=0;iii < firstpart;++iii) {
        ++nel;
        x = v[iii];
        delta = x - mu;
        mu += delta / nel;
        sd += delta * (x - mu);
        vret[iii] = sqrt(sd / (nel - 1));
    }
    if (firstpart < numel) {
        jjj = 0;
        for (int iii=firstpart;iii < numel;++iii) {
            //++nel;
            //x = v[iii];
            //delta = x - mu;
            //mu += delta / nel;
            //sd += delta * (x - mu);

            //--nel;
            //x = v[jjj];
            //delta = x - mu;
            //mu -= delta / nel;
            //sd -= delta * (x - mu);
            nextx = v[iii];
            prevx = v[jjj];
            diffmu = nextx - prevx;
            prevmu = mu;
            mu += (diffmu / nel);
            sd += diffmu*(nextx + prevx - (prevmu + mu));
            ++jjj;
            vret[iii] = sqrt(sd / (nel - 1));
        }
    }
    return vret;
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_narm(NumericVector v,int window=1000) {
    // simulate checking for na
    double nel,mu,sd,delta;
    double x;
    int jjj;
    double addx,remx,diffmu,prevmu;
    
    int numel=v.size();
    NumericVector vret = NumericVector(numel);
    const int firstpart = MIN(numel,window);

    nel = 0.0;
    sd = 0.0;
    mu = 0.0;
    for (int iii=0;iii < firstpart;++iii) {
        x = v[iii];
        if (!ISNAN(x)) {
            ++nel;
            delta = x - mu;
            mu += delta / nel;
            sd += delta * (x - mu);
        }
        vret[iii] = sqrt(sd / (nel - 1));
    }
    if (firstpart < numel) {
        jjj = 0;
        for (int iii=firstpart;iii < numel;++iii) {
            addx = v[iii];
            remx = v[jjj];
            ++jjj;
            if (!ISNAN(addx)) {
                if (!ISNAN(remx)) {
                    diffmu = addx - remx;
                    prevmu = mu;
                    mu += (diffmu / nel);
                    sd += diffmu*(addx + remx - (prevmu + mu));
                } else {
                    ++nel;
                    delta = x - mu;
                    mu += delta / nel;
                    sd += delta * (x - mu);
                }
            } else if (!ISNAN(remx)) {
                --nel;
                delta = x - mu;
                mu -= delta / nel;
                sd -= delta * (x - mu);
            }
            vret[iii] = sqrt(sd / (nel - 1));
        }
    }
    return vret;
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_onecheck(NumericVector v,int window=1000,bool na_rm=false) {
    if (na_rm) {
        return ref_running_sd_narm(v,window);
    }
    return ref_running_sd(v,window);
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_intnel(NumericVector v,int window=1000) {
    // integer number of elements and double conversion? no time.
    double nel,mu,sd,delta;
    double x;
    int jjj;
    int inel;
    double nextx,prevx,diffmu,prevmu;
    
    int numel=v.size();
    NumericVector vret = NumericVector(numel);
    const int firstpart = MIN(numel,window);

    inel = 0;
    nel = 0.0;
    sd = 0.0;
    mu = 0.0;
    for (int iii=0;iii < firstpart;++iii) {
        ++inel;
        nel = double(inel);
        x = v[iii];
        delta = x - mu;
        mu += delta / nel;
        sd += delta * (x - mu);
        vret[iii] = sqrt(sd / (nel - 1));
    }
    if (firstpart < numel) {
        jjj = 0;
        for (int iii=firstpart;iii < numel;++iii) {
            nextx = v[iii];
            prevx = v[jjj];
            jjj++;
            diffmu = nextx - prevx;
            prevmu = mu;
            mu += (diffmu / nel);
            sd += diffmu*(nextx + prevx - (prevmu + mu));
            vret[iii] = sqrt(sd / (nel - 1));
        }
    }
    return vret;
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_objecty(NumericVector v,int window=1000) {
    // integer number of elements and double conversion? no time.
    //
    Welford<double,false,false,false> frets = Welford<double,false,false,false>(2);
    double x;
    double nextx,prevx;
    int jjj;
    
    int numel=v.size();
    NumericVector vret = NumericVector(numel);
    const int firstpart = MIN(numel,window);

    for (int iii=0;iii < firstpart;++iii) {
        x = v[iii];
        frets.add_one(x,1.0);
        vret[iii] = frets.sd(false,1.0);
    }
    if (firstpart < numel) {
        jjj = 0;
        for (int iii=firstpart;iii < numel;++iii) {
            nextx = v[iii];
            prevx = v[jjj];
            ++jjj;
            frets.swap_one(nextx,1.0,prevx,1.0);
            vret[iii] = frets.sd(false,1.0);
        }
    }
    return vret;
}

// check if the currying layers are a problem; apparently not, the slowness is in runQMCurryZero, right?
//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_fooz(NumericVector v,int window=1000) {
    NumericVector dummy_wts;
    return runQM<NumericVector,ret_stdev,NumericVector,double,false,false,false,false>(v,dummy_wts,2,window,10000,0,0,0.0,FALSE,FALSE);
}

//' @export
//' @rdname runningmoments
// [[Rcpp::export]]
NumericVector ref_running_sd_barz(NumericVector v,int window=1000) {
    double nextv, prevv;
    double nextw;
    NumericVector dummy_wts;
    int ord=2;
    bool na_rm=false;
    bool check_wts=false;

    // only for retwhat==ret_sharpese, but cannot define outside its scope.
    // no bigs.
    double sigma,skew,exkurt,sr;

    int iii,jjj,lll,mmm,ppp,qqq,tr_iii,tr_jjj;
    int numel = v.size();

    // preallocated with zeros; should
    // probably be NA?
    int ncols=1;
    NumericVector xret(numel);

    // as an invariant, we will start the computation
    // with vret, which is initialized as the summed
    // means on [jjj,iii]
    tr_iii = - 1;
    tr_jjj = - window;

    Welford<double,false,false,false> frets = Welford<double,false,false,false>(2);
    Welford<double,false,false,true> trets = Welford<double,false,false,true>(2);
    if (na_rm) {
        trets = quasiWeightedThing<NumericVector,NumericVector,double,false,false,true>(v,dummy_wts,ord,
                                                                                        0,       //bottom
                                                                                        0,     //top
                                                                                        check_wts);


        // now run through lll index//FOLDUP
        for (lll=0;lll < numel;++lll) {
            tr_iii++;
            if ((tr_iii < numel) && (tr_iii >= 0)) {
                // add on nextv:
                nextv = double(v[tr_iii]);
                trets.add_one(nextv,1);
            }
            // remove prevv:
            if ((tr_jjj < numel) && (tr_jjj >= 0)) {
                prevv = double(v[tr_jjj]);
                trets.rem_one(prevv,1);
            }
            tr_jjj++;

            // fill in the value in the output.
            // 2FIX: give access to v, not v[lll]...
            xret[lll] = trets.sd(false,1.0);
        }//UNFOLD
    } else {
        frets = quasiWeightedThing<NumericVector,NumericVector,double,false,false,false>(v,dummy_wts,ord,
                                                                                            0,       //bottom
                                                                                            0,     //top
                                                                                            check_wts);
        // now run through lll index//FOLDUP
        for (lll=0;lll < numel;++lll) {
            tr_iii++;
            if ((tr_iii < numel) && (tr_iii >= 0)) {
                // add on nextv:
                nextv = double(v[tr_iii]);
                frets.add_one(nextv,1);
            }
            // remove prevv:
            if ((tr_jjj < numel) && (tr_jjj >= 0)) {
                prevv = double(v[tr_jjj]);
                frets.rem_one(prevv,1);
            }
            tr_jjj++;

            // fill in the value in the output.
            // 2FIX: give access to v, not v[lll]...
            xret[lll] = frets.sd(false,1.0);
        }//UNFOLD
    }
    return xret;
}

//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
