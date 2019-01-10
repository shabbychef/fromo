# Copyright 2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of fromo.
#
# fromo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fromo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with fromo.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2016.03.28
# Copyright: Steven E. Pav, 2016-2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)

# slow version of t_running; requires a helper function.
	slow_op <- function(v,func,outsize=1,missfill=NA,
											time=NULL,time_deltas=NULL,window=Inf,wts=NULL,lb_time=NULL,
											na_rm=FALSE,min_df=0,lookahead=0,variable_win=FALSE,wts_as_delta=TRUE,...) {
		if (is.null(time)) {
			if (is.null(time_deltas) && !is.null(wts) && wts_as_delta) {
				time_deltas <- wts
			} else {
				stop('bad input')
			}
			time <- cumsum(time_deltas)
		}
		if (is.null(lb_time)) {
			lb_time <- time
		}
		lb_time <- lb_time + lookahead
		if (variable_win) {
			tprev <- c(-Inf,lb_time[1:(length(lb_time)-1)])
		} else {
			tprev <- lb_time - window
		}
		# fix weights.
		sapply(seq_along(lb_time),
			function(idx) {
				tf <- lb_time[idx]
				t0 <- tprev[idx]
				takeus <- (t0 < time) & (time <= tf)
				if (na_rm) {
					takeus <- takeus & !is.na(v)
					if (!is.null(wts)) {
						takeus <- takeus & !is.na(wts)
					}
				}
				if (any(takeus)) {
					vsub <- v[takeus]
					if (is.null(wts)) {
						retv <- func(vsub,...)
					} else {
						subwts <- wts[takeus]
						retv <- func(vsub,wts=subwts,...)
					}
				} else {
					retv <- rep(missfill,outsize)
				}
				retv
			})
	}

	slow_t_running_sum <- function(v,...) {
		func <- function(v,wts=NULL,...) { 
			if (is.null(wts)) {
				return(prod(sd3(v,...)[c(2:3)])) 
			}
			return(sum(v*wts))
		}
		as.numeric(slow_op(v=v,func=func,missfill=0,...))
	}
	slow_t_running_mean <- function(v,...) {
		func <- function(v,...) { sd3(v,...)[2] }
		as.numeric(slow_op(v=v,func=func,...))
	}
	slow_t_running_sd <- function(v,...) {
		func <- function(v,...) { sd3(v,...)[1] }
		matrix(slow_op(v=v,func=func,...),ncol=1)
	}
	slow_t_running_skew <- function(v,...) {
		func <- function(v,...) { 
			return(skew4(v,...)[1]) 
		}
		matrix(slow_op(v=v,func=func,...),ncol=1)
	}
	slow_t_running_kurt <- function(v,...) {
		func <- function(v,...) { kurt5(v,...)[1] }
		matrix(slow_op(v=v,func=func,...),ncol=1)
	}

	slow_t_running_sd3 <- function(v,...) {
		func <- function(v,...) { sd3(v,...) }
		t(slow_op(v=v,func=func,outsize=3,...))
	}
	slow_t_running_skew4 <- function(v,...) {
		func <- function(v,...) { 
			if (length(v) > 2) {
				return(skew4(v,...))
			} else {
				return(c(NA,sd3(v,...)))
			}
		}
		t(slow_op(v=v,func=func,outsize=4,...))
	}
	slow_t_running_kurt5 <- function(v,...) {
		func <- function(v,...) { 
			if (length(v) > 3) {
				return(kurt5(v,...))
			} else if (length(v) > 2) {
				return(c(NA,skew4(v,...)))
			} else {
				return(c(NA,NA,sd3(v,...)))
			}
		}
		t(slow_op(v=v,func=func,outsize=5,...))
	}

	reference_sd <- function(x,wts=NULL,na_rm=FALSE,normalize_wts=FALSE,min_df=0,used_df=1) {
		if (na_rm) {
			isok <- !is.na(x)
			if (!is.null(wts)) {
				isok <- isok & !is.na(wts) & wts >= 0
			}
			x <- x[isok]
			if (!is.null(wts)) {
				wts <- wts[isok]
			}
		}
		if (length(x) < min_df) {
			return(NA)
		}
		if (!is.null(wts)) {
			wsum <- sum(wts)
			mu <- sum(x*wts) / wsum
			deno <- wsum - used_df * ifelse(normalize_wts,wsum / length(x),1)
			vv <- sum(wts * (x - mu)^2) / deno
		} else {
			wsum <- length(x)
			mu <- sum(x) / wsum
			vv <- sum((x - mu)^2) / (wsum - used_df)
		}
		return(sqrt(vv))
	}
	# not quite the same as slow_t_running_sd above
	reference_t_running_sd <- function(v,...) {
		matrix(slow_op(v=v,func=reference_sd,...),ncol=1)
	}
#UNFOLD

context("first moments")#FOLDUP
test_that("sd, skew, kurt are correct",{#FOLDUP
	set.char.seed("c4007dba-2010-481e-abe5-f07d3ce94eb4")
	x <- rnorm(1000)

	expect_error(sid <- sd3(x),NA)
	expect_error(ske <- skew4(x),NA)
	expect_error(krt <- kurt5(x),NA)

	expect_equal(length(sid),3)
	expect_equal(length(ske),4)
	expect_equal(length(krt),5)

	# compare computations to gold standard
	# length
	expect_equal(sid[3],length(x))
	expect_equal(sid[3],ske[4])
	expect_equal(sid[3],krt[5])

	# mean
	expect_equal(sid[2],mean(x),tolerance=1e-9)
	expect_equal(sid[2],ske[3],tolerance=1e-9)
	expect_equal(sid[2],krt[4],tolerance=1e-9)

	# standard dev
	expect_equal(sid[1],ske[2],tolerance=1e-9)
	expect_equal(sid[1],krt[3],tolerance=1e-9)
	expect_equal(sid[1],sd(x),tolerance=1e-9)

	# skew
	expect_equal(ske[1],krt[2],tolerance=1e-9)

	if (require(moments)) {
		na_rm <- TRUE
		dumb_count  <- sum(sign(abs(x)+1),na.rm=na_rm) 
		dumb_mean   <- mean(x,na.rm=na_rm) 
		dumb_sd     <- sd(x,na.rm=na_rm) 
		dumb_skew   <- moments::skewness(x,na.rm=na_rm) 
		dumb_exkurt <- moments::kurtosis(x,na.rm=na_rm) - 3.0 

		dumb_cmom2 <- moments::moment(x,central=TRUE,na.rm=na_rm,order=2) 
		dumb_cmom3 <- moments::moment(x,central=TRUE,na.rm=na_rm,order=3)
		dumb_cmom4 <- moments::moment(x,central=TRUE,na.rm=na_rm,order=4)
		dumb_cmom5 <- moments::moment(x,central=TRUE,na.rm=na_rm,order=5)
		dumb_cmom6 <- moments::moment(x,central=TRUE,na.rm=na_rm,order=6)

		# skew
		expect_equal(ske[1],dumb_skew,tolerance=1e-9)
		# kurtosis
		expect_equal(krt[1],dumb_exkurt,tolerance=1e-9)

		# oops. problems with centered moments in terms of the used_df; need a
		# better test...
		cmoms <- cent_moments(x,max_order=6,used_df=0)
		dumbv <- c(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count)
		expect_equal(max(abs(cmoms-dumbv)),0,tolerance=1e-9)

		if (require(PDQutils)) {
			cumuls <- cent_cumulants(x,max_order=length(cmoms)-1)
			dumbv0 <- c(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count)
			dumbv1 <- PDQutils::moment2cumulant(c(0,rev(dumbv0)[3:length(dumbv0)]))
			dumbv <- c(rev(dumbv1[2:length(dumbv1)]),dumb_mean,dumb_count)
			
			expect_equal(max(abs(cumuls-dumbv)),0,tolerance=1e-12)
		}
	}
	if (require(e1071)) {
		dumb_skew   <- e1071::skewness(x,type=3)
		equiv_skew  <- ske[1] * ((ske[4]-1)/(ske[4]))^(3/2)
		expect_equal(dumb_skew,equiv_skew,tolerance=1e-12)
	}

	# 2FIX: add cent_moments and std_moments
	# 2FIX: check NA

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("unit weighted sd, skew, kurt are correct",{#FOLDUP
	set.char.seed("b652ccd2-478b-44d4-90e2-2ca2bad99d25")
	x <- rnorm(1000)
	ones <- rep(1,length(x))

	expect_equal(sd3(x),sd3(x,wts=ones),tolerance=1e-9)
	expect_equal(skew4(x),skew4(x,wts=ones),tolerance=1e-9)
	expect_equal(kurt5(x),kurt5(x,wts=ones),tolerance=1e-9)
	# 2FIX: probably normalize_wts=FALSE should be the default???? for speed?
	expect_equal(running_sd(x),running_sd(x,wts=ones,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_skew(x),running_skew(x,wts=ones,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_kurt(x),running_kurt(x,wts=ones,normalize_wts=TRUE),tolerance=1e-9)

	expect_equal(running_sd(x),running_sd(x,wts=ones,normalize_wts=FALSE),tolerance=1e-9)
	expect_equal(running_skew(x),running_skew(x,wts=ones,normalize_wts=FALSE),tolerance=1e-9)
	expect_equal(running_kurt(x),running_kurt(x,wts=ones,normalize_wts=FALSE),tolerance=1e-9)

	# 2FIX: add more.
	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("normalize weights works",{#FOLDUP
	set.char.seed("2694ae87-62d4-4154-9c32-864f9a6e648d")
	x <- rnorm(25)
	wts <- runif(length(x))

	expect_equal(sd3(x,wts=wts,normalize_wts=TRUE),
							 sd3(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(skew4(x,wts=wts,normalize_wts=TRUE),
							 skew4(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(kurt5(x,wts=wts,normalize_wts=TRUE),
							 kurt5(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)

	expect_equal(running_sd(x,wts=wts,normalize_wts=TRUE),
							 running_sd(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_skew(x,wts=wts,normalize_wts=TRUE),
							 running_skew(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_kurt(x,wts=wts,normalize_wts=TRUE),
							 running_kurt(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_sd3(x,wts=wts,normalize_wts=TRUE),
							 running_sd3(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_skew4(x,wts=wts,normalize_wts=TRUE),
							 running_skew4(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_kurt5(x,wts=wts,normalize_wts=TRUE),
							 running_kurt5(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_sharpe(x,wts=wts,normalize_wts=TRUE),
							 running_sharpe(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_sharpe(x,wts=wts,normalize_wts=TRUE,compute_se=TRUE),
							 running_sharpe(x,wts=2*wts,normalize_wts=TRUE,compute_se=TRUE),tolerance=1e-9)
	expect_equal(running_centered(x,wts=wts,normalize_wts=TRUE),
							 running_centered(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_apx_median(x,wts=wts,normalize_wts=TRUE),
							 running_apx_median(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_tstat(x,wts=wts,normalize_wts=TRUE),
							 running_tstat(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_zscored(x,wts=wts,normalize_wts=TRUE),
							 running_zscored(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_scaled(x,wts=wts,normalize_wts=TRUE),
							 running_scaled(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	expect_equal(running_apx_quantiles(x,p=ptiles,wts=wts,normalize_wts=TRUE),
							 running_apx_quantiles(x,p=ptiles,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_cent_moments(x,wts=wts,normalize_wts=TRUE),
							 running_cent_moments(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_std_moments(x,wts=wts,normalize_wts=TRUE),
							 running_std_moments(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
	expect_equal(running_cumulants(x,wts=wts,normalize_wts=TRUE),
							 running_cumulants(x,wts=2*wts,normalize_wts=TRUE),tolerance=1e-9)
})#UNFOLD
test_that("weight scaling what you expect",{#FOLDUP
	set.char.seed("efaa75ac-bb9e-4e4a-a375-7028f099366e")
	x <- rnorm(50)
	wts <- runif(length(x))

	expect_error(sid_1 <- sd3(x,wts=wts,normalize_wts=FALSE,sg_df=0),NA)
	expect_error(ske_1 <- skew4(x,wts=wts,normalize_wts=FALSE,sg_df=0),NA)
	expect_error(krt_1 <- kurt5(x,wts=wts,normalize_wts=FALSE,sg_df=0),NA)

	expect_error(sid_2 <- sd3(x,wts=2*wts,normalize_wts=FALSE,sg_df=0),NA)
	expect_error(ske_2 <- skew4(x,wts=2*wts,normalize_wts=FALSE,sg_df=0),NA)
	expect_error(krt_2 <- kurt5(x,wts=2*wts,normalize_wts=FALSE,sg_df=0),NA)

	expect_equal(sid_1 * c(1,1,2),sid_2,tolerance=1e-9)
	expect_equal(ske_1 * c(1,1,1,2),ske_2,tolerance=1e-9)
	expect_equal(krt_1 * c(1,1,1,1,2),krt_2,tolerance=1e-9)
})#UNFOLD
test_that("weighted sd, skew, kurt are correct",{#FOLDUP
	set.char.seed("4e17d837-69c1-41d1-906f-c82224d7ce41")
	x <- rnorm(1000)
	wts <- runif(length(x))

	expect_error(sid <- sd3(x,wts=wts,normalize_wts=TRUE),NA)
	expect_error(ske <- skew4(x,wts=wts,normalize_wts=TRUE),NA)
	expect_error(krt <- kurt5(x,wts=wts,normalize_wts=TRUE),NA)
	# 2FIX: add more here to check correctness ... 

	expect_equal(length(sid),3)
	expect_equal(length(ske),4)
	expect_equal(length(krt),5)

	# compare computations to gold standard
	# length
	expect_equal(sid[3],length(x))
	expect_equal(sid[3],ske[4])
	expect_equal(sid[3],krt[5])

	# mean
	expect_equal(sid[2],weighted.mean(x,w=wts),tolerance=1e-9)
	expect_equal(sid[2],ske[3],tolerance=1e-9)
	expect_equal(sid[2],krt[4],tolerance=1e-9)

	# standard dev
	expect_equal(sid[1],ske[2],tolerance=1e-9)
	expect_equal(sid[1],krt[3],tolerance=1e-9)
	wsd <- sqrt(sum(((x - weighted.mean(x,w=wts))^2) * (wts / mean(wts))) / (length(x) - 1))
	# 2FIX!!!
	expect_equal(sid[1],wsd,tolerance=1e-9)

	# skew
	expect_equal(ske[1],krt[2],tolerance=1e-9)
 
	na_rm <- TRUE
	dumb_count  <- length(x)
	dumb_mean   <- weighted.mean(x,w=wts)
	dumb_sd     <- sqrt(sum(((x - weighted.mean(x,w=wts))^2) * (wts / mean(wts))) / (length(x) - 1))

	wcmom <- function(vec,wts,ord) {
		wz <- wts / mean(wts)
		mean(wz * ((x - weighted.mean(x,w=wz))^ord))
	}
	dumb_wcmom2  <- wcmom(x,wts,2)
	dumb_wcmom3  <- wcmom(x,wts,3)
	dumb_wcmom4  <- wcmom(x,wts,4)
	dumb_wcmom5  <- wcmom(x,wts,5)
	dumb_wcmom6  <- wcmom(x,wts,6)

	cmoms <- cent_moments(x,wts=wts,max_order=6,used_df=0,normalize_wts=TRUE)
	dumbv <- c(dumb_wcmom6,dumb_wcmom5,dumb_wcmom4,dumb_wcmom3,dumb_wcmom2,dumb_mean,dumb_count)
	expect_equal(cmoms,dumbv,tolerance=1e-9)

	dumb_skew <- dumb_wcmom3 / (dumb_wcmom2^(3/2))
	dumb_exkurt <- (dumb_wcmom4 / (dumb_wcmom2^(2))) - 3

	# skew
	expect_equal(ske[1],dumb_skew,tolerance=1e-9)
	# kurtosis
	expect_equal(krt[1],dumb_exkurt,tolerance=1e-9)
})#UNFOLD
#UNFOLD

tomat <- function(cbound) {
	dumbv <- as.matrix(cbound)
	attr(dumbv,'dimnames') <- NULL 
	dumbv
}

context("running ops are correct")
test_that("running ops are correct",{#FOLDUP
	skip_on_cran()
	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	set.char.seed("7ffe0035-2d0c-4586-a1a5-6321c7cf8694")
	for (xlen in c(20,100)) {
		for (xmu in c(1e3,1e6)) {
			toler <- xmu ^ (1/3)
			x <- rnorm(xlen,mean=xmu)
			for (window in c(15,50,Inf)) {
				for (restart_period in c(20,1000)) {
					for (na_rm in c(FALSE,TRUE)) {
						dumb_count <- sapply(seq_along(x),function(iii) { sum(sign(abs(x[max(1,iii-window+1):iii])+1),na.rm=na_rm) },simplify=TRUE)
						dumb_sum <- sapply(seq_along(x),function(iii) { sum(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
						dumb_mean <- sapply(seq_along(x),function(iii) { mean(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
						dumb_sd <- sapply(seq_along(x),function(iii) { sd(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
						dumb_skew <- sapply(seq_along(x),function(iii) { moments::skewness(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
						dumb_exkurt <- sapply(seq_along(x),function(iii) { moments::kurtosis(x[max(1,iii-window+1):iii],na.rm=na_rm) - 3.0 },simplify=TRUE)

						dumb_cmom2 <- sapply(seq_along(x),function(iii) { moments::moment(x[max(1,iii-window+1):iii],central=TRUE,na.rm=na_rm,order=2) },simplify=TRUE)
						dumb_cmom3 <- sapply(seq_along(x),function(iii) { moments::moment(x[max(1,iii-window+1):iii],central=TRUE,na.rm=na_rm,order=3) },simplify=TRUE)
						dumb_cmom4 <- sapply(seq_along(x),function(iii) { moments::moment(x[max(1,iii-window+1):iii],central=TRUE,na.rm=na_rm,order=4) },simplify=TRUE)
						dumb_cmom5 <- sapply(seq_along(x),function(iii) { moments::moment(x[max(1,iii-window+1):iii],central=TRUE,na.rm=na_rm,order=5) },simplify=TRUE)
						dumb_cmom6 <- sapply(seq_along(x),function(iii) { moments::moment(x[max(1,iii-window+1):iii],central=TRUE,na.rm=na_rm,order=6) },simplify=TRUE)

						# SD
						expect_error(fastv <- running_sd(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(dumb_sd)
						expect_equal(dumbv[2:xlen],fastv[2:xlen],tolerance=1e-7 * toler)

						expect_error(fastv <- running_sd3(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(cbind(dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[2:xlen,],fastv[2:xlen,],tolerance=1e-7 * toler)

						# skew
						expect_error(fastv <- running_skew(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(dumb_skew)
						expect_equal(dumbv[3:xlen],fastv[3:xlen],tolerance=1e-6 * toler)

						expect_error(fastv <- running_skew4(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(cbind(dumb_skew,dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[3:xlen,],fastv[3:xlen,],tolerance=1e-7 * toler)

						# excess kurtosis
						expect_error(fastv <- running_kurt(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(dumb_exkurt)
						expect_equal(dumbv[4:xlen],fastv[4:xlen],tolerance=1e-6 * toler)

						expect_error(fastv <- running_kurt5(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(cbind(dumb_exkurt,dumb_skew,dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[4:xlen,],fastv[4:xlen,],tolerance=1e-6 * toler)

						# higher order moments

						expect_error(fastv <- running_cent_moments(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(cbind(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count))
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-6 * toler)

						expect_error(fastv <- running_cent_moments(x,window=window,max_order=6L,max_order_only=TRUE,used_df=0L,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(dumb_cmom6)
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-7 * toler)

						expect_error(fastv <- running_std_moments(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- tomat(cbind(dumb_cmom6 / (dumb_cmom2^3),dumb_cmom5 / (dumb_cmom2^2.5),dumb_cmom4 / (dumb_cmom2^2.0),dumb_cmom3 / (dumb_cmom2^1.5),sqrt(dumb_cmom2),dumb_mean,dumb_count))
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-7 * toler)

						# running sum and mean
						expect_error(fastv <- running_sum(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- dumb_sum
						expect_equal(dumbv[2:xlen],fastv[2:xlen],tolerance=1e-7 * toler)

						expect_error(fastv <- running_mean(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
						dumbv <- dumb_mean
						expect_equal(dumbv[2:xlen],fastv[2:xlen],tolerance=1e-7 * toler)

						if (require(PDQutils)) {
							# cumulants
							expect_error(fastv <- running_cumulants(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm),NA)
							pre_dumbv <- cbind(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count)
							dumbv <- t(sapply(seq_along(x),function(iii) { 
															rv <- rev(PDQutils::moment2cumulant(c(0,rev(pre_dumbv[iii,1:(ncol(pre_dumbv)-2)]))))
															rv <- rv[-length(rv)]
															c(rv,pre_dumbv[iii,ncol(pre_dumbv) + (-1:0)])
								},simplify='matrix'))
							expect_equal(max(abs(dumbv[6:xlen,] - fastv[6:xlen,])),0,tolerance=1e-8 * toler)

							# quantiles
							expect_error(fastv <- running_apx_quantiles(x,ptiles,max_order=ncol(dumbv)-1,used_df=0L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
							dumbq <- t(sapply(seq_along(x),function(iii) { 
								PDQutils::qapx_cf(ptiles,raw.cumulants=rev(dumbv[iii,1:(ncol(dumbv)-1)]))
							}, simplify=TRUE))
							expect_equal(max(abs(dumbq[8:xlen,] - fastv[8:xlen,])),0,tolerance=1e-8 * toler)
						}
					}
				}
			}
		}
	}
})#UNFOLD
test_that("running adjustments are correct",{#FOLDUP
	skip_on_cran()

	set.char.seed("967d2149-fbff-4d82-b227-ca3e1034bddb")
	for (xlen in c(20,100)) {
		x <- rnorm(xlen)
		for (window in c(5,50,Inf)) {
			for (restart_period in c(10,1000)) {
				for (na_rm in c(FALSE,TRUE)) {
					dumb_count <- sapply(seq_along(x),function(iii) { sum(sign(abs(x[max(1,iii-window+1):iii])+1),na.rm=na_rm) },simplify=TRUE)
					dumb_mean <- sapply(seq_along(x),function(iii) { mean(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
					dumb_sd <- sapply(seq_along(x),function(iii) { sd(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)

					expect_error(fastv <- running_centered(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
					# the dumb value:
					dumbv <- x - dumb_mean;
					expect_equal(max(abs(dumbv - fastv)),0,tolerance=1e-12)

					expect_error(fastv <- running_scaled(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
					# the dumb value:
					dumbv <- x / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					expect_error(fastv <- running_zscored(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
					# the dumb value:
					dumbv <- (x - dumb_mean) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					expect_error(fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
					# the dumb value:
					dumbv <- dumb_mean / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					expect_error(fastv <- running_tstat(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
					# the dumb value:
					dumbv <- (dumb_mean * sqrt(dumb_count)) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					expect_error(fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm,compute_se=TRUE),NA)
					# the dumb value:
					dumb_sr <- dumb_mean / dumb_sd
					expect_equal(max(abs(dumb_sr[2:length(x)] - fastv[2:length(x),1])),0,tolerance=1e-12)

					if (require(moments)) {
						dumb_skew <- sapply(seq_along(x),function(iii) { moments::skewness(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
						dumb_exkurt <- sapply(seq_along(x),function(iii) { moments::kurtosis(x[max(1,iii-window+1):iii],na.rm=na_rm) - 3.0 },simplify=TRUE)
						dumb_merse <- sqrt((1 + 0.25 * (2+dumb_exkurt) * dumb_sr^2 - dumb_skew * dumb_sr) / dumb_count)
						expect_equal(max(abs(dumb_merse[5:length(x)] - fastv[5:length(x),2])),0,tolerance=1e-9)
					}
				}
			}
		}
	}
})#UNFOLD


context("weighted running ops are correct")
test_that("running weights work correctly",{#FOLDUP
	skip_on_cran()

	set.char.seed("b82d252c-681b-4b98-9bb3-ffd17feeb4a1")
	na_rm <- FALSE

	restart_period <- 1000
	for (xlen in c(20,50)) {
		x <- rnorm(xlen)
		for (wts in list(rep(1L,xlen), runif(xlen,min=2,max=7))) { 
			for (window in c(5,30,Inf)) { # FOLDUP

				# 2FIX: add to this!
				slow_count <- sapply(seq_along(x),function(iii) { sum(sign(abs(x[max(1,iii-window+1):iii])+1),na.rm=na_rm) },simplify=TRUE)
				slow_sumwt <- sapply(seq_along(x),function(iii) { sum(wts[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
				slow_mean <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * x[mydx],na.rm=na_rm) / slow_sumwt[iii]
							 },simplify=TRUE)
				slow_var <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * (x[mydx] - slow_mean[iii])^2,na.rm=na_rm) / (slow_sumwt[iii] - 1)
							 },simplify=TRUE)
				slow_sd <- sqrt(slow_var)
				# the normalize version;
				slow_nvar <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															(slow_count[iii]/slow_sumwt[iii]) * sum(mywts * (x[mydx] - slow_mean[iii])^2,na.rm=na_rm) / (slow_count[iii] - 1)
							 },simplify=TRUE)
				slow_nsd <- sqrt(slow_nvar)

				slow_cent3 <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * (x[mydx] - slow_mean[iii])^3,na.rm=na_rm) / (slow_sumwt[iii])
							 },simplify=TRUE)
				slow_cent4 <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * (x[mydx] - slow_mean[iii])^4,na.rm=na_rm) / (slow_sumwt[iii])
							 },simplify=TRUE)

				expect_error(fastv <- running_mean(x,wts=wts,min_df=0,window=window,na_rm=na_rm),NA)
				expect_equal(fastv,slow_mean,tolerance=1e-8)

				expect_error(fastv <- running_centered(x,wts=wts,window=window,restart_period=restart_period,na_rm=na_rm),NA)
				slowv <- x - slow_mean;
				expect_equal(as.numeric(fastv),slowv,tolerance=1e-8)

				for (nw in c(TRUE,FALSE)) {
					if (nw) {
						use_sd <- slow_nsd
						use_df <- slow_count
					} else {
						use_sd <- slow_sd
						use_df <- slow_sumwt
					}

					expect_error(fastv <- running_sd(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(as.numeric(fastv),use_sd,tolerance=1e-8)

					expect_error(fastv <- running_sd3(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(as.numeric(fastv[,1]),use_sd,tolerance=1e-8)
					expect_equal(as.numeric(fastv[,2]),slow_mean,tolerance=1e-8)
					expect_equal(as.numeric(fastv[,3]),use_df,tolerance=1e-8)

					expect_error(fastv <- running_scaled(x,wts=wts,window=window,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- x / use_sd
					expect_equal(as.numeric(fastv[2:length(x)]),slowv[2:length(x)],tolerance=1e-8)

					expect_error(fastv <- running_zscored(x,wts=wts,window=window,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- (x - slow_mean) / use_sd
					expect_equal(slowv[2:length(x)],fastv[2:length(x)],tolerance=1e-12)

					expect_error(fastv <- running_sharpe(x,wts=wts,window=window,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- slow_mean / use_sd
					expect_equal(slowv[2:length(x)],fastv[2:length(x)],tolerance=1e-12)

					expect_error(fastv <- running_tstat(x,wts=wts,window=window,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- (slow_mean * sqrt(use_df)) / use_sd
					expect_equal(slowv[2:length(x)],fastv[2:length(x)],tolerance=1e-12)

					expect_error(fastv <- running_cent_moments(x,wts=wts,window=window,max_order=3L,max_order_only=TRUE,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- slow_cent3 
					expect_equal(slowv[3:length(x)],fastv[3:length(x)],tolerance=1e-12)

					expect_error(fastv <- running_cent_moments(x,wts=wts,window=window,max_order=4L,max_order_only=TRUE,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					slowv <- slow_cent4 
					expect_equal(slowv[4:length(x)],fastv[4:length(x)],tolerance=1e-12)
				}
			}# UNFOLD
		}
	}
})#UNFOLD

context("t_running for trivial case")
test_that("vs running ops",{#FOLDUP
	skip_on_cran()

	set.char.seed("712463ec-f266-4de7-89d2-ce3c824327b0")
	na_rm <- FALSE
	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (xlen in c(20,50)) {
		x <- rnorm(xlen)
		times <- seq_along(x)
		for (wts in list(NULL,rep(1L,xlen), runif(xlen,min=1.2,max=3.5))) {
			# 2FIX? Inf window?
			for (window in c(5,30,Inf)) { # FOLDUP
				# to avoid roundoff issues on double times.
				t_window <- window - 0.1

				expect_error(box <- running_sum(x,wts=wts,window=window,na_rm=na_rm),NA)
				expect_error(tbox <- t_running_sum(x,time=times,wts=wts,window=t_window,na_rm=na_rm),NA)
				expect_equal(box,tbox,tolerance=1e-8)

				expect_error(box <- running_mean(x,wts=wts,min_df=0,window=window,na_rm=na_rm),NA)
				expect_error(tbox <- t_running_mean(x,time=times,wts=wts,min_df=0,window=t_window,na_rm=na_rm),NA)
				expect_equal(box,tbox,tolerance=1e-8)

				for (nw in c(TRUE,FALSE)) { 
					expect_error(box <- running_sd(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_error(tbox <- t_running_sd(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_skew(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_skew(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_kurt(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_kurt(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_sd3(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_sd3(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)
					
					expect_error(box <- running_skew4(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_skew4(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_kurt5(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_kurt5(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_centered(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_centered(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_scaled(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_scaled(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_zscored(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_zscored(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_tstat(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_tstat(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					for (cse in c(TRUE,FALSE)) {
						expect_error(box <- running_sharpe(x,wts=wts,window=window,na_rm=na_rm,compute_se=cse,normalize_wts=nw),NA)
						# the 0.1 is to avoid roundoff issues on the double times.
						expect_error(tbox <- t_running_sharpe(x,time=times,wts=wts,window=t_window,na_rm=na_rm,compute_se=cse,normalize_wts=nw),NA)
						expect_equal(box,tbox,tolerance=1e-8)
					}

					expect_error(box <- running_apx_median(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_apx_median(x,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)

					expect_error(box <- running_apx_quantiles(x,ptiles,max_order=3,wts=wts,window=window,na_rm=na_rm,normalize_wts=nw),NA)
					# the 0.1 is to avoid roundoff issues on the double times.
					expect_error(tbox <- t_running_apx_quantiles(x,ptiles,max_order=3,time=times,wts=wts,window=t_window,na_rm=na_rm,normalize_wts=nw),NA)
					expect_equal(box,tbox,tolerance=1e-8)
				}
			}# UNFOLD
		}
	}
})#UNFOLD

context("t_running vs slow version")
test_that("check em",{#FOLDUP
	skip_on_cran()

	set.char.seed("91b0bd37-0b8e-49d6-8333-039a7d7f7dd5")
	na_rm <- FALSE
	for (xlen in c(40,90)) {# FOLDUP
		x <- rnorm(xlen)
		for (times in list(NULL,cumsum(runif(length(x),min=0.2,max=0.4)))) {
			for (wts in list(NULL,rep(1L,xlen),runif(xlen,min=1.1,max=2.1))) { 
				wts_as_delta <- is.null(times) & !is.null(wts)
				if (!is.null(times) || (wts_as_delta && !is.null(wts))) {
					for (window in c(11.5,20.5,Inf)) { # FOLDUP
						for (lb_time in list(NULL,3+cumsum(runif(10,min=0.4,max=1.1)))) {
							slow <- slow_t_running_sum(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta)
							expect_error(fast <- t_running_sum(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta),NA)
							expect_equal(fast,slow,tolerance=1e-8)
							
							slow <- slow_t_running_mean(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta)
							expect_error(fast <- t_running_mean(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta),NA)
							expect_equal(fast,slow,tolerance=1e-8)

							for (nw in c(TRUE,FALSE)) { 
								slow <- slow_t_running_sd(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw)
								expect_error(fast <- t_running_sd(x,time=times,wts=wts,window=window,lb_time=lb_time,min_df=1,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw),NA)
								expect_equal(fast,slow,tolerance=1e-8)

								slow <- slow_t_running_skew(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw)
								expect_error(fast <- t_running_skew(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw),NA)
								expect_equal(fast,slow,tolerance=1e-8)

								slow <- slow_t_running_kurt(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw)
								expect_error(fast <- t_running_kurt(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw),NA)
								expect_equal(fast,slow,tolerance=1e-8)

								slow <- slow_t_running_sd3(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw)
								expect_error(fast <- t_running_sd3(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,wts_as_delta=wts_as_delta,normalize_wts=nw),NA)
								# ignore the df computation in slow when empty
								slow[fast[,3]==0,3] <- 0
								slow[is.na(fast[,1]),1] <- NA
								expect_equal(fast,slow,tolerance=1e-8)

								slow <- slow_t_running_skew4(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,normalize_wts=nw)
								expect_error(fast <- t_running_skew4(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,normalize_wts=nw),NA)
								# ignore the df computation in slow when empty
								okrow <- !is.na(fast[,4]) & fast[,4] > 3 & row(fast)[,4] > 3
								expect_equal(fast[okrow,],slow[okrow,],tolerance=1e-8)

								slow <- slow_t_running_kurt5(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,normalize_wts=nw)
								expect_error(fast <- t_running_kurt5(x,time=times,wts=wts,window=window,lb_time=lb_time,na_rm=na_rm,normalize_wts=nw),NA)
								okrow <- !is.na(fast[,5]) & fast[,5] > 4 & row(fast)[,5] > 4
								expect_equal(fast[okrow,],slow[okrow,],tolerance=1e-8)
							}
						}
					}# UNFOLD
				}
			}
		}
	}# UNFOLD
})#UNFOLD

context("t_running_sd")
# t_running_sd is a bellwether for the other methods
# as it goes, so goes the other Welford based functions
test_that("check it",{#FOLDUP
	skip_on_cran()

	set.char.seed("79f60eda-7799-46e6-9096-6817b2d4473b")

	na_rm <- FALSE
	for (xlen in c(20,50)) {# FOLDUP
		x <- rnorm(xlen)
		for (times in list(NULL,cumsum(runif(length(x),min=0.2,max=0.4)))) {
			for (wts in list(NULL,rep(1L,xlen),runif(xlen,min=1.2,max=2.1))) { 
				wts_as_delta <- is.null(times) & !is.null(wts)
				if (!is.null(times) || (wts_as_delta && !is.null(wts))) {
					for (window in c(11.5,20.5,Inf)) { # FOLDUP
						for (lb_time in list(NULL,cumsum(runif(20,min=0.2,max=1)))) {
							for (nw in c(TRUE,FALSE)) { 
								expect_error(slow <- reference_t_running_sd(x,time=times,wts=wts,wts_as_delta=TRUE,window=window,lb_time=lb_time,na_rm=na_rm,min_df=1,normalize_wts=nw),NA)
								expect_error(fast <- t_running_sd(x,time=times,wts=wts,wts_as_delta=TRUE,used_df=1,window=window,lb_time=lb_time,min_df=1,na_rm=na_rm,normalize_wts=nw),NA)
								expect_equal(fast,slow,tolerance=1e-7)
							}
						}
					}# UNFOLD
				}
			}
		}
	}# UNFOLD
})#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
