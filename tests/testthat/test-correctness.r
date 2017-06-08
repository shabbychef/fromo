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
#UNFOLD

context("first moments")#FOLDUP
test_that("sd, skew, kurt are correct",{#FOLDUP
	set.char.seed("c4007dba-2010-481e-abe5-f07d3ce94eb4")
	x <- rnorm(1000)

	sid <- sd3(x)
	ske <- skew4(x)
	krt <- kurt5(x)

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

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("weight scaling what you expect",{#FOLDUP
	set.char.seed("efaa75ac-bb9e-4e4a-a375-7028f099366e")
	x <- rnorm(50)
	wts <- runif(length(x))

	sid_1 <- sd3(x,wts=wts,normalize_wts=FALSE,sg_df=0)
	ske_1 <- skew4(x,wts=wts,normalize_wts=FALSE,sg_df=0)
	krt_1 <- kurt5(x,wts=wts,normalize_wts=FALSE,sg_df=0)

	sid_2 <- sd3(x,wts=2*wts,normalize_wts=FALSE,sg_df=0)
	ske_2 <- skew4(x,wts=2*wts,normalize_wts=FALSE,sg_df=0)
	krt_2 <- kurt5(x,wts=2*wts,normalize_wts=FALSE,sg_df=0)

	expect_equal(sid_1 * c(1,1,2),sid_2,tolerance=1e-9)
	expect_equal(ske_1 * c(1,1,1,2),ske_2,tolerance=1e-9)
	expect_equal(krt_1 * c(1,1,1,1,2),krt_2,tolerance=1e-9)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("weighted sd, skew, kurt are correct",{#FOLDUP
	set.char.seed("4e17d837-69c1-41d1-906f-c82224d7ce41")
	x <- rnorm(1000)
	wts <- runif(length(x))

	sid <- sd3(x,wts=wts,normalize_wts=TRUE)
	ske <- skew4(x,wts=wts,normalize_wts=TRUE)
	krt <- kurt5(x,wts=wts,normalize_wts=TRUE)
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

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD
#
#
tomat <- function(cbound) {
	dumbv <- as.matrix(cbound)
	attr(dumbv,'dimnames') <- NULL 
	dumbv
}

context("running ops")# FOLDUP
test_that("running ops are correct",{#FOLDUP
	# hey, Volkswagon called while you were out:
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
						fastv <- running_sd(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(dumb_sd)
						expect_equal(dumbv[2:xlen],fastv[2:xlen],tolerance=1e-7 * toler)

						fastv <- running_sd3(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(cbind(dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[2:xlen,],fastv[2:xlen,],tolerance=1e-7 * toler)

						# skew
						fastv <- running_skew(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(dumb_skew)
						expect_equal(dumbv[3:xlen],fastv[3:xlen],tolerance=1e-6 * toler)

						fastv <- running_skew4(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(cbind(dumb_skew,dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[3:xlen,],fastv[3:xlen,],tolerance=1e-7 * toler)

						# excess kurtosis
						fastv <- running_kurt(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(dumb_exkurt)
						expect_equal(dumbv[4:xlen],fastv[4:xlen],tolerance=1e-6 * toler)

						fastv <- running_kurt5(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(cbind(dumb_exkurt,dumb_skew,dumb_sd,dumb_mean,dumb_count))
						expect_equal(dumbv[4:xlen,],fastv[4:xlen,],tolerance=1e-6 * toler)

						# higher order moments

						fastv <- running_cent_moments(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(cbind(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count))
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-6 * toler)

						fastv <- running_cent_moments(x,window=window,max_order=6L,max_order_only=TRUE,used_df=0L,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(dumb_cmom6)
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-7 * toler)

						fastv <- running_std_moments(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm)
						dumbv <- tomat(cbind(dumb_cmom6 / (dumb_cmom2^3),dumb_cmom5 / (dumb_cmom2^2.5),dumb_cmom4 / (dumb_cmom2^2.0),dumb_cmom3 / (dumb_cmom2^1.5),sqrt(dumb_cmom2),dumb_mean,dumb_count))
						expect_equal(dumbv[6:xlen,],fastv[6:xlen,],tolerance=1e-7 * toler)

						# running sum and mean
						# do sums twice
						fastv <- running_sum(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- dumb_sum
						expect_equal(dumbv[2:xlen],fastv[2:xlen,],tolerance=1e-8 * toler)

						fastv <- running_mean(x,window=window,restart_period=restart_period,na_rm=na_rm)
						dumbv <- dumb_mean
						expect_equal(dumbv[2:xlen],fastv[2:xlen,],tolerance=1e-8 * toler)

						if (require(PDQutils)) {
							# cumulants
							fastv <- running_cumulants(x,window=window,max_order=6L,used_df=0L,restart_period=restart_period,na_rm=na_rm)
							pre_dumbv <- cbind(dumb_cmom6,dumb_cmom5,dumb_cmom4,dumb_cmom3,dumb_cmom2,dumb_mean,dumb_count)
							dumbv <- t(sapply(seq_along(x),function(iii) { 
															rv <- rev(PDQutils::moment2cumulant(c(0,rev(pre_dumbv[iii,1:(ncol(pre_dumbv)-2)]))))
															rv <- rv[-length(rv)]
															c(rv,pre_dumbv[iii,ncol(pre_dumbv) + (-1:0)])
								},simplify='matrix'))
							expect_equal(max(abs(dumbv[6:xlen,] - fastv[6:xlen,])),0,tolerance=1e-8 * toler)

							# quantiles
							fastv <- running_apx_quantiles(x,ptiles,max_order=ncol(dumbv)-1,used_df=0L,window=window,restart_period=restart_period,na_rm=na_rm)
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

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running adjustments are correct",{#FOLDUP
	# hey, Volkswagon called while you were out:
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

					fastv <- running_centered(x,window=window,restart_period=restart_period,na_rm=na_rm)
					# the dumb value:
					dumbv <- x - dumb_mean;
					expect_equal(max(abs(dumbv - fastv)),0,tolerance=1e-12)

					fastv <- running_scaled(x,window=window,restart_period=restart_period,na_rm=na_rm)
					# the dumb value:
					dumbv <- x / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_zscored(x,window=window,restart_period=restart_period,na_rm=na_rm)
					# the dumb value:
					dumbv <- (x - dumb_mean) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm)
					# the dumb value:
					dumbv <- dumb_mean / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_tstat(x,window=window,restart_period=restart_period,na_rm=na_rm)
					# the dumb value:
					dumbv <- (dumb_mean * sqrt(dumb_count)) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm,compute_se=TRUE)
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

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running weights work correctly",{#FOLDUP
	# hey, Volkswagon called while you were out:
	skip_on_cran()

	set.char.seed("b82d252c-681b-4b98-9bb3-ffd17feeb4a1")
	na_rm <- FALSE

	for (xlen in c(20,50)) {
		x <- rnorm(xlen)
		wtlist <- list(rep(1L,xlen), 2+5*runif(xlen))
		for (wts in wtlist) {
			for (window in c(5,30,Inf)) { # FOLDUP

				# 2FIX: add to this!
				dumb_count <- sapply(seq_along(x),function(iii) { sum(sign(abs(x[max(1,iii-window+1):iii])+1),na.rm=na_rm) },simplify=TRUE)
				dumb_sumwt <- sapply(seq_along(x),function(iii) { sum(wts[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
				dumb_mean <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * x[mydx],na.rm=na_rm) / dumb_sumwt[iii]
							 },simplify=TRUE)
				dumb_var <- sapply(seq_along(x),function(iii) { 
															mydx <- max(1,iii-window+1):iii
															mywts <- wts[mydx]
															sum(mywts * (x[mydx] - dumb_mean[iii])^2,na.rm=na_rm) / (dumb_sumwt[iii] - 1)
							 },simplify=TRUE)
				dumb_sd <- sqrt(dumb_var)

				# ack! need weights!
				fastv <- running_mean(x,window=window,na_rm=na_rm)
				fastv <- running_sd(x,wts=wts,window=window,na_rm=na_rm,normalize_wts=TRUE)


				#fastv <- running_centered(x,window=window,restart_period=restart_period,na_rm=na_rm)
				## the dumb value:
				#dumbv <- x - dumb_mean;
				#expect_equal(max(abs(dumbv - fastv)),0,tolerance=1e-12)

				#fastv <- running_scaled(x,window=window,restart_period=restart_period,na_rm=na_rm)
				## the dumb value:
				#dumbv <- x / dumb_sd
				#expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

				#fastv <- running_zscored(x,window=window,restart_period=restart_period,na_rm=na_rm)
				## the dumb value:
				#dumbv <- (x - dumb_mean) / dumb_sd
				#expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

				#fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm)
				## the dumb value:
				#dumbv <- dumb_mean / dumb_sd
				#expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

				#fastv <- running_tstat(x,window=window,restart_period=restart_period,na_rm=na_rm)
				## the dumb value:
				#dumbv <- (dumb_mean * sqrt(dumb_count)) / dumb_sd
				#expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

				#fastv <- running_sharpe(x,window=window,restart_period=restart_period,na_rm=na_rm,compute_se=TRUE)
				## the dumb value:
				#dumb_sr <- dumb_mean / dumb_sd
				#expect_equal(max(abs(dumb_sr[2:length(x)] - fastv[2:length(x),1])),0,tolerance=1e-12)

				#if (require(moments)) {
					#dumb_skew <- sapply(seq_along(x),function(iii) { moments::skewness(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
					#dumb_exkurt <- sapply(seq_along(x),function(iii) { moments::kurtosis(x[max(1,iii-window+1):iii],na.rm=na_rm) - 3.0 },simplify=TRUE)
					#dumb_merse <- sqrt((1 + 0.25 * (2+dumb_exkurt) * dumb_sr^2 - dumb_skew * dumb_sr) / dumb_count)
					#expect_equal(max(abs(dumb_merse[5:length(x)] - fastv[5:length(x),2])),0,tolerance=1e-9)
				#}
				
			}# UNFOLD
		}
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# UNFOLD
context("monoid nonsense")# FOLDUP
test_that("join/unjoin",{#FOLDUP
	set.char.seed("1325a51e-1584-4f89-9ea3-f15223a223d9")

	x1 <- rnorm(1e3,mean=1)
	x2 <- rnorm(1e3,mean=1)
	max_ord <- 6L
	rs1 <- cent_sums(x1,max_ord)
	rs2 <- cent_sums(x2,max_ord)
	rs3 <- cent_sums(c(x1,x2),max_ord)
	rs3alt <- join_cent_sums(rs1,rs2)
	expect_lt(max(abs(rs3 - rs3alt)),1e-7)

	rs1alt <- unjoin_cent_sums(rs3,rs2)
	rs2alt <- unjoin_cent_sums(rs3,rs1)
	expect_lt(max(abs(rs1 - rs1alt)),1e-7)
	expect_lt(max(abs(rs2 - rs2alt)),1e-7)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("cosums are sane",{#FOLDUP
	set.char.seed("0020a8c0-ff6a-447c-a9bf-c6cc7160195f")

	x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	max_ord <- 2L
	rs1 <- cent_comoments(x1,max_ord,used_df=1L)
	expect_equal(rs1[1,1],nrow(x1))
	expect_equal(rs1[1,1 + (1:ncol(x1))],colMeans(x1),tolerance=1e-7)
	expect_equal(rs1[1 + (1:ncol(x1)),1],colMeans(x1),tolerance=1e-7)
	expect_equal(rs1[1 + (1:ncol(x1)),1 + (1:ncol(x1))],cov(x1),tolerance=1e-7)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("join/unjoin cosums",{#FOLDUP
	set.char.seed("9ecdda29-aaae-4f88-9fe7-4418846ca54c")

	x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	x2 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	max_ord <- 2L
	rs1 <- cent_cosums(x1,max_ord)
	rs2 <- cent_cosums(x2,max_ord)
	rs3 <- cent_cosums(rbind(x1,x2),max_ord)
	rs3alt <- join_cent_cosums(rs1,rs2)
	expect_lt(max(abs(rs3 - rs3alt)),1e-7)

	rs1alt <- unjoin_cent_cosums(rs3,rs2)
	rs2alt <- unjoin_cent_cosums(rs3,rs1)
	expect_lt(max(abs(rs1 - rs1alt)),1e-7)
	expect_lt(max(abs(rs2 - rs2alt)),1e-7)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# UNFOLD
# 2FIX: check the effects of NA

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
