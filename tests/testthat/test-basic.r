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
# Created: 2016.03.25
# Copyright: Steven E. Pav, 2016-2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

# code runs at all
context("code runs: sd, skew, kurt")
test_that("sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("569dd47d-f9e5-40e4-b2ac-e5dbb4771a53")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')
	wts <- 100 * runif(length(x))

	for (na_rm in c(FALSE,TRUE)) {
		for (fnc in list(function(x) { x },as.integer,function(x) { as.logical(x > 0) })) {
			cls_x <- fnc(x)
			expect_error(sd3(cls_x,na_rm=na_rm),NA)
			expect_error(skew4(cls_x,na_rm=na_rm),NA)
			expect_error(kurt5(cls_x,na_rm=na_rm),NA)
			expect_error(cent_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(std_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(cent_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(std_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)

			# weights!?
			for (wfnc in list(function(x) { x },as.integer,function(x) { as.logical(x > 10) })) {
				cls_wts <- wfnc(wts)
				expect_error(sd3(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(skew4(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(kurt5(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(kurt5(cls_x,na_rm=na_rm,wts=cls_wts,check_wts=TRUE),NA)
				
				expect_error(cent_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(std_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(cent_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(std_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
			}
		}
	}
	# check some things about weights;
	expect_error(sd3(x,wts=q))
	expect_error(sd3(x,wts=rep(-1,length(x)),check_wts=TRUE))

	for (checkme in list(sd3,skew4,kurt5,cent_moments,std_moments,cent_cumulants,std_cumulants)) {
		expect_error(checkme(q))
		expect_error(checkme(x,wts=q))
		expect_error(checkme(x,wts=rep(-1,length(x)),check_wts=TRUE))
	}
})#UNFOLD
context("code runs: cosum and comoment")
test_that("cosum and comoment run without error",{#FOLDUP
	set.char.seed("ff34b509-a113-41c9-8517-aa72792c42f7")
	x <- matrix(rnorm(30*4),ncol=4)
	y <- matrix(as.integer(x),ncol=ncol(x))
	z <- matrix(as.logical(y),ncol=ncol(x))
	q <- matrix(letters[1:24],ncol=4)

	for (na_omit in c(FALSE,TRUE)) {
		expect_error(cent_cosums(x,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(x,max_order=2L,used_df=1L,na_omit=na_omit),NA)

		expect_error(cent_cosums(y,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(y,max_order=2L,used_df=1L,na_omit=na_omit),NA)

		expect_error(cent_cosums(z,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(z,max_order=2L,used_df=1L,na_omit=na_omit),NA)
	}

	expect_error(cent_cosums(x,max_order=4L))
	expect_error(cent_cosums(q,max_order=2L))
	expect_error(cent_comoments(q,max_order=2L))
})#UNFOLD
context("code runs: running_sd etc")
test_that("running sd, skew, kurt run without error",{#FOLDUP
	skip_on_cran()

	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (window in c(50,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			expect_error(running_sum(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_mean(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sd(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sd3(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew4(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt5(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
			expect_error(running_std_moments(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cumulants(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)

			expect_error(running_sum(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_mean(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sd(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sd3(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew4(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt5(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
			expect_error(running_std_moments(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cumulants(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_apx_quantiles(y,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)

			expect_error(running_sd(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sd3(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_skew4(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_kurt5(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cent_moments(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
			expect_error(running_std_moments(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_cumulants(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_apx_quantiles(z,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
		}
	}
	for (min_df in c(2L,10L)) {
		expect_error(running_mean(x,window=window,min_df=min_df),NA)
		expect_error(running_sd(x,window=window,min_df=min_df),NA)
		expect_error(running_skew(x,window=window,min_df=min_df),NA)
		expect_error(running_kurt(x,window=window,min_df=min_df),NA)
		expect_error(running_sd3(x,window=window,min_df=min_df),NA)
		expect_error(running_skew4(x,window=window,min_df=min_df),NA)
		expect_error(running_kurt5(x,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,min_df=min_df,max_order_only=TRUE),NA)
		expect_error(running_std_moments(x,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cumulants(x,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,min_df=min_df),NA)

		expect_error(running_mean(y,window=window,min_df=min_df),NA)
		expect_error(running_sd(y,window=window,min_df=min_df),NA)
		expect_error(running_skew(y,window=window,min_df=min_df),NA)
		expect_error(running_kurt(y,window=window,min_df=min_df),NA)
		expect_error(running_sd3(y,window=window,min_df=min_df),NA)
		expect_error(running_skew4(y,window=window,min_df=min_df),NA)
		expect_error(running_kurt5(y,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(y,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(y,max_order=5L,window=window,min_df=min_df,max_order_only=TRUE),NA)
		expect_error(running_std_moments(y,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cumulants(y,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_apx_quantiles(y,p=ptiles,max_order=5L,window=window,min_df=min_df),NA)

		expect_error(running_mean(z,window=window,min_df=min_df),NA)
		expect_error(running_sd(z,window=window,min_df=min_df),NA)
		expect_error(running_skew(z,window=window,min_df=min_df),NA)
		expect_error(running_kurt(z,window=window,min_df=min_df),NA)
		expect_error(running_sd3(z,window=window,min_df=min_df),NA)
		expect_error(running_skew4(z,window=window,min_df=min_df),NA)
		expect_error(running_kurt5(z,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(z,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cent_moments(z,max_order=5L,window=window,min_df=min_df,max_order_only=TRUE),NA)
		expect_error(running_std_moments(z,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_cumulants(z,max_order=5L,window=window,min_df=min_df),NA)
		expect_error(running_apx_quantiles(z,p=ptiles,max_order=5L,window=window,min_df=min_df),NA)
	}

	expect_error(running_sum(q))
	expect_error(running_mean(q))
	expect_error(running_sd(q))
	expect_error(running_skew(q))
	expect_error(running_kurt(q))
	expect_error(running_sd3(q))
	expect_error(running_skew4(q))
	expect_error(running_kurt5(q))
	expect_error(running_cent_moments(q,max_order=5L))
	expect_error(running_cent_moments(q,max_order=5L,max_order_only=TRUE))
	expect_error(running_std_moments(q,max_order=5L))
	expect_error(running_cumulants(q,max_order=5L))
	expect_error(running_apx_quantiles(q,p=ptiles,max_order=5L))
	expect_error(running_apx_quantiles(x,p=q,max_order=5L))
	expect_error(running_apx_median(q,p=ptiles,max_order=5L))
})#UNFOLD
context("code runs: running_foo and weights")
test_that("running foo and weights",{#FOLDUP
	skip_on_cran()

	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	nel <- 20
	xall <- list(rnorm(nel),
							 as.integer(rnorm(nel,sd=100)))
	wall <- list(rep(1.0,nel),
							 rep(2L,nel),
							 NULL)

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (x in xall) {
		for (wts in wall) {
			for (window in c(5,21,Inf,NULL)) {
				for (na_rm in c(FALSE,TRUE)) {
						#running_sum(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_mean(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						expect_error(running_sd(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_skew(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_kurt(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_sd3(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_skew4(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_kurt5(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_cent_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_cent_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
						expect_error(running_std_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_cumulants(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_apx_quantiles(x,wts=wts,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_apx_median(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				}
			}
		}
	}
})#UNFOLD
context("code runs: hits heywood")
test_that("hit heywood branch",{#FOLDUP
	skip_on_cran()

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	set.char.seed("3d318f1d-9921-4a20-84fc-c5ffc722d52c")
	xvals <- list(rnorm(1e5,mean=1e10))
	x <- rnorm(1e4)
	x[x < 1.0] <- NA
	xvals[[length(xvals)+1]] <- x
	x <- rnorm(1e4)
	x[x < 1.5] <- NA
	xvals[[length(xvals)+1]] <- x

	for (x in xvals) {
		window <- 500L
		restart_period <- 100000L
		na_rm <- TRUE
		# no heywood branch for these?
		expect_error(running_sum(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_mean(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)

		expect_error(running_sd3(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_skew4(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_kurt5(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,max_order_only=TRUE,na_rm=na_rm),NA)
		expect_error(running_std_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cumulants(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_median(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)

	}
})#UNFOLD
context("code runs: NA restart")
test_that("NA restart period?",{#FOLDUP
	skip_on_cran()

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	set.char.seed("3d318f1d-9921-4a20-84fc-c5ffc722d52c")
	xvals <- list(rnorm(1e2,mean=1e10))
	x <- rnorm(1e3)
	x[x < 1.0] <- NA
	xvals[[length(xvals)+1]] <- x
	x <- rnorm(1e3)
	x[x < 1.5] <- NA
	xvals[[length(xvals)+1]] <- x

	for (x in xvals) {
		window <- 50L
		restart_period <- NA_integer_
		na_rm <- TRUE
		expect_error(running_sum(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_mean(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)

		expect_error(running_sd3(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_skew4(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_kurt5(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,max_order_only=TRUE,na_rm=na_rm),NA)
		expect_error(running_std_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cumulants(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_median(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)

	}
})#UNFOLD
context("code runs: running adjustments")
test_that("running adjustments",{#FOLDUP
	skip_on_cran()

	set.char.seed("e36291f6-70e0-4412-9c50-bc46b6ab8639")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (window in c(50,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			expect_error(running_centered(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_scaled(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_zscored(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
			expect_error(running_tstat(x,window=window,restart_period=50L,na_rm=na_rm),NA)

			expect_error(running_centered(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_scaled(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_zscored(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
			expect_error(running_tstat(y,window=window,restart_period=50L,na_rm=na_rm),NA)

			expect_error(running_centered(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_scaled(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_zscored(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
			expect_error(running_tstat(z,window=window,restart_period=50L,na_rm=na_rm),NA)
		}
	}
	window <- 10L

	for (na_rm in c(FALSE,TRUE)) {
		expect_error(running_centered(x,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_scaled(x,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_zscored(x,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
		expect_error(running_tstat(x,window=window,restart_period=50L,na_rm=na_rm),NA)

		expect_error(running_centered(y,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_scaled(y,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_zscored(y,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
		expect_error(running_tstat(y,window=window,restart_period=50L,na_rm=na_rm),NA)

		expect_error(running_centered(z,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_scaled(z,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_zscored(z,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm),NA)
		expect_error(running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
		expect_error(running_tstat(z,window=window,restart_period=50L,na_rm=na_rm),NA)
	}

	for (min_df in c(2L,10L)) {
		expect_error(running_centered(x,window=window,min_df=min_df),NA)
		expect_error(running_scaled(x,window=window,min_df=min_df),NA)
		expect_error(running_zscored(x,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(x,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(x,window=window,min_df=min_df,compute_se=TRUE),NA)
		expect_error(running_tstat(x,window=window,min_df=min_df),NA)

		expect_error(running_centered(y,window=window,min_df=min_df),NA)
		expect_error(running_scaled(y,window=window,min_df=min_df),NA)
		expect_error(running_zscored(y,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(y,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(y,window=window,min_df=min_df,compute_se=TRUE),NA)
		expect_error(running_tstat(y,window=window,min_df=min_df),NA)

		expect_error(running_centered(z,window=window,min_df=min_df),NA)
		expect_error(running_scaled(z,window=window,min_df=min_df),NA)
		expect_error(running_zscored(z,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(z,window=window,min_df=min_df),NA)
		expect_error(running_sharpe(z,window=window,min_df=min_df,compute_se=TRUE),NA)
		expect_error(running_tstat(z,window=window,min_df=min_df),NA)
	}

	expect_error(running_centered(q))
	expect_error(running_scaled(q))
	expect_error(running_zscored(q))
	expect_error(running_sharpe(q))
	expect_error(running_sharpe(q,compute_se=TRUE))
	expect_error(running_tstat(q))

	expect_error(running_tstat(x,window='FOO'))
	expect_error(running_tstat(x,window=-20L))
	expect_error(running_tstat(x,window=20L,restart_period='FOO'))
})#UNFOLD
context("code runs: t_running_sd etc")
test_that("t_running sd, skew, kurt run without error",{#FOLDUP
	skip_on_cran()

	set.char.seed("1cfb2b84-7eae-42f6-90ad-2552e74a5ad0")
	x <- rnorm(100)
	times <- seq_along(x)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (thingy in list(x,y,z)) { 
		for (window in c(50,Inf)) {
			for (na_rm in c(FALSE,TRUE)) {
				expect_error(t_running_sum(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_mean(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_sd(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_skew(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_kurt(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_sd3(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_skew4(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_kurt5(thingy,time=times,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_cent_moments(thingy,time=times,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_cent_moments(thingy,time=times,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
				expect_error(t_running_std_moments(thingy,time=times,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_cumulants(thingy,time=times,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(t_running_apx_quantiles(thingy,time=times,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			}
		}
	}

	for (thingy in list(x,y,z)) { 
		for (min_df in c(2L,10L)) {
			expect_error(t_running_sum(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_mean(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_sd(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_skew(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_kurt(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_sd3(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_skew4(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_kurt5(thingy,time=times,window=window,min_df=min_df),NA)
			expect_error(t_running_cent_moments(thingy,time=times,max_order=5L,window=window,min_df=min_df),NA)
			expect_error(t_running_cent_moments(thingy,time=times,max_order=5L,window=window,min_df=min_df,max_order_only=TRUE),NA)
			expect_error(t_running_std_moments(thingy,time=times,max_order=5L,window=window,min_df=min_df),NA)
			expect_error(t_running_cumulants(thingy,time=times,max_order=5L,window=window,min_df=min_df),NA)
			expect_error(t_running_apx_quantiles(thingy,time=times,p=ptiles,max_order=5L,window=window,min_df=min_df),NA)
		}
	}


	expect_error(t_running_sum(q))
	expect_error(t_running_mean(q))
	expect_error(t_running_sd(q))
	expect_error(t_running_skew(q))
	expect_error(t_running_kurt(q))
	expect_error(t_running_sd3(q))
	expect_error(t_running_skew4(q))
	expect_error(t_running_kurt5(q))
	expect_error(t_running_cent_moments(q,max_order=5L))
	expect_error(t_running_cent_moments(q,max_order=5L,max_order_only=TRUE))
	expect_error(t_running_std_moments(q,max_order=5L))
	expect_error(t_running_cumulants(q,max_order=5L))
	expect_error(t_running_apx_quantiles(q,p=ptiles,max_order=5L))
	expect_error(t_running_apx_quantiles(x,p=q,max_order=5L))
	expect_error(t_running_apx_median(q,p=ptiles,max_order=5L))
})#UNFOLD
#test_that("running foo and weights",{#FOLDUP
	## hey, Volkswagon called while you were out:
	#skip_on_cran()

	#set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	#nel <- 20
	#xall <- list(rnorm(nel),
							 #as.integer(rnorm(nel,sd=100)))
	#wall <- list(rep(1.0,nel),
							 #rep(2L,nel),
							 #NULL)

	#ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	#for (x in xall) {
		#for (wts in wall) {
			#for (window in c(5,21,Inf,NULL)) {
				#for (na_rm in c(FALSE,TRUE)) {
						##running_sum(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						##running_mean(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_sd(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_skew(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_kurt(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_sd3(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_skew4(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_kurt5(x,wts=wts,window=window,restart_period=50L,na_rm=na_rm)
						#running_cent_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
						#running_cent_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE)
						#running_std_moments(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
						#running_cumulants(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
						#running_apx_quantiles(x,wts=wts,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
						#running_apx_median(x,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
				#}
			#}
		#}
	#}

	## sentinel
	#expect_true(TRUE)
#})#UNFOLD




# 2FIX: check the effects of NA

context("weighted estimation")# FOLDUP
test_that("running adjustments",{#FOLDUP
	skip_on_cran()

	set.char.seed("e36291f6-70e0-4412-9c50-bc46b6ab8639")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (window in c(50,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			running_centered(x,window=window,restart_period=50L,na_rm=na_rm)
			running_scaled(x,window=window,restart_period=50L,na_rm=na_rm)
			running_zscored(x,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
			running_tstat(x,window=window,restart_period=50L,na_rm=na_rm)

			running_centered(y,window=window,restart_period=50L,na_rm=na_rm)
			running_scaled(y,window=window,restart_period=50L,na_rm=na_rm)
			running_zscored(y,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
			running_tstat(y,window=window,restart_period=50L,na_rm=na_rm)

			running_centered(z,window=window,restart_period=50L,na_rm=na_rm)
			running_scaled(z,window=window,restart_period=50L,na_rm=na_rm)
			running_zscored(z,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
			running_tstat(z,window=window,restart_period=50L,na_rm=na_rm)
		}
	}
	window <- 10L

	for (na_rm in c(FALSE,TRUE)) {
		running_centered(x,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(x,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(x,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
		running_tstat(x,window=window,restart_period=50L,na_rm=na_rm)

		running_centered(y,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(y,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(y,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
		running_tstat(y,window=window,restart_period=50L,na_rm=na_rm)

		running_centered(z,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(z,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(z,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE)
		running_tstat(z,window=window,restart_period=50L,na_rm=na_rm)
	}

	for (min_df in c(2L,10L)) {
		running_centered(x,window=window,min_df=min_df)
		running_scaled(x,window=window,min_df=min_df)
		running_zscored(x,window=window,min_df=min_df)
		running_sharpe(x,window=window,min_df=min_df)
		running_sharpe(x,window=window,min_df=min_df,compute_se=TRUE)
		running_tstat(x,window=window,min_df=min_df)

		running_centered(y,window=window,min_df=min_df)
		running_scaled(y,window=window,min_df=min_df)
		running_zscored(y,window=window,min_df=min_df)
		running_sharpe(y,window=window,min_df=min_df)
		running_sharpe(y,window=window,min_df=min_df,compute_se=TRUE)
		running_tstat(y,window=window,min_df=min_df)

		running_centered(z,window=window,min_df=min_df)
		running_scaled(z,window=window,min_df=min_df)
		running_zscored(z,window=window,min_df=min_df)
		running_sharpe(z,window=window,min_df=min_df)
		running_sharpe(z,window=window,min_df=min_df,compute_se=TRUE)
		running_tstat(z,window=window,min_df=min_df)
	}

	expect_error(running_centered(q))
	expect_error(running_scaled(q))
	expect_error(running_zscored(q))
	expect_error(running_sharpe(q))
	expect_error(running_sharpe(q,compute_se=TRUE))
	expect_error(running_tstat(q))

	expect_error(running_tstat(x,window='FOO'))
	expect_error(running_tstat(x,window=-20L))
	expect_error(running_tstat(x,window=20L,restart_period='FOO'))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# UNFOLD

context("order of parameters")# FOLDUP
# named parameters can be given in any order?

# UNFOLD


#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
