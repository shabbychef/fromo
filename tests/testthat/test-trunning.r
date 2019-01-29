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
# Created: 2019.01.15
# Copyright: Steven E. Pav, 2016-2019
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
#UNFOLD

# code runs at all
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
	#expect_error(t_running_apx_quantiles(x,p=q,max_order=5L))
	expect_error(t_running_apx_median(q,p=ptiles,max_order=5L))
})#UNFOLD

context("code runs: t_running_foo and weights")
test_that("t_running foo and weights",{#FOLDUP
	skip_on_cran()

	set.char.seed("95aa6ed1-14dc-473a-a7aa-8dc9869b1dbc")
	nel <- 30

	xna <- rnorm(2*nel)
	xna[xna < 0] <- NA

	xall <- list(rnorm(nel),
							 xna,
							 as.integer(rnorm(nel,sd=100)))

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	rp <- 5L

	for (thingy in xall) {
		for (times in list(NULL,cumsum(runif(length(thingy),min=0.2,max=0.4)))) {
			wna <- runif(length(thingy),min=1.0,max=3.0)
			wna[wna < 1.4] <- NA
			wall <- list(rep(1.0,length(thingy)),
									 runif(length(thingy),min=0.9,max=3.5),
									 wna,
									 NULL)
			for (wts in wall) {
				wts_as_delta <- is.null(times) & !is.null(wts)
				can_test <- (!wts_as_delta || is.null(wts) || !any(is.na(wts))) && (!is.null(times) || wts_as_delta)
				if (can_test) { 
					for (window in c(5.5,21.23,Inf,NULL)) {
						for (na_rm in c(FALSE,TRUE)) {
							for (lb_time in list(NULL,3+cumsum(runif(10,min=0.4,max=1.1)))) {
								expect_error(t_running_sum(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_mean(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_sd(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_skew(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_kurt(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_sd3(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_skew4(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_kurt5(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_cent_moments(thingy,time=times,wts=wts,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_cent_moments(thingy,time=times,wts=wts,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,max_order_only=TRUE),NA)
								expect_error(t_running_std_moments(thingy,time=times,wts=wts,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_cumulants(thingy,time=times,wts=wts,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_apx_quantiles(thingy,time=times,wts=wts,p=ptiles,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_apx_median(thingy,time=times,wts=wts,max_order=5L,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_sharpe(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
								expect_error(t_running_sharpe(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,compute_se=TRUE),NA)
								expect_error(t_running_tstat(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)

								# NOTA BENE: these three do not accept an lb_time b/c they are associated with the thingy;
								# so we *do* expect these to error.
								expect_error(t_running_centered(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
								expect_error(t_running_scaled(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
								expect_error(t_running_zscored(thingy,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
							}

							for (lookahead in c(0,8.3)) {
								expect_error(t_running_centered(thingy,time=times,wts=wts,window=window,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
								expect_error(t_running_scaled(thingy,time=times,wts=wts,window=window,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
								expect_error(t_running_zscored(thingy,time=times,wts=wts,window=window,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
							}
						}
					}
				}
			}
		}
	}

})#UNFOLD
test_that("t_running foo logical weights",{# FOLDUP
	skip_on_cran()

	set.char.seed("422008f4-5782-4f01-a357-5b254a16660e")
	nel <- 30

	xna <- rnorm(2*nel)
	xna[xna < 0] <- NA

	rp <- 5L
	times <- cumsum(runif(length(xna),min=0.5,max=1.5))
	window <- 17.5
	na_rm <- TRUE
	lb_time <- NULL
	wts <- rnorm(length(xna)) > 0
	expect_error(t_running_sum(xna,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=FALSE,restart_period=rp,na_rm=na_rm),NA)
	expect_error(t_running_mean(xna,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=FALSE,restart_period=rp,na_rm=na_rm),NA)
	expect_error(t_running_sd(xna,time=times,wts=wts,window=window,lb_time=lb_time,wts_as_delta=FALSE,restart_period=rp,na_rm=na_rm),NA)

})# UNFOLD
context("code runs: t_running_foo NA weights")
test_that("t_running foo  NA weights",{#FOLDUP
	skip_on_cran()

	set.char.seed("82da516f-6185-4f65-8576-f8c02e00d61e")
	nel <- 25

	xall <- list(rnorm(nel))

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	rp <- 5L
	wts_as_delta <- FALSE

	for (thingy in xall) {
		wts <- runif(length(thingy),min=0.9,max=1.2)
		wts[wts < 1] <- NA
		for (times in list(cumsum(runif(length(thingy),min=0.2,max=0.4)))) {
			for (na_rm in c(FALSE,TRUE)) {
				for (lb_time in list(max(times) + c(1:2),3+cumsum(runif(10,min=0.4,max=1.1)))) {
					expect_error(t_running_sum(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_mean(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_sd(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_skew(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_kurt(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_sd3(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_skew4(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_kurt5(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_cent_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_cent_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,max_order_only=TRUE),NA)
					expect_error(t_running_std_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_cumulants(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_apx_quantiles(thingy,time=times,wts=wts,variable_win=TRUE,p=ptiles,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_apx_median(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_sharpe(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
					expect_error(t_running_sharpe(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,compute_se=TRUE),NA)
					expect_error(t_running_tstat(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)

					# NOTA BENE: these three do not accept an lb_time b/c they are associated with the thingy;
					# so we *do* expect these to error.
					expect_error(t_running_centered(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
					expect_error(t_running_scaled(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
					expect_error(t_running_zscored(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
				}

				for (lookahead in c(0,8.3)) {
					expect_error(t_running_centered(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
					expect_error(t_running_scaled(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
					expect_error(t_running_zscored(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
				}
			}
		}
	}

})#UNFOLD
context("code runs: t_running_foo odds and ends")
test_that("t_running foo  variable_win",{#FOLDUP
	skip_on_cran()

	set.char.seed("6e33dc8f-7a0b-4128-956f-f9ac960549ee")
	nel <- 25

	xna <- rnorm(2*nel)
	xna[xna < 0] <- NA

	xall <- list(rnorm(nel),
							 xna,
							 as.integer(rnorm(nel,sd=100)))

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	rp <- 5L
	wts_as_delta <- FALSE

	for (thingy in xall) {
		for (times in list(cumsum(runif(length(thingy),min=0.2,max=0.4)))) {
			wall <- list(rep(1.0,length(thingy)), NULL)
			for (wts in wall) {
				for (na_rm in c(FALSE,TRUE)) {
					for (lb_time in list(max(times) + c(1:2),3+cumsum(runif(10,min=0.4,max=1.1)))) {
						expect_error(t_running_sum(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_mean(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_sd(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_skew(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_kurt(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_sd3(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_skew4(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_kurt5(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_cent_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_cent_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,max_order_only=TRUE),NA)
						expect_error(t_running_std_moments(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_cumulants(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_apx_quantiles(thingy,time=times,wts=wts,variable_win=TRUE,p=ptiles,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_apx_median(thingy,time=times,wts=wts,variable_win=TRUE,max_order=5L,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_sharpe(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_sharpe(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm,compute_se=TRUE),NA)
						expect_error(t_running_tstat(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)

						# NOTA BENE: these three do not accept an lb_time b/c they are associated with the thingy;
						# so we *do* expect these to error.
						expect_error(t_running_centered(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
						expect_error(t_running_scaled(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
						expect_error(t_running_zscored(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm))
					}

					for (lookahead in c(0,8.3)) {
						expect_error(t_running_centered(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
						expect_error(t_running_scaled(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
						expect_error(t_running_zscored(thingy,time=times,wts=wts,variable_win=TRUE,wts_as_delta=wts_as_delta,restart_period=rp,lookahead=lookahead,na_rm=na_rm),NA)
					}
				}
			}
		}
	}


})#UNFOLD
test_that("t_running foo  infinite recompute",{#FOLDUP
	skip_on_cran()

	set.char.seed("b3efd3a7-9a46-4d37-b6a1-d864a5f9a0f0")
	nel <- 25

	xna <- rnorm(2*nel)
	xna[xna < 0] <- NA

	xall <- list(rnorm(nel),
							 xna,
							 as.integer(rnorm(nel,sd=100)))

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)
	rp <- 5L
	wts_as_delta <- FALSE

	for (thingy in xall) {
		for (times in list(cumsum(runif(length(thingy),min=0.2,max=0.4)))) {
			wall <- list(rep(1.0,length(thingy)), NULL)
			for (wts in wall) {
				for (na_rm in c(FALSE,TRUE)) {
					for (lb_time in list(max(times) + c(1:2),3+cumsum(runif(10,min=0.4,max=1.1)))) {
						expect_error(vers1 <- t_running_sd(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=rp,na_rm=na_rm),NA)
						expect_error(vers2 <- t_running_sd(thingy,time=times,wts=wts,variable_win=TRUE,lb_time=lb_time,wts_as_delta=wts_as_delta,restart_period=NA_integer_,na_rm=na_rm),NA)
						expect_equal(vers1,vers2)
					}
				}
			}
		}
	}


})#UNFOLD

context("t_running foo input params")
test_that("window as integer or double",{#FOLDUP
	skip_on_cran()

	set.char.seed("c0b02fb4-cd35-4403-8a13-950b5f2236f2")
	nel <- 40
	thingy <- rnorm(nel)
	times <- seq_along(thingy)
	iwin <- 50L
	dwin <- as.numeric(iwin)
	charwin <- 'window'

	expect_error(rd <- t_running_sum(thingy,time=times,window=dwin),NA)
	expect_error(ri <- t_running_sum(thingy,time=times,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(t_running_sum(thingy,times=times,window=charwin))
	
	expect_error(rd <- t_running_mean(thingy,time=times,window=dwin),NA)
	expect_error(ri <- t_running_mean(thingy,time=times,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(t_running_mean(thingy,times=times,window=charwin))

	expect_error(rd <- t_running_sd(thingy,time=times,window=dwin),NA)
	expect_error(ri <- t_running_sd(thingy,time=times,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(t_running_sd(thingy,times=times,window=charwin))

	expect_error(rd <- t_running_skew(thingy,time=times,window=dwin),NA)
	expect_error(ri <- t_running_skew(thingy,time=times,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(t_running_skew(thingy,times=times,window=charwin))

	expect_error(rd <- t_running_kurt(thingy,time=times,window=dwin),NA)
	expect_error(ri <- t_running_kurt(thingy,time=times,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(t_running_kurt(thingy,times=times,window=charwin))


})#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
