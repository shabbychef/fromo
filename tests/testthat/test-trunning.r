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
context("trunning_foo check heywood cases")
test_that("hit heywood branch",{#FOLDUP
	skip_on_cran()
	# not sure what I had planned here, and how it tests for heywoods.
	set.char.seed("3d318f1d-9921-4a20-84fc-c5ffc722d52c")
	xvals <- list(rnorm(1e5,mean=1e10))
	x <- rnorm(1e4)
	x[x < 1.0] <- NA
	xvals[[length(xvals)+1]] <- x
	x <- rnorm(1e4)
	x[x < 1.5] <- NA
	xvals[[length(xvals)+1]] <- x

	for (x in xvals) {
		window <- 500
		restart_period <- 100000L
		na_rm <- TRUE
		times <- seq_along(x)

		expect_error(y <- t_running_sd3(x,time=times,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- y[window:nrow(y),1]
		expect_false(any(ys < 0 | is.na(ys)))

		expect_error(y <- t_running_skew4(x,time=times,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- y[window:nrow(y),2]
		expect_false(any(ys < 0 | is.na(ys)))

		expect_error(y <- t_running_kurt5(x,time=times,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- 3 + y[window:nrow(y),1]
		expect_false(any(ys < 0 | is.na(ys)))
		ys <- y[window:nrow(y),3]
		expect_false(any(ys < 0 | is.na(ys)))
		# not sure where I would expect Heywoods in the other ones.
	}
})#UNFOLD

context("running x y code")
test_that("runs without error",{#FOLDUP
	skip_on_cran()

	set.char.seed("c13a49ac-14e9-462b-b81e-d3a5fc3be491")
	nel <- 20
	xna <- rnorm(nel)
	xna[xna < -0.5] <- NA
	xall <- list(rnorm(nel),
							 xna,
							 as.integer(rnorm(nel,sd=100)))

	wna <- runif(nel,min=1,max=3)
	wna[wna < 1.5] <- NA
	wall <- list(rep(1.0,nel),
							 runif(nel,min=0.9,max=3.5),
							 wna,
							 as.integer(ceiling(runif(nel,min=2,max=100))),
							 as.logical(ceiling(pmax(0,rnorm(nel)))),
							 NULL)

	for (x_thingy in xall) {
		y_thingy <- x_thingy + 1
		for (wts in wall) {
			for (window in c(5,21,Inf,NULL)) {
				for (na_rm in c(FALSE,TRUE)) {
					for (rp in c(1L,40L)) {
						times <- seq_along(x_thingy)
						expect_error(t_running_correlation(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_covariance(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_covariance_3(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_regression_slope(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_regression_intercept(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_regression_fit(x_thingy,y_thingy,wts=wts,time=times,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(t_running_regression_diagnostics(x_thingy,y_thingy,time=times,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
					}
				}
			}
		}
	}
})#UNFOLD
test_that("covariance correctness",{#FOLDUP
	skip_on_cran()

	set.char.seed("b81052f6-d7a5-4f3f-95cf-bc8d0030fcf8")
	window <- 30
	maxplus <- 50
	nel <- window + maxplus
	xvec <- rnorm(nel)
	beta_0 <- 0.33
	beta_1 <- 5
	sigma <- 0.5
	yvec <- beta_0 + beta_1 * xvec + rnorm(length(xvec),sd=sigma)
	times <- seq_along(xvec)
	expect_error(rho <- t_running_correlation(xvec,yvec,time=times,window=window-0.5),NA)
	expect_equal(rho[window], cor(xvec[1:window],yvec[1:window]), tolerance=1e-12)
	expect_equal(rho[1+window], cor(xvec[2:(1+window)],yvec[2:(1+window)]), tolerance=1e-12)

	expect_error(rho <- t_running_covariance(xvec,yvec,time=times,window=window-0.5,used_df=1),NA)
	expect_equal(rho[window], cov(xvec[1:window],yvec[1:window]), tolerance=1e-12)
	expect_equal(rho[1+window], cov(xvec[2:(1+window)],yvec[2:(1+window)]), tolerance=1e-12)

	expect_error(beta_1 <- t_running_regression_slope(xvec,yvec,time=times,window=window),NA)
	expect_error(beta_0 <- t_running_regression_intercept(xvec,yvec,time=times,window=window),NA)
	mod0 <- lm(yvec[1:window] ~ xvec[1:window])
	ses <- sqrt(diag(vcov(mod0)))

	expect_equal(beta_0[window], coefficients(mod0)[[1]], tolerance=1e-12)
	expect_equal(beta_1[window], coefficients(mod0)[[2]], tolerance=1e-12)
	expect_error(beta_ff <- t_running_regression_fit(xvec,yvec,time=times,window=window),NA)
	expect_equal(beta_ff[,1,drop=FALSE], beta_0, tolerance=1e-12)
	expect_equal(beta_ff[,2,drop=FALSE], beta_1, tolerance=1e-12)
	expect_error(beta_dd <- t_running_regression_diagnostics(xvec,yvec,time=times,window=window),NA)
	expect_equal(beta_dd[,1,drop=FALSE], beta_0, tolerance=1e-12)
	expect_equal(beta_dd[,2,drop=FALSE], beta_1, tolerance=1e-12)
	# also compare against diag(vcov(mod0))
	expect_equal(beta_dd[window,3], summary(mod0)$sigma, tolerance=1e-12)
	expect_equal(beta_dd[window,4], as.numeric(ses[1]), tolerance=1e-12)
	expect_equal(beta_dd[window,5], as.numeric(ses[2]), tolerance=1e-12)

	# and a little bit forward
	for (offs in c(5, maxplus - 1)) {
		mod1 <- lm(yvec[offs + (1:window)] ~ xvec[offs + (1:window)])
		ses <- sqrt(diag(vcov(mod1)))
		expect_equal(beta_0[offs+window], coefficients(mod1)[[1]], tolerance=1e-12)
		expect_equal(beta_1[offs+window], coefficients(mod1)[[2]], tolerance=1e-12)
		expect_equal(beta_dd[offs+window,3], summary(mod1)$sigma, tolerance=1e-12)
		expect_equal(beta_dd[offs+window,4], as.numeric(ses[1]), tolerance=1e-12)
		expect_equal(beta_dd[offs+window,5], as.numeric(ses[2]), tolerance=1e-12)
	}
})#UNFOLD
test_that("compare to running",{#FOLDUP
	# the idea here is that 
	skip_on_cran()

	set.char.seed("b81052f6-d7a5-4f3f-95cf-bc8d0030fcf8")
	window <- 30
	maxplus <- 50
	nel <- window + maxplus
	xvec <- rnorm(nel)
	beta_0 <- 0.33
	beta_1 <- 5
	sigma <- 0.5
	yvec <- beta_0 + beta_1 * xvec + rnorm(length(xvec),sd=sigma)
	times <- seq_along(xvec)
	expect_error(trho <- t_running_correlation(xvec,yvec,time=times,window=window),NA)
	expect_error(rho <- running_correlation(xvec,yvec,window=window,check_negative_moments=FALSE),NA)
	expect_equal(rho, trho, tolerance=1e-12)

	expect_error(tsxy <- t_running_covariance(xvec,yvec,time=times,window=window),NA)
	expect_error(sxy <- running_covariance(xvec,yvec,window=window,check_negative_moments=FALSE),NA)
	expect_equal(sxy, tsxy, tolerance=1e-12)

	expect_error(tbeta_0 <- t_running_regression_intercept(xvec,yvec,time=times,window=window),NA)
	expect_error(beta_0 <- running_regression_intercept(xvec,yvec,window=window,check_negative_moments=FALSE),NA)
	expect_equal(beta_0, tbeta_0, tolerance=1e-12)

	expect_error(tbeta_1 <- t_running_regression_slope(xvec,yvec,time=times,window=window),NA)
	expect_error(beta_1 <- running_regression_slope(xvec,yvec,window=window,check_negative_moments=FALSE),NA)
	expect_equal(beta_1, tbeta_1, tolerance=1e-12)

	expect_error(tbeta_dd <- t_running_regression_diagnostics(xvec,yvec,time=times,window=window),NA)
	expect_error(beta_dd <- running_regression_diagnostics(xvec,yvec,window=window,check_negative_moments=FALSE),NA)
	expect_equal(beta_dd, tbeta_dd, tolerance=1e-12)

})#UNFOLD
test_that("covariance weighting correctness",{#FOLDUP
	skip_on_cran()

	set.char.seed("f3d1652c-2d83-4670-998b-f96d18d18374")
	window <- 50
	nel <- window
	xvec <- rnorm(nel)
	beta_0 <- 0.33
	beta_1 <- 5
	sigma <- 0.5
	yvec <- beta_0 + beta_1 * xvec + rnorm(length(xvec),sd=sigma)
	wts <- sample(c(1,2,3),nel,replace=TRUE)
	times <- seq_along(xvec)
	expect_error(rho <- t_running_correlation(xvec,yvec,time=times,window=window,wts=wts),NA)
	bigx <- rep(xvec,wts)
	bigy <- rep(yvec,wts)
	bigtimes <- seq_along(bigx)
	expect_error(rho2 <- t_running_correlation(bigx,bigy,time=bigtimes,window=length(bigx)),NA)
	expect_equal(rho[window],rho2[length(rho2)], tolerance=1e-12)

	expect_error(rho <- t_running_covariance(xvec,yvec,time=times,wts=wts,window=window,normalize_wts=FALSE,used_df=1),NA)
	expect_error(rho2 <- t_running_covariance(bigx,bigy,time=bigtimes,window=length(bigx),used_df=1),NA)
	expect_equal(rho[window],rho2[length(rho2)], tolerance=1e-12)

	expect_error(beta_0 <- t_running_regression_intercept(xvec,yvec,time=times,wts=wts,window=window),NA)
	expect_error(beta_1 <- t_running_regression_slope(xvec,yvec,time=times,wts=wts,window=window),NA)
	expect_error(beta_02 <- t_running_regression_intercept(bigx,bigy,time=bigtimes,window=length(bigx)),NA)
	expect_error(beta_12 <- t_running_regression_slope(bigx,bigy,time=bigtimes,window=length(bigx)),NA)
	expect_equal(beta_0[window],beta_02[length(bigx)], tolerance=1e-12)
	expect_equal(beta_1[window],beta_12[length(bigx)], tolerance=1e-12)

	expect_error(beta_dd <- t_running_regression_diagnostics(xvec,yvec,time=times,wts=wts,window=window,normalize_wts=FALSE),NA)
	expect_error(beta_dd2 <- t_running_regression_diagnostics(bigx,bigy,time=bigtimes,window=length(bigx)),NA)
	expect_equal(beta_dd[window,,drop=TRUE],beta_dd2[length(bigx),,drop=TRUE], tolerance=1e-12)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
