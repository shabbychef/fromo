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

context("code runs at all")#FOLDUP
test_that("sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("569dd47d-f9e5-40e4-b2ac-e5dbb4771a53")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (na_rm in c(FALSE,TRUE)) {
		sd3(x,na_rm=na_rm)
		skew4(x,na_rm=na_rm)
		kurt5(x,na_rm=na_rm)
		cent_moments(x,max_order=5L,used_df=1L,na_rm=na_rm)
		std_moments(x,max_order=5L,used_df=1L,na_rm=na_rm)

		sd3(y,na_rm=na_rm)
		skew4(y,na_rm=na_rm)
		kurt5(y,na_rm=na_rm)
		cent_moments(y,max_order=5L,used_df=1L,na_rm=na_rm)
		std_moments(y,max_order=5L,used_df=1L,na_rm=na_rm)

		sd3(z,na_rm=na_rm)
		skew4(z,na_rm=na_rm)
		kurt5(z,na_rm=na_rm)
		cent_moments(z,max_order=5L,used_df=1L,na_rm=na_rm)
		std_moments(z,max_order=5L,used_df=1L,na_rm=na_rm)
	}

	expect_error(sd3(q))
	expect_error(skew4(q))
	expect_error(kurt5(q))
	expect_error(cent_moments(q))
	expect_error(std_moments(q))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (window in c(50,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			running_sd3(x,window=window,restart_period=50L,na_rm=na_rm)
			running_skew4(x,window=window,restart_period=50L,na_rm=na_rm)
			running_kurt5(x,window=window,restart_period=50L,na_rm=na_rm)
			running_cent_moments(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
			running_std_moments(x,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)

			running_sd3(y,window=window,restart_period=50L,na_rm=na_rm)
			running_skew4(y,window=window,restart_period=50L,na_rm=na_rm)
			running_kurt5(y,window=window,restart_period=50L,na_rm=na_rm)
			running_cent_moments(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
			running_std_moments(y,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)

			running_sd3(z,window=window,restart_period=50L,na_rm=na_rm)
			running_skew4(z,window=window,restart_period=50L,na_rm=na_rm)
			running_kurt5(z,window=window,restart_period=50L,na_rm=na_rm)
			running_cent_moments(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
			running_std_moments(z,max_order=5L,window=window,restart_period=50L,na_rm=na_rm)
		}
	}
	for (min_df in c(2L,10L)) {
		running_sd3(x,window=window,min_df=min_df)
		running_skew4(x,window=window,min_df=min_df)
		running_kurt5(x,window=window,min_df=min_df)
		running_cent_moments(x,max_order=5L,window=window,min_df=min_df)
		running_std_moments(x,max_order=5L,window=window,min_df=min_df)

		running_sd3(y,window=window,min_df=min_df)
		running_skew4(y,window=window,min_df=min_df)
		running_kurt5(y,window=window,min_df=min_df)
		running_cent_moments(y,max_order=5L,window=window,min_df=min_df)
		running_std_moments(y,max_order=5L,window=window,min_df=min_df)

		running_sd3(z,window=window,min_df=min_df)
		running_skew4(z,window=window,min_df=min_df)
		running_kurt5(z,window=window,min_df=min_df)
		running_cent_moments(z,max_order=5L,window=window,min_df=min_df)
		running_std_moments(z,max_order=5L,window=window,min_df=min_df)
	}

	expect_error(running_sd3(q))
	expect_error(running_skew4(q))
	expect_error(running_kurt5(q))
	expect_error(running_cent_moments(q,max_order=5L))
	expect_error(running_std_moments(q,max_order=5L))

	# make sure the Heywood branch gets hit
	x <- rnorm(1e5,mean=1e10)
	window <- 500L
	restart_period <- 100000L
	running_sd3(x,window=window,restart_period=restart_period)
	running_skew4(x,window=window,restart_period=restart_period)
	running_kurt5(x,window=window,restart_period=restart_period)
	running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period)
	running_std_moments(x,max_order=5L,window=window,restart_period=restart_period)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running adjustments",{#FOLDUP
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
			running_tstat(x,window=window,restart_period=50L,na_rm=na_rm)

			running_centered(y,window=window,restart_period=50L,na_rm=na_rm)
			running_scaled(y,window=window,restart_period=50L,na_rm=na_rm)
			running_zscored(y,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm)
			running_tstat(y,window=window,restart_period=50L,na_rm=na_rm)

			running_centered(z,window=window,restart_period=50L,na_rm=na_rm)
			running_scaled(z,window=window,restart_period=50L,na_rm=na_rm)
			running_zscored(z,window=window,restart_period=50L,na_rm=na_rm)
			running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm)
			running_tstat(z,window=window,restart_period=50L,na_rm=na_rm)
		}
	}
	window <- 10L

	for (na_rm in c(FALSE,TRUE)) {
		running_centered(x,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(x,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(x,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(x,window=window,restart_period=50L,na_rm=na_rm)
		running_tstat(x,window=window,restart_period=50L,na_rm=na_rm)

		running_centered(y,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(y,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(y,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(y,window=window,restart_period=50L,na_rm=na_rm)
		running_tstat(y,window=window,restart_period=50L,na_rm=na_rm)

		running_centered(z,window=window,restart_period=50L,na_rm=na_rm)
		running_scaled(z,window=window,restart_period=50L,na_rm=na_rm)
		running_zscored(z,window=window,restart_period=50L,na_rm=na_rm)
		running_sharpe(z,window=window,restart_period=50L,na_rm=na_rm)
		running_tstat(z,window=window,restart_period=50L,na_rm=na_rm)
	}

	for (min_df in c(2L,10L)) {
		running_centered(x,window=window,min_df=min_df)
		running_scaled(x,window=window,min_df=min_df)
		running_zscored(x,window=window,min_df=min_df)
		running_sharpe(x,window=window,min_df=min_df)
		running_tstat(x,window=window,min_df=min_df)

		running_centered(y,window=window,min_df=min_df)
		running_scaled(y,window=window,min_df=min_df)
		running_zscored(y,window=window,min_df=min_df)
		running_sharpe(y,window=window,min_df=min_df)
		running_tstat(y,window=window,min_df=min_df)

		running_centered(z,window=window,min_df=min_df)
		running_scaled(z,window=window,min_df=min_df)
		running_zscored(z,window=window,min_df=min_df)
		running_sharpe(z,window=window,min_df=min_df)
		running_tstat(z,window=window,min_df=min_df)
	}

	expect_error(running_centered(q))
	expect_error(running_scaled(q))
	expect_error(running_zscored(q))
	expect_error(running_sharpe(q))
	expect_error(running_tstat(q))

	expect_error(running_tstat(x,window='FOO'))
	expect_error(running_tstat(x,window=-20L))
	expect_error(running_tstat(x,window=20L,restart_period='FOO'))

	# sentinel
	expect_true(TRUE)
})#UNFOLD

# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
