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
# Copyright: Steven E. Pav, 2016-2019
# Author: Steven E. Pav
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
#UNFOLD

context("running_foo run without error")
test_that("running sd, skew, kurt run without error",{#FOLDUP
	skip_on_cran()

	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (thingy in list(x,y,z)) { 
		for (window in c(50,Inf)) {
			for (na_rm in c(FALSE,TRUE)) {
				expect_error(running_sum(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_mean(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_sd(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_skew(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_kurt(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_sd3(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_skew4(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_kurt5(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				for (mol in c(1L,2L,5L)) {
					expect_error(running_cent_moments(thingy,max_order=mol,window=window,restart_period=50L,na_rm=na_rm),NA)
					expect_error(running_cent_moments(thingy,max_order=mol,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
				}
				expect_error(running_std_moments(thingy,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_cumulants(thingy,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_apx_quantiles(thingy,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
			}
		}
	}
	for (thingy in list(x,y,z)) { 
		for (min_df in c(2L,10L)) {
			expect_error(running_mean(thingy,window=window,min_df=min_df),NA)
			expect_error(running_sd(thingy,window=window,min_df=min_df),NA)
			expect_error(running_skew(thingy,window=window,min_df=min_df),NA)
			expect_error(running_kurt(thingy,window=window,min_df=min_df),NA)
			expect_error(running_sd3(thingy,window=window,min_df=min_df),NA)
			expect_error(running_skew4(thingy,window=window,min_df=min_df),NA)
			expect_error(running_kurt5(thingy,window=window,min_df=min_df),NA)
			for (mol in c(1L,2L,5L)) {
				expect_error(running_cent_moments(thingy,max_order=mol,window=window,min_df=min_df),NA)
				expect_error(running_cent_moments(thingy,max_order=mol,window=window,min_df=min_df,max_order_only=TRUE),NA)
			}
			expect_error(running_std_moments(thingy,max_order=5L,window=window,min_df=min_df),NA)
			expect_error(running_cumulants(thingy,max_order=5L,window=window,min_df=min_df),NA)
			expect_error(running_apx_quantiles(thingy,p=ptiles,max_order=5L,window=window,min_df=min_df),NA)
		}
	}

	# will not work on character input
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
	# this crashes R when I run it:
	#expect_error(running_apx_quantiles(x,p=q,max_order=5L))
	expect_error(running_apx_median(q,p=ptiles,max_order=5L))
})#UNFOLD
context("running_foo weighted OK")
test_that("running foo and weights",{#FOLDUP
	skip_on_cran()

	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
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
							 NULL)

	ptiles <- c(0.1,0.25,0.5,0.75,0.9)

	for (thingy in xall) {
		for (wts in wall) {
			for (window in c(5,21,Inf,NULL)) {
				for (na_rm in c(FALSE,TRUE)) {
						expect_error(running_sum(thingy,wts=wts,window=window,restart_period=2L,na_rm=na_rm),NA)

						expect_error(running_sum(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_mean(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_sd(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_skew(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_kurt(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_sd3(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_skew4(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_kurt5(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						for (mol in c(1L,2L,5L)) {
							expect_error(running_cent_moments(thingy,wts=wts,max_order=mol,window=window,restart_period=50L,na_rm=na_rm),NA)
							expect_error(running_cent_moments(thingy,wts=wts,max_order=mol,window=window,restart_period=50L,na_rm=na_rm,max_order_only=TRUE),NA)
						}
						expect_error(running_std_moments(thingy,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_cumulants(thingy,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_apx_quantiles(thingy,wts=wts,p=ptiles,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_apx_median(thingy,wts=wts,max_order=5L,window=window,restart_period=50L,na_rm=na_rm),NA)

						#expect_error(running_centered(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)

						for (lookahead in c(0L,8L)) {
							expect_error(running_centered(thingy,wts=wts,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
							expect_error(running_scaled(thingy,wts=wts,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
							expect_error(running_zscored(thingy,wts=wts,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
						}
						expect_error(running_sharpe(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
						expect_error(running_sharpe(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
						expect_error(running_tstat(thingy,wts=wts,window=window,restart_period=50L,na_rm=na_rm),NA)
				}
			}
		}
	}


})#UNFOLD

context("running_sum every way")
test_that("running sum",{#FOLDUP
	skip_on_cran()

	set.char.seed("20c032c9-deb2-4000-9d73-b7d8395b67b2")
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

	for (thingy in xall) {
		for (wts in wall) {
			for (window in c(5,21,Inf,NULL)) {
				for (na_rm in c(FALSE,TRUE)) {
					for (rp in c(1L,40L)) {
						expect_error(running_sum(thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
					}
				}
			}
		}
	}
})#UNFOLD
test_that("constant input",{#FOLDUP
	x <- rep(1,8)
	for (maxord in c(2,4,5)) {
		expect_error(yy <- running_cent_moments(x,window=6,max_order=maxord),NA)
		expect_true(!any(is.na(yy[nrow(yy),,drop=TRUE]) ))
	}
	expect_error(y3 <- running_sd3(x,window=5),NA)
	expect_true(!is.na(y3[nrow(y3),1,drop=TRUE]))
	expect_error(y4 <- running_skew4(x,window=5),NA)
	expect_true(!is.na(y4[nrow(y4),2,drop=TRUE]))
	expect_error(y5 <- running_kurt5(x,window=5),NA)
	expect_true(!is.na(y5[nrow(y5),3,drop=TRUE]))
	# if skewness is divided by stdev, there will be 0/0 condition here, 
	# so drop these
	# expect_error(y4 <- running_skew4(x,window=5),NA)
	# expect_true(!any(is.na(y4[nrow(y4),,drop=TRUE]) ))
	# expect_error(y5 <- running_kurt5(x,window=5),NA)
	# expect_true(!any(is.na(y5[nrow(y5),,drop=TRUE]) ))
})#UNFOLD

context("running foo input params")
test_that("window as integer or double",{#FOLDUP
	skip_on_cran()

	set.char.seed("da774aca-bbde-4350-87b2-d21ed0c84124")
	nel <- 40
	thingy <- rnorm(nel)
	iwin <- 50L
	dwin <- as.numeric(iwin)
	charwin <- 'window'

	expect_error(rd <- running_sum(thingy,window=dwin),NA)
	expect_error(ri <- running_sum(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_sum(thingy,window=charwin))

	expect_error(rd <- running_mean(thingy,window=dwin),NA)
	expect_error(ri <- running_mean(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_mean(thingy,window=charwin))

	expect_error(rd <- running_sd(thingy,window=dwin),NA)
	expect_error(ri <- running_sd(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_sd(thingy,window=charwin))

	expect_error(rd <- running_skew(thingy,window=dwin),NA)
	expect_error(ri <- running_skew(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_skew(thingy,window=charwin))

	expect_error(rd <- running_kurt(thingy,window=dwin),NA)
	expect_error(ri <- running_kurt(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_kurt(thingy,window=charwin))

	expect_error(rd <- running_sd3(thingy,window=dwin),NA)
	expect_error(ri <- running_sd3(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_sd3(thingy,window=charwin))

	expect_error(rd <- running_skew4(thingy,window=dwin),NA)
	expect_error(ri <- running_skew4(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_skew4(thingy,window=charwin))

	expect_error(rd <- running_kurt5(thingy,window=dwin),NA)
	expect_error(ri <- running_kurt5(thingy,window=iwin),NA)
	expect_equal(rd,ri,tolerance=1e-12)
	expect_error(running_kurt5(thingy,window=charwin))
})#UNFOLD
test_that("bad input",{#FOLDUP
	skip_on_cran()

	set.char.seed("e6bc7b67-8ba5-4c48-aeed-180a27d3303c")
	nel <- 10
	thingy <- rnorm(nel)
	# for some reason, these now cause R to abort? auuggggh!

	expect_error(rd <- running_sum(thingy,window='bad idea'))
	expect_error(rd <- running_mean(thingy,window=5,restart_period='dumb'))
	expect_error(rd <- running_mean(thingy,window=5,min_df='dumb'))
	expect_error(rd <- running_mean(thingy,window=5,na_rm='dumb'))
	expect_error(rd <- running_sd3(thingy,window=5,used_df='dumb'))
	expect_error(rd <- running_sd3(thingy,window=5,check_wts='dumb'))
	expect_error(rd <- running_sd3(thingy,window=5,normalize_wts='dumb'))
})#UNFOLD
context("running_foo check heywood cases")
test_that("hit heywood branch",{#FOLDUP
	skip_on_cran()

	# not sure what I had planned here, and how it tests for heywoods.
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

		expect_error(y <- running_sd3(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- y[window:nrow(y),1]
		expect_false(any(ys < 0 | is.na(ys)))

		expect_error(y <- running_skew4(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- y[window:nrow(y),2]
		expect_false(any(ys < 0 | is.na(ys)))

		expect_error(y <- running_kurt5(x,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		ys <- 3 + y[window:nrow(y),1]
		expect_false(any(ys < 0 | is.na(ys)))
		ys <- y[window:nrow(y),3]
		expect_false(any(ys < 0 | is.na(ys)))

		# not sure where I would expect Heywoods in here.
		# expect_error(y <- running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		# expect_error(running_cent_moments(x,max_order=5L,window=window,restart_period=restart_period,max_order_only=TRUE,na_rm=na_rm),NA)
		# expect_error(running_std_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		# expect_error(running_cumulants(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		# expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		# expect_error(running_apx_median(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
	}
})#UNFOLD
test_that("hit corner cases",{#FOLDUP
	# created with the intention of causing Heywoods
	epsi <- 1e-16
	xx <- c(1,1+epsi,epsi,1,1,1)
	expect_error(y <- fromo::running_sd3(xx,3,check_negative_moments=FALSE),NA)
	expect_true(is.nan(y[nrow(y),1]))
	expect_error(y <- fromo::running_sd3(xx,3,check_negative_moments=TRUE),NA)
	expect_false(is.nan(y[nrow(y),1]))
})#UNFOLD
context("running_foo: NA restart")
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
		for (mol in c(1L,2L,5L)) {
			expect_error(running_cent_moments(x,max_order=mol,window=window,restart_period=restart_period,max_order_only=TRUE,na_rm=na_rm),NA)
		}
		expect_error(running_std_moments(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_cumulants(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_quantiles(x,p=ptiles,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)
		expect_error(running_apx_median(x,max_order=5L,window=window,restart_period=restart_period,na_rm=na_rm),NA)

	}
})#UNFOLD
context("running foo NANANANANA")
test_that("nananana",{#FOLDUP
	skip_on_cran()

	x <- rep(NA,20)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')
	restart_period <- NA_integer_

	for (thingy in list(x,y,z)) { 
		for (wts in list(NULL,rep(NA_real_,length(thingy)))) {
			for (window in c(15)) {
				for (na_rm in c(FALSE,TRUE)) {
					expect_error(running_sum(thingy,window=window,wts=wts,restart_period=restart_period,na_rm=na_rm),NA)
					expect_error(running_mean(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm),NA)

					for (nw in c(FALSE,TRUE)) { 
						expect_error(running_sd(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_skew(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_kurt(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)

						expect_error(running_sd3(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_skew4(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_kurt5(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_cent_moments(thingy,max_order=5L,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)

						expect_error(running_centered(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)

						expect_error(running_scaled(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_zscored(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_sharpe(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_sharpe(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,compute_se=TRUE,na_rm=na_rm,normalize_wts=nw),NA)
						expect_error(running_tstat(thingy,window=window,min_df=10L,wts=wts,restart_period=restart_period,na_rm=na_rm,normalize_wts=nw),NA)
					}

				}
			}
		}
	}
})#UNFOLD
context("running adjustments run")
test_that("running adjustments",{#FOLDUP
	skip_on_cran()

	set.char.seed("e36291f6-70e0-4412-9c50-bc46b6ab8639")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (thingy in list(x,y,z)) { 
		for (window in c(50,Inf)) {
			for (na_rm in c(FALSE,TRUE)) {
				expect_error(running_centered(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_scaled(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_zscored(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_sharpe(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
				expect_error(running_sharpe(thingy,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
				expect_error(running_tstat(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
			}
		}
	}

	window <- 10L
	for (thingy in list(x,y,z)) { 
		for (na_rm in c(FALSE,TRUE)) {
			expect_error(running_centered(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_scaled(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_zscored(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
			expect_error(running_sharpe(thingy,window=window,restart_period=50L,na_rm=na_rm,compute_se=TRUE),NA)
			expect_error(running_tstat(thingy,window=window,restart_period=50L,na_rm=na_rm),NA)
		}
	}

	window <- 10L
	for (thingy in list(x,y,z)) { 
		for (na_rm in c(FALSE,TRUE)) {
			for (lookahead in c(0,-5,10)) {
				expect_error(running_centered(thingy,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
				expect_error(running_scaled(thingy,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
				expect_error(running_zscored(thingy,window=window,restart_period=50L,lookahead=lookahead,na_rm=na_rm),NA)
			}
		}
	}

	for (thingy in list(x,y,z)) { 
		for (min_df in c(2L,10L)) {
			expect_error(running_centered(thingy,window=window,min_df=min_df),NA)
			expect_error(running_scaled(thingy,window=window,min_df=min_df),NA)
			expect_error(running_zscored(thingy,window=window,min_df=min_df),NA)
			expect_error(running_sharpe(thingy,window=window,min_df=min_df),NA)
			expect_error(running_sharpe(thingy,window=window,min_df=min_df,compute_se=TRUE),NA)
			expect_error(running_tstat(thingy,window=window,min_df=min_df),NA)
		}
	}

	expect_error(running_centered(q))
	expect_error(running_scaled(q))
	expect_error(running_zscored(q))
	expect_error(running_sharpe(q))
	expect_error(running_sharpe(q,compute_se=TRUE))
	expect_error(running_tstat(q))

	#expect_error(running_tstat(x,window='FOO'))
	#expect_error(running_tstat(x,window=-20L))
	#expect_error(running_tstat(x,window=20L,restart_period='FOO'))
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
						expect_error(running_correlation(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_covariance(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_covariance_3(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_regression_slope(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_regression_intercept(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_regression_fit(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
						expect_error(running_regression_diagnostics(x_thingy,y_thingy,wts=wts,window=window,restart_period=rp,na_rm=na_rm),NA)
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
	expect_error(rho <- running_correlation(xvec,yvec,window=window),NA)
	expect_equal(rho[window], cor(xvec[1:window],yvec[1:window]), tolerance=1e-12)
	expect_equal(rho[1+window], cor(xvec[2:(1+window)],yvec[2:(1+window)]), tolerance=1e-12)

	expect_error(rho <- running_covariance(xvec,yvec,window=window,used_df=1),NA)
	expect_equal(rho[window], cov(xvec[1:window],yvec[1:window]), tolerance=1e-12)
	expect_equal(rho[1+window], cov(xvec[2:(1+window)],yvec[2:(1+window)]), tolerance=1e-12)

	expect_error(beta_1 <- running_regression_slope(xvec,yvec,window=window),NA)
	expect_error(beta_0 <- running_regression_intercept(xvec,yvec,window=window),NA)
	mod0 <- lm(yvec[1:window] ~ xvec[1:window])
	ses <- sqrt(diag(vcov(mod0)))

	expect_equal(beta_0[window], coefficients(mod0)[[1]], tolerance=1e-12)
	expect_equal(beta_1[window], coefficients(mod0)[[2]], tolerance=1e-12)
	expect_error(beta_ff <- running_regression_fit(xvec,yvec,window=window),NA)
	expect_equal(beta_ff[,1,drop=FALSE], beta_0, tolerance=1e-12)
	expect_equal(beta_ff[,2,drop=FALSE], beta_1, tolerance=1e-12)
	expect_error(beta_dd <- running_regression_diagnostics(xvec,yvec,window=window),NA)
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
	expect_error(rho <- running_correlation(xvec,yvec,wts=wts,window=window),NA)
	bigx <- rep(xvec,wts)
	bigy <- rep(yvec,wts)
	expect_error(rho2 <- running_correlation(bigx,bigy,window=length(bigx)),NA)
	expect_equal(rho[window],rho2[length(rho2)], tolerance=1e-12)

	expect_error(rho <- running_covariance(xvec,yvec,wts=wts,window=window,normalize_wts=FALSE,used_df=1),NA)
	expect_error(rho2 <- running_covariance(bigx,bigy,window=length(bigx),used_df=1),NA)
	expect_equal(rho[window],rho2[length(rho2)], tolerance=1e-12)

	expect_error(beta_0 <- running_regression_intercept(xvec,yvec,wts=wts,window=window),NA)
	expect_error(beta_1 <- running_regression_slope(xvec,yvec,wts=wts,window=window),NA)
	expect_error(beta_02 <- running_regression_intercept(bigx,bigy,window=length(bigx)),NA)
	expect_error(beta_12 <- running_regression_slope(bigx,bigy,window=length(bigx)),NA)
	expect_equal(beta_0[window],beta_02[length(bigx)], tolerance=1e-12)
	expect_equal(beta_1[window],beta_12[length(bigx)], tolerance=1e-12)

	expect_error(beta_dd <- running_regression_diagnostics(xvec,yvec,wts=wts,window=window,normalize_wts=FALSE),NA)
	expect_error(beta_dd2 <- running_regression_diagnostics(bigx,bigy,window=length(bigx)),NA)
	expect_equal(beta_dd[window,,drop=TRUE],beta_dd2[length(bigx),,drop=TRUE], tolerance=1e-12)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
