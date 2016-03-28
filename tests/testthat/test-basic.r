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

		sd3(y,na_rm=na_rm)
		skew4(y,na_rm=na_rm)
		kurt5(y,na_rm=na_rm)
		cent_moments(y,max_order=5L,used_df=1L,na_rm=na_rm)

		sd3(z,na_rm=na_rm)
		skew4(z,na_rm=na_rm)
		kurt5(z,na_rm=na_rm)
		cent_moments(z,max_order=5L,used_df=1L,na_rm=na_rm)
	}

	expect_error(sd3(q))
	expect_error(skew4(q))
	expect_error(kurt5(q))
	expect_error(cent_moments(q))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	x <- rnorm(300)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (winsize in c(100,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			run_sd3(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_skew4(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_kurt5(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_centered(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_scaled(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_zscored(x,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_tscored(x,winsize=winsize,recoper=50L,na_rm=na_rm)

			run_sd3(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_skew4(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_kurt5(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_centered(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_scaled(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_zscored(y,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_tscored(y,winsize=winsize,recoper=50L,na_rm=na_rm)

			run_sd3(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_skew4(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_kurt5(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_centered(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_scaled(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_zscored(z,winsize=winsize,recoper=50L,na_rm=na_rm)
			run_tscored(z,winsize=winsize,recoper=50L,na_rm=na_rm)

		}
	}
	expect_error(run_sd3(q))
	expect_error(run_skew4(q))
	expect_error(run_kurt5(q))
	expect_error(run_centered(q))
	expect_error(run_scaled(q))
	expect_error(run_zscored(q))
	expect_error(run_tscored(q))


	# sentinel
	expect_true(TRUE)
})#UNFOLD
# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
