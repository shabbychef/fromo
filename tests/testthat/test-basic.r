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
	sd3(x)
	skew4(x)
	kurt5(x)

	y <- as.integer(x)
	sd3(y)
	skew4(y)
	kurt5(y)

	z <- as.logical(y)
	sd3(z)
	skew4(z)
	kurt5(z)

	q <- c('a','b','c')
	expect_error(sd3(q))
	expect_error(skew4(q))
	expect_error(kurt5(q))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("running sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("7097f6ae-eac7-4e3a-b2cc-e9d4a01d43f7")
	x <- rnorm(200)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')

	for (winsize in c(20,Inf)) {
		for (na_rm in c(FALSE,TRUE)) {
			run_sd3(x,winsize=winsize,na_rm=na_rm)
			run_skew4(x,winsize=winsize,na_rm=na_rm)
			run_kurt5(x,winsize=winsize,na_rm=na_rm)

			run_sd3(y,winsize=winsize,na_rm=na_rm)
			run_skew4(y,winsize=winsize,na_rm=na_rm)
			run_kurt5(y,winsize=winsize,na_rm=na_rm)

			run_sd3(z,winsize=winsize,na_rm=na_rm)
			run_skew4(z,winsize=winsize,na_rm=na_rm)
			run_kurt5(z,winsize=winsize,na_rm=na_rm)
		}
	}

	expect_error(run_sd3(q))
	expect_error(run_skew4(q))
	expect_error(run_kurt5(q))

	# sentinel
	expect_true(TRUE)
})#UNFOLD
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
		# skew
		expect_equal(ske[1],skewness(x),tolerance=1e-9)
		# kurtosis
		expect_equal(krt[1],kurtosis(x) - 3.0,tolerance=1e-9)
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
