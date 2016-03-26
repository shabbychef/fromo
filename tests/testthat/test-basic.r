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
test_that("they run",{#FOLDUP
	set.char.seed("569dd47d-f9e5-40e4-b2ac-e5dbb4771a53")
	x <- rnorm(100)
	sd3(x)
	skew4(x)
	kurt5(x)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
test_that("sanity",{#FOLDUP
	set.char.seed("c4007dba-2010-481e-abe5-f07d3ce94eb4")
	x <- rnorm(1000)

	sid <- sd3(x)
	ske <- skew4(x)
	krt <- kurt5(x)

	expect_equal(length(sid),3)
	expect_equal(length(ske),4)
	expect_equal(length(krt),5)
	# length
	expect_equal(length(x),sid[3])
	expect_equal(sid[3],ske[4])
	expect_equal(ske[4],krt[5])

	# compare computations to gold standard
	expect_equal(sid[2],mean(x),tolerance=1e-9)
	expect_equal(sid[2],ske[3],tolerance=1e-9)
	expect_equal(sid[2],krt[4],tolerance=1e-9)

	expect_equal(sid[1],ske[2],tolerance=1e-9)
	expect_equal(sid[1],krt[3],tolerance=1e-9)
	expect_equal(sid[1],sd(x),tolerance=1e-9)

	expect_equal(ske[1],krt[2],tolerance=1e-9)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
