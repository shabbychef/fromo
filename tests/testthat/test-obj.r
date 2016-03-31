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
# Created: 2016.03.31
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
test_that("constructor and such",{#FOLDUP
	set.char.seed("33c133f6-930f-4656-88d4-84493784eee3")

	x <- rnorm(100)
	csums <- cent_sums(x, max_order=5L, na_rm=TRUE)
	xobj <- centsums(csums)

	foo <- sums(xobj)
	for (type in c('central','raw','standardized')) {
		moms <- moments(xobj,type=type)
	}
	show(xobj)

	xobj <- as.centsums(x,order=5L, na.rm=TRUE)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# 2FIX: check the effects of NA
#UNFOLD
context("correctness")#FOLDUP
test_that("monoidal homomorphism",{#FOLDUP
	set.char.seed("3fb6ab54-f5b9-4965-8a76-477dd62d6d7e")

	x <- rnorm(100)
	y <- rnorm(100)
	order <- 5L
	xobj <- as.centsums(x,order=order, na.rm=TRUE)
	yobj <- as.centsums(y,order=order, na.rm=TRUE)
	zobj <- as.centsums(c(x,y),order=order, na.rm=TRUE)

	zalt <- c(xobj,yobj)
	expect_equal(sums(zalt),sums(zobj),tolerance=1e-9)

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
