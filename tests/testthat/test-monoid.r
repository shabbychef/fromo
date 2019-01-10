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
#UNFOLD

context("monoid join/unjoin")
test_that("join/unjoin",{#FOLDUP
	set.char.seed("1325a51e-1584-4f89-9ea3-f15223a223d9")

	x1 <- rnorm(1e3,mean=1)
	x2 <- rnorm(1e3,mean=1)
	max_ord <- 6L
	expect_error(rs1 <- cent_sums(x1,max_ord),NA)
	expect_error(rs2 <- cent_sums(x2,max_ord),NA)
	expect_error(rs3 <- cent_sums(c(x1,x2),max_ord),NA)
	# make sure these don't change? 
	copy_rs1 <- rs1 + 0
	copy_rs2 <- rs2 + 0
	rs3alt <- join_cent_sums(rs1,rs2)
	expect_equal(rs1,copy_rs1,tolerance=1e-7)
	expect_equal(rs2,copy_rs2,tolerance=1e-7)
	expect_equal(rs3,rs3alt,tolerance=1e-7)

	copy_rs1 <- rs1 + 0
	copy_rs2 <- rs2 + 0
	copy_rs3 <- rs3 + 0

	rs1alt <- unjoin_cent_sums(rs3,rs2)
	rs2alt <- unjoin_cent_sums(rs3,rs1)
	expect_equal(rs1,copy_rs1,tolerance=1e-7)
	expect_equal(rs2,copy_rs2,tolerance=1e-7)
	expect_equal(rs3,copy_rs3,tolerance=1e-7)

	expect_equal(rs1,rs1alt,tolerance=1e-7)
	expect_equal(rs2,rs2alt,tolerance=1e-7)
})#UNFOLD
context("monoid cosum")
test_that("cosums are sane",{#FOLDUP
	set.char.seed("0020a8c0-ff6a-447c-a9bf-c6cc7160195f")

	x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	max_ord <- 2L
	expect_error(rs1 <- cent_comoments(x1,max_ord,used_df=1L),NA)
	expect_equal(rs1[1,1],nrow(x1))
	expect_equal(rs1[1,1 + (1:ncol(x1))],colMeans(x1),tolerance=1e-7)
	expect_equal(rs1[1 + (1:ncol(x1)),1],colMeans(x1),tolerance=1e-7)
	expect_equal(rs1[1 + (1:ncol(x1)),1 + (1:ncol(x1))],cov(x1),tolerance=1e-7)
})#UNFOLD
context("monoid join/unjoin cosums")
test_that("join/unjoin cosums",{#FOLDUP
	set.char.seed("9ecdda29-aaae-4f88-9fe7-4418846ca54c")

	x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	x2 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	max_ord <- 2L
	expect_error(rs1 <- cent_cosums(x1,max_ord),NA)
	expect_error(rs2 <- cent_cosums(x2,max_ord),NA)
	expect_error(rs3 <- cent_cosums(rbind(x1,x2),max_ord),NA)
	rs3alt <- join_cent_cosums(rs1,rs2)
	expect_lt(max(abs(rs3 - rs3alt)),1e-7)

	expect_error(rs1alt <- unjoin_cent_cosums(rs3,rs2),NA)
	expect_error(rs2alt <- unjoin_cent_cosums(rs3,rs1),NA)
	expect_lt(max(abs(rs1 - rs1alt)),1e-7)
	expect_lt(max(abs(rs2 - rs2alt)),1e-7)
})#UNFOLD

# 2FIX: check the effects of NA

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
