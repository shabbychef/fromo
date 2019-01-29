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
	for (max_ord in c(1L,2L,3L,6L)) {
		for (narm in c(FALSE,TRUE)) { 
			expect_error(rs1 <- cent_sums(x1,max_ord,na_rm=narm),NA)
			expect_equal(length(rs1),max_ord+1)
			expect_error(rs2 <- cent_sums(x2,max_ord,na_rm=narm),NA)
			expect_equal(length(rs2),max_ord+1)
			expect_error(rs3 <- cent_sums(c(x1,x2),max_ord,na_rm=narm),NA)
			# make sure these don't change? 
			copy_rs1 <- rs1 + 0
			copy_rs2 <- rs2 + 0
			expect_error(rs3alt <- join_cent_sums(rs1,rs2),NA)
			expect_equal(rs1,copy_rs1,tolerance=1e-7)
			expect_equal(rs2,copy_rs2,tolerance=1e-7)
			expect_equal(rs3,rs3alt,tolerance=1e-7)

			copy_rs1 <- rs1 + 0
			copy_rs2 <- rs2 + 0
			copy_rs3 <- rs3 + 0

			expect_error(rs1alt <- unjoin_cent_sums(rs3,rs2),NA)
			expect_error(rs2alt <- unjoin_cent_sums(rs3,rs1),NA)
			expect_equal(rs1,copy_rs1,tolerance=1e-7)
			expect_equal(rs2,copy_rs2,tolerance=1e-7)
			expect_equal(rs3,copy_rs3,tolerance=1e-7)

			expect_equal(rs1,rs1alt,tolerance=1e-7)
			expect_equal(rs2,rs2alt,tolerance=1e-7)

			# now an empty guy; this should return empty.
			x0 <- c()
			expect_error(rs0 <- cent_sums(x0,max_ord,na_rm=narm),NA)
			expect_equal(length(rs0),max_ord+1)
			expect_error(rs1 <- cent_sums(x1,max_ord,na_rm=narm),NA)
			copy_rs1 <- rs1 + 0
			expect_error(rs1alt <- join_cent_sums(rs1,rs0),NA)
			expect_equal(rs1,rs1alt,tolerance=1e-7)
		}
	}
})#UNFOLD
test_that("weighted join/unjoin",{#FOLDUP
	set.char.seed("1100eb6f-4108-414b-94d8-e31524e93461")

	x1 <- rnorm(1e3,mean=1)
	w1 <- runif(length(x1),min=0.5,max=1.5)
	x2 <- rnorm(1e3,mean=1)
	w2 <- runif(length(x2),min=0.5,max=1.5)
	for (max_ord in c(1L,2L,3L,6L)) {
		expect_error(rs1 <- cent_sums(x1,max_ord,wts=w1,normalize_wts=FALSE),NA)
		expect_equal(length(rs1),max_ord+1)
		expect_error(rs2 <- cent_sums(x2,max_ord,wts=w2,normalize_wts=FALSE),NA)
		expect_equal(length(rs2),max_ord+1)
		expect_error(rs3 <- cent_sums(c(x1,x2),max_ord,wts=c(w1,w2),normalize_wts=FALSE),NA)
		# make sure these don't change? 
		copy_rs1 <- rs1 + 0
		copy_rs2 <- rs2 + 0
		expect_error(rs3alt <- join_cent_sums(rs1,rs2),NA)
		expect_equal(rs1,copy_rs1,tolerance=1e-7)
		expect_equal(rs2,copy_rs2,tolerance=1e-7)
		expect_equal(rs3,rs3alt,tolerance=1e-7)

		copy_rs1 <- rs1 + 0
		copy_rs2 <- rs2 + 0
		copy_rs3 <- rs3 + 0

		expect_error(rs1alt <- unjoin_cent_sums(rs3,rs2),NA)
		expect_error(rs2alt <- unjoin_cent_sums(rs3,rs1),NA)
		expect_equal(rs1,copy_rs1,tolerance=1e-7)
		expect_equal(rs2,copy_rs2,tolerance=1e-7)
		expect_equal(rs3,copy_rs3,tolerance=1e-7)

		expect_equal(rs1,rs1alt,tolerance=1e-7)
		expect_equal(rs2,rs2alt,tolerance=1e-7)

		# now an empty guy; this should return empty.
		x0 <- c()
		expect_error(rs0 <- cent_sums(x0,max_ord,wts=c()),NA)
		expect_equal(length(rs0),max_ord+1)
		expect_error(rs1 <- cent_sums(x1,max_ord,wts=w1,normalize_wts=FALSE),NA)
		copy_rs1 <- rs1 + 0
		expect_error(rs1alt <- join_cent_sums(rs1,rs0),NA)
		expect_equal(rs1,rs1alt,tolerance=1e-7)
	}
	# just hit code to make sure it runs
	for (max_ord in c(1L,2L,3L,6L)) {
		for (narm in c(FALSE,TRUE)) {
			for (cw in c(FALSE,TRUE)) {
				for (nw in c(FALSE,TRUE)) {
					expect_error(rs1 <- cent_sums(x1,max_ord,wts=w1,na_rm=narm,
																				check_wts=cw,normalize_wts=nw),NA)
					expect_equal(length(rs1),max_ord+1)
				}
			}
		}
	}

})#UNFOLD
test_that("join commutativity",{#FOLDUP
	set.char.seed("f33946b3-216e-4977-9535-447b55214197")

	for (max_ord in c(1L,2L,3L,6L)) {
		x1 <- rnorm(1e3,mean=1)
		x2 <- rnorm(1e3,mean=1)
		expect_error(rs1 <- cent_sums(x1,max_ord),NA)
		expect_error(rs2 <- cent_sums(x2,max_ord),NA)
		expect_error(rs3_a <- join_cent_sums(rs1,rs2),NA)
		expect_error(rs3_b <- join_cent_sums(rs2,rs1),NA)
		expect_equal(rs3_a,rs3_b,tolerance=1e-7)

		# on an empty
		x1 <- c()
		x2 <- rnorm(1e3,mean=1)
		expect_error(rs1 <- cent_sums(x1,max_ord),NA)
		expect_error(rs2 <- cent_sums(x2,max_ord),NA)
		expect_error(rs3_a <- join_cent_sums(rs1,rs2),NA)
		expect_error(rs3_b <- join_cent_sums(rs2,rs1),NA)
		expect_equal(rs3_a,rs3_b,tolerance=1e-7)
	}

	# should this be an error?
	#rs1[1] <- -1
	#expect_error(rs3_a <- join_cent_sums(rs1,rs2))
	#expect_error(rs3_b <- join_cent_sums(rs2,rs1))
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
test_that("cosums and na_omit ",{#FOLDUP
	set.char.seed("9b87264c-e543-44ed-b4df-7e70cfbaf59c")

	x1 <- matrix(rnorm(1e2*5,mean=1),ncol=5)
	x1[1,1] <- NA
	x2 <- x1[2:nrow(x1),,drop=FALSE]

	max_ord <- 2L
	expect_error(rs1 <- cent_comoments(x1,max_ord,used_df=1L,na_omit=TRUE),NA)
	expect_error(rs2 <- cent_comoments(x2,max_ord,used_df=1L,na_omit=TRUE),NA)
	expect_equal(rs1,rs2)
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
	expect_error(rs3alt <- join_cent_cosums(rs1,rs2),NA)
	expect_equal(rs3,rs3alt,tolerance=1e-7)

	expect_error(rs1alt <- unjoin_cent_cosums(rs3,rs2),NA)
	expect_error(rs2alt <- unjoin_cent_cosums(rs3,rs1),NA)
	expect_equal(rs1,rs1alt,tolerance=1e-7)
	expect_equal(rs2,rs2alt,tolerance=1e-7)
})#UNFOLD
test_that("commutativity of join cosums",{#FOLDUP
	max_ord <- 2L

	set.char.seed("eb9064dc-9463-4a2a-b824-ec33defed3b6")
	x1 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	x2 <- matrix(rnorm(1e3*5,mean=1),ncol=5)
	expect_error(rs1 <- cent_cosums(x1,max_ord),NA)
	expect_error(rs2 <- cent_cosums(x2,max_ord),NA)

	expect_error(rs3_a <- join_cent_cosums(rs1,rs2),NA)
	expect_error(rs3_b <- join_cent_cosums(rs2,rs1),NA)
	expect_equal(rs3_a,rs3_b,tolerance=1e-7)

})#UNFOLD

# 2FIX: check the effects of NA

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
