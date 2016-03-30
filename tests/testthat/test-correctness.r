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
THOROUGHNESS <- getOption('test.thoroughness',1.0)
#UNFOLD

context("basic correctness")#FOLDUP
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
test_that("running ops are correct",{#FOLDUP
	set.char.seed("967d2149-fbff-4d82-b227-ca3e1034bddb")
	for (xlen in c(20,100)) {
		x <- rnorm(xlen)
		for (window in c(5,50,Inf)) {
			for (recoper in c(10,1000)) {
				for (na_rm in c(FALSE,TRUE)) {
					dumb_count <- sapply(seq_along(x),function(iii) { sum(sign(abs(x[max(1,iii-window+1):iii])+1),na.rm=na_rm) },simplify=TRUE)
					dumb_mean <- sapply(seq_along(x),function(iii) { mean(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)
					dumb_sd <- sapply(seq_along(x),function(iii) { sd(x[max(1,iii-window+1):iii],na.rm=na_rm) },simplify=TRUE)

					fastv <- running_sd3(x,window=window,recoper=recoper,na_rm=na_rm)
					dumbv <- cbind(dumb_sd,dumb_mean,dumb_count)
					expect_equal(max(abs(dumbv[2:xlen,] - fastv[2:xlen,])),0,tolerance=1e-12)

					fastv <- running_centered(x,window=window,recoper=recoper,na_rm=na_rm)
					# the dumb value:
					dumbv <- x - dumb_mean;
					expect_equal(max(abs(dumbv - fastv)),0,tolerance=1e-12)

					fastv <- running_scaled(x,window=window,recoper=recoper,na_rm=na_rm)
					# the dumb value:
					dumbv <- x / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_zscored(x,window=window,recoper=recoper,na_rm=na_rm)
					# the dumb value:
					dumbv <- (x - dumb_mean) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_sharpe(x,window=window,recoper=recoper,na_rm=na_rm)
					# the dumb value:
					dumbv <- dumb_mean / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)

					fastv <- running_tstat(x,window=window,recoper=recoper,na_rm=na_rm)
					# the dumb value:
					dumbv <- (dumb_mean * sqrt(dumb_count)) / dumb_sd
					expect_equal(max(abs(dumbv[2:length(x)] - fastv[2:length(x)])),0,tolerance=1e-12)
				}
			}
		}
	}

	# sentinel
	expect_true(TRUE)
})#UNFOLD
# 2FIX: check the effects of NA
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
