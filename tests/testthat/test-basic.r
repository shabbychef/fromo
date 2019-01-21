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

# code runs at all
context("code runs: sd, skew, kurt")
test_that("sd, skew, kurt run without error",{#FOLDUP
	set.char.seed("569dd47d-f9e5-40e4-b2ac-e5dbb4771a53")
	x <- rnorm(100)
	y <- as.integer(x)
	z <- as.logical(y)
	q <- c('a','b','c')
	wts <- 100 * runif(length(x))

	for (na_rm in c(FALSE,TRUE)) {
		for (fnc in list(function(x) { x },as.integer,function(x) { as.logical(x > 0) })) {
			cls_x <- fnc(x)
			expect_error(sd3(cls_x,na_rm=na_rm),NA)
			expect_error(skew4(cls_x,na_rm=na_rm),NA)
			expect_error(kurt5(cls_x,na_rm=na_rm),NA)
			expect_error(cent_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(std_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(cent_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)
			expect_error(std_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm),NA)

			# weights!?
			for (wfnc in list(function(x) { x },as.integer,function(x) { as.logical(x > 10) })) {
				cls_wts <- wfnc(wts)
				expect_error(sd3(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(skew4(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(kurt5(cls_x,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(kurt5(cls_x,na_rm=na_rm,wts=cls_wts,check_wts=TRUE),NA)
				
				expect_error(cent_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(std_moments(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(cent_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
				expect_error(std_cumulants(cls_x,max_order=5L,used_df=1L,na_rm=na_rm,wts=cls_wts),NA)
			}
		}
	}
	# check some things about weights;
	expect_error(sd3(x,wts=q))
	expect_error(sd3(x,wts=rep(-1,length(x)),check_wts=TRUE))

	for (checkme in list(sd3,skew4,kurt5,cent_moments,std_moments,cent_cumulants,std_cumulants)) {
		expect_error(checkme(q))
		expect_error(checkme(x,wts=q))
		expect_error(checkme(x,wts=rep(-1,length(x)),check_wts=TRUE))
	}
})#UNFOLD
context("code runs: cosum and comoment")
test_that("cosum and comoment run without error",{#FOLDUP
	set.char.seed("ff34b509-a113-41c9-8517-aa72792c42f7")
	x <- matrix(rnorm(30*4),ncol=4)
	y <- matrix(as.integer(x),ncol=ncol(x))
	z <- matrix(as.logical(y),ncol=ncol(x))
	q <- matrix(letters[1:24],ncol=4)

	for (na_omit in c(FALSE,TRUE)) {
		expect_error(cent_cosums(x,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(x,max_order=2L,used_df=1L,na_omit=na_omit),NA)

		expect_error(cent_cosums(y,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(y,max_order=2L,used_df=1L,na_omit=na_omit),NA)

		expect_error(cent_cosums(z,max_order=2L,na_omit=na_omit),NA)
		expect_error(cent_comoments(z,max_order=2L,used_df=1L,na_omit=na_omit),NA)
	}

	expect_error(cent_cosums(x,max_order=4L))
	expect_error(cent_cosums(q,max_order=2L))
	expect_error(cent_comoments(q,max_order=2L))
})#UNFOLD
# 2FIX: check the effects of NA
#context("order of parameters")# FOLDUP
# named parameters can be given in any order?

# UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
