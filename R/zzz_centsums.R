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

# Created: 2016.03.30
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <steven@corecast.io>
# Comments: Steven E. Pav
# Copyright 2016-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

#' @title centsums Class.
#'
#' @description 
#'
#' An S4 class to store (centered) sums of data, and to support operations on 
#' the same.
#'
#' @details
#'
#' A \code{centsums} object contains a vector value of the data count,
#' the mean, and the \eqn{k}th centered sum, for \eqn{k} up to some
#' maximum order.
#'
#' @slot sums a numeric vector of the sums.
#' @slot order the maximum order.
#'
#' @return An object of class \code{centsums}.
#' @keywords moments
#'
#' @examples 
#' obj <- new("centsums",sums=c(1000,1.234,0.235),order=2)
#'
#' @template etc
#' @template ref-romo
#' @name centsums-class
#' @rdname centsums-class
#' @exportClass centsums
#' @export
setClass("centsums", 
				 representation(sums="numeric",order="numeric"),
				 prototype(sums=c(0.0,0.0,0.0),
									 order=2),
				 validity=function(object) {
					 # ... 
					 # http://www.cyclismo.org/tutorial/R/s4Classes.html
					 if ((!is.null(object@order)) && (length(object@sum) != (1+object@order))) { return("bad dimensionality or order given.") }
					 return(TRUE)
				 }
)
# constructor method documentation
#  
#' @param .Object a \code{centums} object, or proto-object.
#' @rdname centsums-class
#' @aliases initialize,centsums-class
setMethod('initialize',
					signature('centsums'),
					function(.Object,sums,order=NA_real_) {
						if (is.null(order)) {
							order <- length(sums) - 1
						}
					 	.Object@sums <- sums
					 	.Object@order <- order

						.Object
					})

#'
#' @param sums a numeric vector.
#' @param order the order, defaulting to \code{length(sums)+1}.
#' @name centsums
#' @rdname centsums-class
#' @export
centsums <- function(sums,order=NULL) {
	if (is.null(order)) {
		order <- length(sums) + 1
	}
	retv <- new("centsums", sums=sums, order=order)
	invisible(retv)
}

#' @title Coerce to a centsums object.
#'
#' @description 
#'
#' Convert data to a \code{centsums} object.
#'
#' @details
#'
#' Computes the raw sums on data, and stuffs the results into a 
#' \code{centsums} object.
#'
#' @usage
#'
#' as.centsums(x, order=3, na.rm=TRUE)
#'
#' @param x a numeric, array, or matrix.
#' @param na.rm whether to remove \code{NA}.
#' @inheritParams centsums
#' @return A centsums object.
#' @template etc
#' @examples 
#' set.seed(123)
#' x <- rnorm(1000)
#' cs <- as.centsums(x, order=5)
#' @rdname as.centsums
#' @export as.centsums
as.centsums <- function(x, order=3, na.rm=TRUE) {
	UseMethod("as.centsums", x)
}
#' @rdname as.centsums
#' @export
#' @method as.centsums default
#' @aliases as.centsums
as.centsums.default <- function(x, order=3, na.rm=TRUE) {
	sums <- cent_sums(x, max_order=order, na_rm=na.rm)
	invisible(centsums(sums,order=order))
}

#' @title Accessor methods.
#'
#' @description
#'
#' Access slot data from a \code{centsums} object.
#'
#' @param x a \code{centsums} object.
#' @param type the type of moment to compute.
#' @template etc
#' @name accessor
#' @rdname accessor-methods
#' @aliases sums
#' @exportMethod sums
setGeneric('sums', signature="x", function(x) standardGeneric('sums'))
#' @rdname accessor-methods
#' @aliases sums,centsums-method
setMethod('sums', 'centsums', function(x) x@sums )

#' @rdname accessor-methods
#' @aliases moments
#' @exportMethod moments
setGeneric('moments', function(x,type=c('central','raw','standardized')) standardGeneric('moments'))
#' @rdname accessor-methods
#' @aliases moments,centsums-method
setMethod('moments', signature(x='centsums'),
	function(x,type=c('central','raw','standardized')) {
		# add used_df
		type <- match.arg(type)
		c_sums <- x@sums
		cmoments <- c(c_sums[1],c_sums[2:length(c_sums)] / c_sums[1])

		retv <- switch(type,
			raw={
				retv <- cent2raw(cmoments)
			},
			central={ 
				retv <- cmoments[2:length(cmoments)]
				retv[1] <- 0
			},
			standardized={
				retv <- cmoments[2:length(cmoments)]
				retv[1] <- 0.0
				if (length(v) > 1) {
					if (length(v) > 2) {
						sigma2 <- retv[2]
						retv[3:length(retv)] <- retv[3:length(retv)] / (sigma2 ^ ((3:length(retv))/2.0))
					}
					retv[2] <- 1.0
				}
				retv
			})
			retv
	})

# do not export this.
.join2 <- function(x,y) {
	x@sums <- join_cent_sums(x@sums,y@sums)
	x
}

#' @title concatenate centsums objects.
#' @description 
#'
#' Concatenate centsums objects.
#'
#' @param ... \code{centsums} objects
#' @rdname centsums-concat
#' @method c centsums
#' @export
#' @usage \\method{c}{centsums}(...)
c.centsums <- function(...) { 
	Reduce(.join2,list(...))
	x
} 

# show#FOLDUP
# 2FIX: add documentation and export
#' @title Show a centsums object.
#'
#' @description 
#'
#' Displays the centsums object.
#'
#' @usage
#'
#' show(object)
#'
#' @param object a \code{centsums} object.
#' @examples 
#' set.seed(123)
#' x <- rnorm(1000)
#' obj <- as.centsums(x, order=5)
#' obj
#' @template etc
#' @name show
#' @rdname show-methods
#' @exportMethod show
#' @aliases show
NULL
#' @rdname show-methods
#' @aliases show,centsums-method
setMethod('show', signature('centsums'), 
					function(object) {
						cat('class:', class(object), '\n')
						cat(' moms:', moments(object,'central'), '\n')
					})
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
