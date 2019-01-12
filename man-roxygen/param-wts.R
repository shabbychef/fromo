#' @param wts an optional vector of weights. Weights are \sQuote{replication}
#' weights, meaning a value of 2 is shorthand for having two observations
#' with the corresponding \code{v} value. If \code{NULL}, corresponds to
#' equal unit weights, the default. Note that weights are typically only meaningfully defined
#' up to a multiplicative constant, meaning the units of weights are
#' immaterial, with the exception that methods which check for minimum df will,
#' in the weighted case, check against the sum of weights. For this reason,
#' weights less than 1 could cause \code{NA} to be returned unexpectedly due
#' to the minimum condition. When weights are \code{NA}, the same rules for checking \code{v}
#' are applied. That is, the observation will not contribute to the moment
#' if the weight is \code{NA} when \code{na_rm} is true. When there is no
#' checking, an \code{NA} value will cause the output to be \code{NA}.
#' @param check_wts a boolean for whether the code shall check for negative
#' weights, and throw an error when they are found. Default false for speed.
#' @param normalize_wts a boolean for whether the weights should be
#' renormalized to have a mean value of 1. This mean is computed over elements
#' which contribute to the moments, so if \code{na_rm} is set, that means non-NA
#' elements of \code{wts} that correspond to non-NA elements of the data
#' vector.
