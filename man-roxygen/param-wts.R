#' @param wts an optional vector of weights. Weights are \sQuote{replication}
#' weights, meaning a value of 2 is shorthand for having two observations
#' with the corresponding \code{v} value. If \code{NULL}, corresponds to
#' equal weights, the default. Note that weights are typically only meaningfully defined
#' up to a multiplicative constant, meaning the units of weights are
#' immaterial. When weights are \code{NA}, the same rules for checking \code{v}
#' are applied.
#' @param check_wts a boolean for whether the code shall check for negative
#' weights, and throw an error when they are found. Default false for speed.
