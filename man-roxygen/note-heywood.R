#' @note
#' As this code may add and remove observations, numerical imprecision
#' may result in negative estimates of squared quantities, like the
#' second or fourth moments.  We do not currently correct for this
#' issue, although it may be somewhat mitigated by setting a smaller
#' \code{restart_period}. In the future we will add a check for
#' this case. Post an issue if you experience this bug.
