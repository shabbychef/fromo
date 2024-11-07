#' @note
#' As this code may add and remove observations, numerical imprecision
#' may result in negative estimates of squared quantities, like the
#' second or fourth moments.  By default we check for this condition
#' in running computations. It may also be mitigated somewhat by setting 
#' a smaller \code{restart_period}. Post an issue if you experience this bug.
