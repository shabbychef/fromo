#' @param check_negative_moments  a boolean flag. Normal computation of running
#' moments can result in negative estimates of even order moments due to loss of
#' numerical precision. With this flag active, the computation checks for negative
#' even order moments and restarts the computation when one is detected. This
#' should eliminate the possibility of negative even order moments. The
#' downside is the speed hit of checking on every output step. Note also the
#' code checks for negative moments of every even order tracked, even if they
#' are not output; that is if the kurtosis, say, is being computed, and a
#' negative variance is detected, then the computation is restarted.
#' Defaults to \code{TRUE} to avoid negative even moments. Set to \code{FALSE}
#' only if you know what you are doing.
