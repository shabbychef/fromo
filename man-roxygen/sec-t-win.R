#' @section Time Windowing :
#'
#' This function supports time (or other counter) based running computation. 
#' Here the input are the data \eqn{x_i}, and optional weights vectors, \eqn{w_i}, defaulting to 1,
#' and a vector of time indices, \eqn{t_i} of the same length as \eqn{x}. The
#' times must be non-decreasing:
#' \deqn{t_1 \le t_2 \le \ldots}{t_1 <= t_2 <= ...}
#' It is assumed that \eqn{t_0 = -\infty}.
#' The window, \eqn{W} is now a time-based window. 
#' An optional set of \emph{lookback times} are also given, \eqn{b_j}, which
#' may have different length than the \eqn{x} and \eqn{w}. 
#' The output will correspond to the lookback times, and should be the same
#' length. The \eqn{j}th output is computed over indices \eqn{i} such that
#' \deqn{b_j - W < t_i \le b_j.}
#'
#' For comparison functions (like Z-score, rescaling, centering), which compare
#' values of \eqn{x_i} to local moments, the lookbacks may not be given, but
#' a lookahead \eqn{L} is admitted. In this case, the \eqn{j}th output is computed over
#' indices \eqn{i} such that
#' \deqn{t_j - W + L < t_i \le t_j + L.}
#'
#' If the times are not given, \sQuote{deltas} may be given instead. If
#' \eqn{\delta_i} are the deltas, then we compute the times as
#' \deqn{t_i = \sum_{1 \le j \le i} \delta_j.}
#' The deltas must be the same length as \eqn{x}.  
#' If times and deltas are not given, but weights are given and the \sQuote{weights as deltas}
#' flag is set true, then the weights are used as the deltas.
#'
#' Some times it makes sense to have the computational window be the space
#' between lookback times. That is, the \eqn{j}th output is to be
#' computed over indices \eqn{i} such that
#' \deqn{b_{j-1} - W < t_i \le b_j.}
#' This can be achieved by setting the \sQuote{variable window} flag true
#' and setting the window to null. This will not make much sense if
#' the lookback times are equal to the times, since each moment computation
#' is over a set of a single index, and most moments are underdefined.
#'
