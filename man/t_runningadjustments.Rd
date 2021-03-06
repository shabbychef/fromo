% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{t_running_centered}
\alias{t_running_centered}
\alias{t_running_scaled}
\alias{t_running_zscored}
\alias{t_running_sharpe}
\alias{t_running_tstat}
\title{Compare data to moments computed over a time sliding window.}
\usage{
t_running_centered(v, time = NULL, time_deltas = NULL, window = NULL,
  wts = NULL, na_rm = FALSE, min_df = 0L, used_df = 1, lookahead = 0,
  restart_period = 100L, variable_win = FALSE, wts_as_delta = TRUE,
  check_wts = FALSE, normalize_wts = TRUE)

t_running_scaled(v, time = NULL, time_deltas = NULL, window = NULL,
  wts = NULL, na_rm = FALSE, min_df = 0L, used_df = 1, lookahead = 0,
  restart_period = 100L, variable_win = FALSE, wts_as_delta = TRUE,
  check_wts = FALSE, normalize_wts = TRUE)

t_running_zscored(v, time = NULL, time_deltas = NULL, window = NULL,
  wts = NULL, na_rm = FALSE, min_df = 0L, used_df = 1, lookahead = 0,
  restart_period = 100L, variable_win = FALSE, wts_as_delta = TRUE,
  check_wts = FALSE, normalize_wts = TRUE)

t_running_sharpe(v, time = NULL, time_deltas = NULL, window = NULL,
  wts = NULL, lb_time = NULL, na_rm = FALSE, compute_se = FALSE,
  min_df = 0L, used_df = 1, restart_period = 100L, variable_win = FALSE,
  wts_as_delta = TRUE, check_wts = FALSE, normalize_wts = TRUE)

t_running_tstat(v, time = NULL, time_deltas = NULL, window = NULL,
  wts = NULL, lb_time = NULL, na_rm = FALSE, compute_se = FALSE,
  min_df = 0L, used_df = 1, restart_period = 100L, variable_win = FALSE,
  wts_as_delta = TRUE, check_wts = FALSE, normalize_wts = TRUE)
}
\arguments{
\item{v}{a vector of data.}

\item{time}{an optional vector of the timestamps of \code{v}. If given, must be
the same length as \code{v}. If not given, we try to infer it by summing the
\code{time_deltas}.}

\item{time_deltas}{an optional vector of the deltas of timestamps. If given, must be
the same length as \code{v}. If not given, and \code{wts} are given and \code{wts_as_delta} is true,
we take the \code{wts} as the time deltas.  The deltas must be positive. We sum them to arrive
at the times.}

\item{window}{the window size, in time units. if given as finite integer or double, passed through.
If \code{NULL}, \code{NA_integer_}, \code{NA_real_} or \code{Inf} are given, 
 and \code{variable_win} is true, then we infer the window from the lookback times: the
 first window is infinite, but the remaining is the deltas between lookback times.
 If \code{variable_win} is false, then these undefined values are equivalent to an
 infinite window.
 If negative, an error will be thrown.}

\item{wts}{an optional vector of weights. Weights are \sQuote{replication}
weights, meaning a value of 2 is shorthand for having two observations
with the corresponding \code{v} value. If \code{NULL}, corresponds to
equal unit weights, the default. Note that weights are typically only meaningfully defined
up to a multiplicative constant, meaning the units of weights are
immaterial, with the exception that methods which check for minimum df will,
in the weighted case, check against the sum of weights. For this reason,
weights less than 1 could cause \code{NA} to be returned unexpectedly due
to the minimum condition. When weights are \code{NA}, the same rules for checking \code{v}
are applied. That is, the observation will not contribute to the moment
if the weight is \code{NA} when \code{na_rm} is true. When there is no
checking, an \code{NA} value will cause the output to be \code{NA}.}

\item{na_rm}{whether to remove NA, false by default.}

\item{min_df}{the minimum df to return a value, otherwise \code{NaN} is returned.
This can be used to prevent \emph{e.g.} Z-scores from being computed on only 3
observations. Defaults to zero, meaning no restriction, which can result in 
infinite Z-scores during the burn-in period.}

\item{used_df}{the number of degrees of freedom consumed, used in the denominator
of the centered moments computation. These are subtracted from the number of
observations.}

\item{lookahead}{for some of the operations, the value is compared to 
mean and standard deviation possibly using 'future' or 'past' information
by means of a non-zero lookahead. Positive values mean data are taken from
the future. This is in time units, and so should be a real.}

\item{restart_period}{the recompute period. because subtraction of elements can cause
loss of precision, the computation of moments is restarted periodically based on 
this parameter. Larger values mean fewer restarts and faster, though less accurate
results.}

\item{variable_win}{if true, and the \code{window} is not a concrete number,
the computation window becomes the time between lookback times.}

\item{wts_as_delta}{if true and the \code{time} and \code{time_deltas} are not
given, but \code{wts} are given, we take \code{wts} as the \code{time_deltas}.}

\item{check_wts}{a boolean for whether the code shall check for negative
weights, and throw an error when they are found. Default false for speed.}

\item{normalize_wts}{a boolean for whether the weights should be
renormalized to have a mean value of 1. This mean is computed over elements
which contribute to the moments, so if \code{na_rm} is set, that means non-NA
elements of \code{wts} that correspond to non-NA elements of the data
vector.}

\item{lb_time}{a vector of the times from which lookback will be performed. The output should
be the same size as this vector. If not given, defaults to \code{time}.}

\item{compute_se}{for \code{running_sharpe}, return an extra column of the
standard error, as computed by Mertens' correction.}
}
\value{
a vector the same size as the input consisting of the adjusted version of the input.
When there are not sufficient (non-nan) elements for the computation, \code{NaN} are returned.
}
\description{
Computes moments over a sliding window, then adjusts the data accordingly, centering, or scaling,
or z-scoring, and so on.
}
\details{
Given the length \eqn{n} vector \eqn{x}, for
a given index \eqn{i}, define \eqn{x^{(i)}}{x^(i)}
as the elements of \eqn{x} defined by the sliding time window (see the section
on time windowing).
Then define \eqn{\mu_i}{mu_i}, \eqn{\sigma_i}{sigma_i}
and \eqn{n_i}{n_i} as, respectively, the sample mean, standard deviation and number of
non-NA elements in \eqn{x^{(i)}}{x^(i)}. 

We compute output vector \eqn{m} the same size as \eqn{x}. 
For the 'centered' version of \eqn{x}, we have \eqn{m_i = x_i - \mu_i}{m_i = x_i - mu_i}.
For the 'scaled' version of \eqn{x}, we have \eqn{m_i = x_i / \sigma_i}{m_i = x_i / sigma_i}.
For the 'z-scored' version of \eqn{x}, we have \eqn{m_i = (x_i - \mu_i) / \sigma_i}{m_i = (x_i - mu_i) / sigma_i}.
For the 't-scored' version of \eqn{x}, we have \eqn{m_i = \sqrt{n_i} \mu_i / \sigma_i}{m_i = sqrt(n_i) mu_i / sigma_i}.

We also allow a 'lookahead' for some of these operations.
If positive, the moments are computed using data from larger indices;
if negative, from smaller indices.
}
\note{
The moment computations provided by fromo are 
numerically robust, but will often \emph{not} provide the
same results as the 'standard' implementations,
due to differences in roundoff. We make every attempt to balance
speed and robustness. User assumes all risk from using
the fromo package.

Note that when weights are given, they are treated as replication weights.
This can have subtle effects on computations which require minimum
degrees of freedom, since the sum of weights will be compared to
that minimum, not the number of data points. Weight values
(much) less than 1 can cause computations to return \code{NA}
somewhat unexpectedly due to this condition, while values greater
than one might cause the computation to spuriously return a value
with little precision.
}
\section{Time Windowing }{


This function supports time (or other counter) based running computation. 
Here the input are the data \eqn{x_i}, and optional weights vectors, \eqn{w_i}, defaulting to 1,
and a vector of time indices, \eqn{t_i} of the same length as \eqn{x}. The
times must be non-decreasing:
\deqn{t_1 \le t_2 \le \ldots}{t_1 <= t_2 <= ...}
It is assumed that \eqn{t_0 = -\infty}.
The window, \eqn{W} is now a time-based window. 
An optional set of \emph{lookback times} are also given, \eqn{b_j}, which
may have different length than the \eqn{x} and \eqn{w}. 
The output will correspond to the lookback times, and should be the same
length. The \eqn{j}th output is computed over indices \eqn{i} such that
\deqn{b_j - W < t_i \le b_j.}

For comparison functions (like Z-score, rescaling, centering), which compare
values of \eqn{x_i} to local moments, the lookbacks may not be given, but
a lookahead \eqn{L} is admitted. In this case, the \eqn{j}th output is computed over
indices \eqn{i} such that
\deqn{t_j - W + L < t_i \le t_j + L.}

If the times are not given, \sQuote{deltas} may be given instead. If
\eqn{\delta_i} are the deltas, then we compute the times as
\deqn{t_i = \sum_{1 \le j \le i} \delta_j.}
The deltas must be the same length as \eqn{x}.  
If times and deltas are not given, but weights are given and the \sQuote{weights as deltas}
flag is set true, then the weights are used as the deltas.

Some times it makes sense to have the computational window be the space
between lookback times. That is, the \eqn{j}th output is to be
computed over indices \eqn{i} such that
\deqn{b_{j-1} - W < t_i \le b_j.}
This can be achieved by setting the \sQuote{variable window} flag true
and setting the window to null. This will not make much sense if
the lookback times are equal to the times, since each moment computation
is over a set of a single index, and most moments are underdefined.
}

\references{
Terriberry, T. "Computing Higher-Order Moments Online."
\url{http://people.xiph.org/~tterribe/notes/homs.html}

J. Bennett, et. al., "Numerically Stable, Single-Pass, 
Parallel Statistics Algorithms," Proceedings of IEEE
International Conference on Cluster Computing, 2009.
\url{https://www.semanticscholar.org/paper/Numerically-stable-single-pass-parallel-statistics-Bennett-Grout/a83ed72a5ba86622d5eb6395299b46d51c901265}

Cook, J. D. "Accurately computing running variance."
\url{http://www.johndcook.com/standard_deviation.html}

Cook, J. D. "Comparing three methods of computing 
standard deviation."
\url{http://www.johndcook.com/blog/2008/09/26/comparing-three-methods-of-computing-standard-deviation}
}
\seealso{
\code{\link{running_centered}}, \code{\link{scale}}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
