#' @note
#' Note that when weights are given, they are treated as replication weights.
#' This can have subtle effects on computations which require minimum
#' degrees of freedom, since the sum of weights will be compared to
#' that minimum, not the number of data points. Weight values
#' less than 1 should be used with caution, then.
