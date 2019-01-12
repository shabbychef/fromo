#' @note
#' Note that when weights are given, they are treated as replication weights.
#' This can have subtle effects on computations which require minimum
#' degrees of freedom, since the sum of weights will be compared to
#' that minimum, not the number of data points. Weight values
#' (much) less than 1 can cause computations to return \code{NA}
#' somewhat unexpectedly due to this condition, while values greater
#' than one might cause the computation to spuriously return a value
#' with little precision.
