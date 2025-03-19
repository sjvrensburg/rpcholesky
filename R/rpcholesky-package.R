## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib rpcholesky, .registration = TRUE
## usethis namespace: end
NULL

#' Accelerated Randomly Pivoted Cholesky Decomposition
#'
#' @param A A positive semidefinite matrix
#' @param k Target rank for the approximation
#' @param b Batch size (an integer or "auto")
#' @param stoptol Early stopping tolerance
#' @param verbose If TRUE, prints progress information
#'
#' @return A list containing low-rank factor G, selected pivot indices, and corresponding rows
#' @export
rpcholesky <- function(A, k, b = "auto", stoptol = 1e-10, verbose = FALSE) {
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  # Call the C++ implementation
  result <- accelerated_rpcholesky(A, k, b, stoptol, verbose)
  return(result)
}
