#' Calculate comprehensive error metrics for low-rank approximation
#'
#' @param A Original positive semidefinite matrix
#' @param approx Approximated matrix or approximation components
#' @param factor Whether approx contains the factors (G) or the complete approximation (G^T * G)
#'
#' @return A named list of various error metrics
#' @export
rpcholesky_error <- function(A, approx, factor = TRUE) {
  if (!is.matrix(A)) A <- as.matrix(A)

  # If approx is a factor (G), compute the approximation
  if (factor) {
    if (is.list(approx) && "G" %in% names(approx)) {
      G <- approx$G
      A_approx <- t(G) %*% G
    } else {
      G <- approx
      A_approx <- t(G) %*% G
    }
  } else {
    A_approx <- approx
  }

  # Calculate various error metrics
  diff <- A - A_approx

  # Frobenius norm error
  frob_norm_orig <- sqrt(sum(A^2))
  frob_error <- sqrt(sum(diff^2))
  rel_frob_error <- frob_error / frob_norm_orig

  # Spectral norm (using largest singular value)
  spec_error <- NA
  tryCatch({
    spec_error <- max(svd(diff, nu = 1, nv = 1)$d)
    rel_spec_error <- spec_error / max(svd(A, nu = 1, nv = 1)$d)
  }, error = function(e) {
    message("Could not compute spectral error for large matrices")
  })

  # Trace error
  trace_orig <- sum(diag(A))
  trace_approx <- sum(diag(A_approx))
  trace_error <- abs(trace_orig - trace_approx)
  rel_trace_error <- trace_error / trace_orig

  # Maximum absolute entry error
  max_error <- max(abs(diff))
  rel_max_error <- max_error / max(abs(A))

  # Return all metrics
  list(
    frobenius_error = frob_error,
    relative_frobenius_error = rel_frob_error,
    spectral_error = spec_error,
    relative_spectral_error = rel_spec_error,
    trace_error = trace_error,
    relative_trace_error = rel_trace_error,
    max_entry_error = max_error,
    relative_max_error = rel_max_error
  )
}

#' Calculate kernel alignment between original and approximated matrices
#'
#' @param K Original kernel matrix
#' @param K_approx Approximated kernel matrix or approximation components
#' @param factor Whether K_approx contains the factors (G) or the complete approximation
#'
#' @return A scalar representing the alignment between the matrices (1 is perfect)
#' @export
kernel_alignment <- function(K, K_approx, factor = TRUE) {
  if (!is.matrix(K)) K <- as.matrix(K)

  # If K_approx is a factor (G), compute the approximation
  if (factor) {
    if (is.list(K_approx) && "G" %in% names(K_approx)) {
      G <- K_approx$G
      K_approx <- t(G) %*% G
    } else {
      G <- K_approx
      K_approx <- t(G) %*% G
    }
  }

  # Compute Frobenius inner product <K1, K2>_F = trace(K1^T K2)
  inner_prod <- sum(K * K_approx)
  norm_K <- sqrt(sum(K^2))
  norm_K_approx <- sqrt(sum(K_approx^2))

  # Return the normalized alignment
  inner_prod / (norm_K * norm_K_approx)
}
