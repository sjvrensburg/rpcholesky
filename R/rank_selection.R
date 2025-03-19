#' Find the minimum rank needed to achieve a target error threshold
#'
#' @param A Matrix to approximate
#' @param error_threshold Maximum relative error to tolerate
#' @param error_type Type of error to use ("frobenius", "spectral", or "trace")
#' @param max_rank Maximum rank to try
#' @param min_rank Minimum rank to consider
#' @param step_size Increment in rank for each iteration
#' @param min_variance_explained Minimum proportion of variance that must be explained (NA to disable)
#' @param binary_search Whether to use binary search to speed up the process
#' @param verbose Whether to print progress
#'
#' @return A list with optimal rank and error metrics
#' @export
find_optimal_rank <- function(A, error_threshold = 0.01,
                              error_type = "frobenius",
                              max_rank = min(nrow(A), 100),
                              min_rank = 1,
                              step_size = 5,
                              min_variance_explained = 0.8,
                              binary_search = TRUE,
                              verbose = FALSE) {

  if (!is.matrix(A)) A <- as.matrix(A)
  n <- nrow(A)

  # Optionally compute eigendecomposition for variance check
  if (!is.na(min_variance_explained)) {
    if (verbose) cat("Computing minimum rank for variance explained...\n")
    eig_info <- estimate_rank_eigenvalues(A, method = "threshold",
                                          threshold = min_variance_explained,
                                          max_rank = max_rank, plot = FALSE)
    min_variance_rank <- eig_info$optimal_rank
    if (verbose) cat("Minimum rank for", min_variance_explained*100, "% variance:", min_variance_rank, "\n")
  } else {
    min_variance_rank <- 0
  }

  # Error evaluation function
  evaluate_error <- function(k) {
    if (verbose) cat("Testing rank", k, "...\n")
    result <- rpcholesky(A, k, b = "auto", verbose = FALSE)
    errors <- rpcholesky_error(A, result)

    if (error_type == "frobenius") {
      return(list(error = errors$relative_frobenius_error, result = result, errors = errors))
    } else if (error_type == "spectral") {
      return(list(error = errors$relative_spectral_error, result = result, errors = errors))
    } else if (error_type == "trace") {
      return(list(error = errors$relative_trace_error, result = result, errors = errors))
    } else {
      stop("Unknown error type")
    }
  }

  if (binary_search) {
    # Binary search for optimal rank
    lower <- max(min_rank, 1)
    upper <- max_rank
    best_below_threshold <- NULL

    # Keep track of best result so far
    best_error <- Inf
    best_result <- NULL

    while (lower <= upper) {
      mid <- floor((lower + upper) / 2)
      eval_result <- evaluate_error(mid)

      # Track best result so far
      if (eval_result$error < best_error) {
        best_error <- eval_result$error
        best_result <- list(rank = mid,
                            error = eval_result$error,
                            result = eval_result$result,
                            errors = eval_result$errors)
      }

      if (eval_result$error <= error_threshold) {
        # Current rank is good, try a lower rank
        best_below_threshold <- list(rank = mid,
                                     error = eval_result$error,
                                     result = eval_result$result,
                                     errors = eval_result$errors)
        upper <- mid - 1
      } else {
        # Current rank is not good enough, try a higher rank
        lower <- mid + 1
      }
    }

    if (!is.null(best_below_threshold)) {
      # If we found a rank below threshold, fine-tune it
      found_rank <- best_below_threshold$rank

      # Additional tuning by checking ranks below the found rank
      if (found_rank > min_rank + 1) {
        for (k in (found_rank - 1):max(min_rank, found_rank - 5)) {
          eval_result <- evaluate_error(k)
          if (eval_result$error <= error_threshold) {
            best_below_threshold <- list(rank = k,
                                         error = eval_result$error,
                                         result = eval_result$result,
                                         errors = eval_result$errors)
          } else {
            break
          }
        }
      }

      # Apply minimum variance constraint if specified
      if (!is.na(min_variance_explained) && best_below_threshold$rank < min_variance_rank) {
        if (verbose) cat("Increasing rank to meet minimum variance explained:", min_variance_rank, "\n")
        eval_result <- evaluate_error(min_variance_rank)
        best_below_threshold <- list(rank = min_variance_rank,
                                     error = eval_result$error,
                                     result = eval_result$result,
                                     errors = eval_result$errors,
                                     note = "Rank increased to meet minimum variance explained")
      }

      best_below_threshold$status <- "found"
      return(best_below_threshold)
    } else {
      # No rank meets the threshold
      if (is.null(best_result)) {
        eval_result <- evaluate_error(max_rank)
        best_result <- list(rank = max_rank,
                            error = eval_result$error,
                            result = eval_result$result,
                            errors = eval_result$errors)
      }

      best_result$status <- "max_rank_reached"

      # If min_variance_rank is higher than max_rank, note this limitation
      if (!is.na(min_variance_explained) && min_variance_rank > max_rank) {
        best_result$note <- paste("Minimum variance explained would require rank", min_variance_rank)
      }

      return(best_result)
    }
  } else {
    # Linear search with step size
    min_error_rank <- NULL
    min_error_value <- Inf
    found_threshold <- FALSE

    for (k in seq(min_rank, max_rank, by = step_size)) {
      eval_result <- evaluate_error(k)

      if (eval_result$error < min_error_value) {
        min_error_value <- eval_result$error
        min_error_rank <- list(rank = k,
                               error = eval_result$error,
                               result = eval_result$result,
                               errors = eval_result$errors)
      }

      if (!found_threshold && eval_result$error <= error_threshold) {
        found_threshold <- TRUE

        # Found a rank that meets threshold, fine-tune it
        best_rank <- k

        # Check smaller ranks with step size 1
        if (k > min_rank) {
          for (j in (k - 1):max(min_rank, k - step_size + 1)) {
            fine_eval <- evaluate_error(j)
            if (fine_eval$error <= error_threshold) {
              best_rank <- j
            } else {
              break
            }
          }
        }

        # Apply minimum variance constraint
        final_rank <- best_rank
        if (!is.na(min_variance_explained) && best_rank < min_variance_rank) {
          if (verbose) cat("Increasing rank to meet minimum variance explained:", min_variance_rank, "\n")
          final_rank <- min_variance_rank
          note <- "Rank increased to meet minimum variance explained"
        } else {
          note <- NULL
        }

        final_eval <- evaluate_error(final_rank)
        result <- list(rank = final_rank,
                       error = final_eval$error,
                       result = final_eval$result,
                       errors = final_eval$errors,
                       status = "found")
        if (!is.null(note)) result$note <- note

        return(result)
      }
    }

    # If we reach here, we couldn't find a rank meeting the error threshold
    if (is.null(min_error_rank)) {
      eval_result <- evaluate_error(max_rank)
      min_error_rank <- list(rank = max_rank,
                             error = eval_result$error,
                             result = eval_result$result,
                             errors = eval_result$errors)
    }

    # Apply minimum variance constraint if applicable
    final_rank <- min_error_rank$rank
    if (!is.na(min_variance_explained) && final_rank < min_variance_rank) {
      if (verbose) cat("Increasing rank to meet minimum variance explained:", min_variance_rank, "\n")
      eval_result <- evaluate_error(min_variance_rank)
      result <- list(rank = min_variance_rank,
                     error = eval_result$error,
                     result = eval_result$result,
                     errors = eval_result$errors,
                     status = "min_variance",
                     note = "Rank increased to meet minimum variance explained")
      return(result)
    }

    min_error_rank$status <- "max_rank_reached"
    return(min_error_rank)
  }
}

#' Estimate optimal rank based on eigenvalue decay
#'
#' @param A Matrix to analyze
#' @param method Method to use: "elbow" or "threshold"
#' @param threshold Fraction of total variance to capture (if method="threshold")
#'                  For "elbow" method, this sets a minimum variance floor
#' @param max_rank Maximum number of eigenvalues to compute
#' @param plot Whether to create a visualization of the eigenvalue decay
#'
#' @return Estimated optimal rank
#' @export
estimate_rank_eigenvalues <- function(A, method = "elbow", threshold = 0.95,
                                      max_rank = min(nrow(A), 100),
                                      plot = TRUE) {
  if (!is.matrix(A)) A <- as.matrix(A)
  n <- nrow(A)

  # Compute the eigendecomposition (or SVD for efficiency with large matrices)
  if (n <= 1000) {
    eig <- eigen(A, symmetric = TRUE)
    eigenvalues <- eig$values[1:min(max_rank, n)]
  } else {
    # For large matrices, use partial SVD for efficiency
    svd_result <- svd(A, nu = 0, nv = 0)
    eigenvalues <- svd_result$d[1:min(max_rank, length(svd_result$d))]
  }

  # Ensure eigenvalues are non-negative (for numerical stability)
  eigenvalues <- pmax(eigenvalues, 0)

  # Normalize eigenvalues
  total_variance <- sum(eigenvalues)
  norm_eigenvalues <- eigenvalues / total_variance
  cum_variance <- cumsum(norm_eigenvalues)

  # Determine optimal rank
  if (method == "threshold") {
    optimal_rank <- which(cum_variance >= threshold)[1]
    if (is.na(optimal_rank)) {
      optimal_rank <- length(eigenvalues)
      warning("Could not reach target threshold with max_rank. Using all available eigenvalues.")
    }
  } else if (method == "elbow") {
    # IMPROVED ELBOW METHOD: Uses maximum distance from line approach
    n_eig <- length(eigenvalues)

    # Create points from (index, cumulative variance) coordinates
    points <- cbind(1:n_eig, cum_variance)

    # Define the line from first to last point
    start_point <- points[1,]
    end_point <- points[n_eig,]

    # Compute the distance from each point to the line
    # Line equation: ax + by + c = 0
    a <- end_point[2] - start_point[2]
    b <- start_point[1] - end_point[1]
    c <- end_point[1] * start_point[2] - start_point[1] * end_point[2]

    # Normalize for proper distance calculation
    norm_factor <- sqrt(a^2 + b^2)
    a <- a / norm_factor
    b <- b / norm_factor
    c <- c / norm_factor

    # Calculate distances
    distances <- abs(a * points[,1] + b * points[,2] + c)

    # Find index with maximum distance (the elbow)
    elbow_idx <- which.max(distances)

    # Apply a minimum variance constraint to prevent too low ranks
    min_variance <- if (is.null(threshold)) 0.7 else threshold

    if (cum_variance[elbow_idx] < min_variance) {
      # Find the smallest rank that explains at least min_variance
      threshold_idx <- which(cum_variance >= min_variance)[1]
      if (!is.na(threshold_idx)) {
        optimal_rank <- threshold_idx
        method_used <- "elbow with minimum variance constraint"
      } else {
        optimal_rank <- n_eig
        method_used <- "maximum available rank (variance constraint not met)"
      }
    } else {
      optimal_rank <- elbow_idx
      method_used <- "elbow"
    }
  } else {
    stop("Unknown method. Use 'elbow' or 'threshold'.")
  }

  # Plotting the eigenvalue decay
  if (plot) {
    par(mfrow = c(1, 2))
    # Plot eigenvalue decay
    plot(1:length(eigenvalues), eigenvalues, type = "b",
         xlab = "Index", ylab = "Eigenvalue",
         main = "Eigenvalue Decay")
    abline(v = optimal_rank, col = "red", lty = 2)

    # Plot cumulative explained variance
    plot(1:length(cum_variance), cum_variance, type = "b",
         xlab = "Number of components", ylab = "Cumulative Proportion of Variance",
         main = "Cumulative Explained Variance")
    abline(v = optimal_rank, col = "red", lty = 2)
    abline(h = threshold, col = "blue", lty = 3)
    par(mfrow = c(1, 1))
  }

  return(list(
    optimal_rank = optimal_rank,
    eigenvalues = eigenvalues,
    cumulative_variance = cum_variance,
    method = ifelse(exists("method_used"), method_used, method),
    threshold = threshold,
    explained_variance = cum_variance[optimal_rank]
  ))
}

#' Profile approximation error across multiple ranks
#'
#' @param A Matrix to approximate
#' @param ranks Vector of ranks to test
#' @param error_types Types of errors to calculate
#' @param verbose Whether to show progress
#' @param plot Whether to create a visualization
#'
#' @return Data frame with errors at each rank
#' @export
profile_rank_errors <- function(A,
                                ranks = seq(5, min(nrow(A), 100), by = 5),
                                error_types = c("frobenius", "trace", "max"),
                                verbose = TRUE,
                                plot = TRUE) {
  if (!is.matrix(A)) A <- as.matrix(A)

  results <- data.frame(rank = ranks)
  for (err_type in error_types) {
    results[[paste0("rel_", err_type, "_error")]] <- NA
  }

  # Calculate errors for each rank
  for (i in seq_along(ranks)) {
    k <- ranks[i]
    if (verbose) cat("Testing rank", k, "...\n")

    # Run rpcholesky
    rpchol_res <- rpcholesky(A, k, b = "auto", verbose = FALSE)

    # Calculate errors
    errors <- rpcholesky_error(A, rpchol_res)

    # Store the relevant errors
    if ("frobenius" %in% error_types) {
      results$rel_frobenius_error[i] <- errors$relative_frobenius_error
    }
    if ("spectral" %in% error_types) {
      results$rel_spectral_error[i] <- errors$relative_spectral_error
    }
    if ("trace" %in% error_types) {
      results$rel_trace_error[i] <- errors$relative_trace_error
    }
    if ("max" %in% error_types) {
      results$rel_max_error[i] <- errors$relative_max_error
    }
  }

  # Create visualization
  if (plot) {
    error_cols <- grep("^rel_", names(results), value = TRUE)

    # Reshape data for ggplot
    plot_data <- reshape2::melt(results,
                                id.vars = "rank",
                                measure.vars = error_cols,
                                variable.name = "error_type",
                                value.name = "error")

    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = rank, y = error, color = error_type)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::scale_y_log10() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Approximation Error vs Rank",
        x = "Rank",
        y = "Relative Error (log scale)",
        color = "Error Type"
      )

    print(p)
  }

  return(results)
}

#' Compare rpcholesky with SVD for low-rank approximation
#'
#' @param A Matrix to approximate
#' @param ranks Vector of ranks to test
#' @param metrics Error metrics to compute
#' @param time_benchmark Whether to benchmark computation time
#' @param plot Whether to create visualizations
#'
#' @return List with comparative results
#' @export
compare_with_svd <- function(A,
                             ranks = seq(5, min(nrow(A), 50), by = 5),
                             metrics = c("error", "time"),
                             time_benchmark = TRUE,
                             plot = TRUE) {
  if (!is.matrix(A)) A <- as.matrix(A)

  comparison <- data.frame(
    rank = ranks,
    rpchol_error = NA,
    svd_error = NA,
    error_ratio = NA
  )

  if (time_benchmark) {
    comparison$rpchol_time <- NA
    comparison$svd_time <- NA
    comparison$speedup <- NA
  }

  # Run comparisons for each rank
  for (i in seq_along(ranks)) {
    k <- ranks[i]

    # Timing and error for rpcholesky
    if (time_benchmark) {
      rpchol_time <- system.time({
        rpchol_result <- rpcholesky(A, k, b = "auto", verbose = FALSE)
        G <- rpchol_result$G
        A_rpchol <- t(G) %*% G
      })
      comparison$rpchol_time[i] <- rpchol_time[3]  # elapsed time
    } else {
      rpchol_result <- rpcholesky(A, k, b = "auto", verbose = FALSE)
      G <- rpchol_result$G
      A_rpchol <- t(G) %*% G
    }

    # Timing and error for SVD
    if (time_benchmark) {
      svd_time <- system.time({
        svd_result <- svd(A, nu = k, nv = k)
        A_svd <- svd_result$u %*% diag(svd_result$d[1:k]) %*% t(svd_result$v)
      })
      comparison$svd_time[i] <- svd_time[3]  # elapsed time
      comparison$speedup[i] <- svd_time[3] / rpchol_time[3]
    } else {
      svd_result <- svd(A, nu = k, nv = k)
      A_svd <- svd_result$u %*% diag(svd_result$d[1:k]) %*% t(svd_result$v)
    }

    # Calculate errors
    frob_norm_orig <- sqrt(sum(A^2))
    rpchol_error <- sqrt(sum((A - A_rpchol)^2)) / frob_norm_orig
    svd_error <- sqrt(sum((A - A_svd)^2)) / frob_norm_orig

    comparison$rpchol_error[i] <- rpchol_error
    comparison$svd_error[i] <- svd_error
    comparison$error_ratio[i] <- rpchol_error / svd_error
  }

  # Create visualizations
  if (plot) {
    if ("error" %in% metrics) {
      # Plot error comparison
      error_data <- reshape2::melt(
        comparison,
        id.vars = "rank",
        measure.vars = c("rpchol_error", "svd_error"),
        variable.name = "method",
        value.name = "error"
      )

      p1 <- ggplot2::ggplot(error_data, ggplot2::aes(x = rank, y = error, color = method)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_y_log10() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Approximation Error Comparison",
          x = "Rank",
          y = "Relative Frobenius Error (log scale)",
          color = "Method"
        )

      print(p1)

      # Plot error ratio
      p2 <- ggplot2::ggplot(comparison, ggplot2::aes(x = rank, y = error_ratio)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Error Ratio (rpcholesky / SVD)",
          x = "Rank",
          y = "Error Ratio"
        )

      print(p2)
    }

    if (time_benchmark && "time" %in% metrics) {
      # Plot time comparison
      time_data <- reshape2::melt(
        comparison,
        id.vars = "rank",
        measure.vars = c("rpchol_time", "svd_time"),
        variable.name = "method",
        value.name = "time"
      )

      p3 <- ggplot2::ggplot(time_data, ggplot2::aes(x = rank, y = time, color = method)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Computation Time Comparison",
          x = "Rank",
          y = "Time (seconds)",
          color = "Method"
        )

      print(p3)

      # Plot speedup
      p4 <- ggplot2::ggplot(comparison, ggplot2::aes(x = rank, y = speedup)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Computational Speedup (SVD time / rpcholesky time)",
          x = "Rank",
          y = "Speedup Factor"
        )

      print(p4)
    }
  }

  return(comparison)
}

#' Analyze the pivot selection behavior of rpcholesky
#'
#' @param A Matrix being approximated
#' @param result Result from rpcholesky
#' @param compute_eigendecomp Whether to compute eigendecomposition
#' @param plot Whether to create visualizations
#'
#' @return Analysis of pivot selection
#' @export
analyze_pivots <- function(A, result, compute_eigendecomp = TRUE, plot = TRUE) {
  if (!is.matrix(A)) A <- as.matrix(A)

  if (is.list(result) && "idx" %in% names(result)) {
    pivots <- result$idx
  } else {
    stop("Result must be the output of rpcholesky function")
  }

  n <- nrow(A)

  # Vector of selected positions
  selected <- rep(FALSE, n)
  selected[pivots + 1] <- TRUE  # Add 1 for R indexing

  # Basic pivot statistics
  pivot_stats <- list(
    total_pivots = length(pivots),
    pivot_fraction = length(pivots) / n,
    min_index = min(pivots),
    max_index = max(pivots),
    median_index = median(pivots)
  )

  # Compute diagonal elements
  diag_values <- diag(A)
  pivot_diag_values <- diag_values[pivots + 1]

  # Compare diagonal statistics
  diag_stats <- list(
    mean_diag_value = mean(diag_values),
    mean_pivot_diag = mean(pivot_diag_values),
    median_diag_value = median(diag_values),
    median_pivot_diag = median(pivot_diag_values),
    min_diag_value = min(diag_values),
    min_pivot_diag = min(pivot_diag_values),
    max_diag_value = max(diag_values),
    max_pivot_diag = max(pivot_diag_values)
  )

  # Eigenvalue analysis if requested
  eigen_stats <- list()
  if (compute_eigendecomp) {
    # Compute eigendecomposition
    eig <- tryCatch({
      eigen(A, symmetric = TRUE)
    }, error = function(e) {
      message("Could not compute full eigendecomposition. Using partial SVD.")
      # Use SVD which is more stable for large matrices
      svd_result <- svd(A)
      list(values = svd_result$d,
           vectors = svd_result$u)
    })

    # Calculate the "importance" of each row/column in the eigenvector basis
    importance <- rowSums(eig$vectors[, 1:length(pivots)]^2)

    # Compare the importance of selected pivots vs non-selected
    selected_importance <- importance[pivots + 1]
    non_selected_importance <- importance[!(seq_len(n) %in% (pivots + 1))]

    eigen_stats <- list(
      mean_selected_importance = mean(selected_importance),
      mean_non_selected_importance = mean(non_selected_importance),
      importance_ratio = mean(selected_importance) / mean(non_selected_importance),
      correlation_with_diag = cor(importance, diag_values)
    )

    # Visualization
    if (plot) {
      # Create data frame for plotting
      plot_data <- data.frame(
        index = 1:n,
        diagonal = diag_values,
        importance = importance,
        selected = selected
      )

      # Diagonal values vs selection
      p1 <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = index, y = diagonal, color = selected)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Diagonal Values and Pivot Selection",
          x = "Index",
          y = "Diagonal Value",
          color = "Selected as Pivot"
        )

      # Importance vs selection
      p2 <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = index, y = importance, color = selected)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Eigenvector Importance and Pivot Selection",
          x = "Index",
          y = "Eigenvector Importance",
          color = "Selected as Pivot"
        )

      # Relationship between diagonal and importance
      p3 <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = diagonal, y = importance, color = selected)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Diagonal vs Eigenvector Importance",
          x = "Diagonal Value",
          y = "Eigenvector Importance",
          color = "Selected as Pivot"
        )

      gridExtra::grid.arrange(p1, p2, p3, ncol = 2)
    }
  }

  # Combine all stats
  return(list(
    pivot_statistics = pivot_stats,
    diagonal_statistics = diag_stats,
    eigenvalue_statistics = eigen_stats,
    pivots = pivots,
    selected = selected
  ))
}
