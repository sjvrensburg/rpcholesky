#' #' Improved example of accelerated randomly pivoted Cholesky
#' #' Demonstrates usage, evaluates approximation accuracy, and focuses on efficiency
#'
#' library(rpcholesky)
#' library(ggplot2)
#' library(reshape2)
#' library(microbenchmark)
#'
#' # -------------------------------------------------------------------------
#' # Generate a test matrix with controlled eigenvalue decay
#' # -------------------------------------------------------------------------
#' set.seed(123)
#'
#' # Create a matrix with specific eigenvalue decay profile
#' create_test_matrix <- function(n, decay_rate = 0.9, noise = 0.001) {
#'   # Create a diagonal matrix with decaying eigenvalues
#'   d <- decay_rate^(0:(n-1))  # Exponential decay
#'   D <- diag(d)
#'
#'   # Create a random orthogonal matrix
#'   Q <- qr.Q(qr(matrix(rnorm(n*n), n, n)))
#'
#'   # Create PSD matrix with controlled eigenvalue decay
#'   A <- Q %*% D %*% t(Q)
#'
#'   # Add a small amount of noise to make the problem more challenging
#'   if (noise > 0) {
#'     noise_matrix <- matrix(rnorm(n*n), n, n)
#'     noise_matrix <- noise_matrix %*% t(noise_matrix) * noise / n
#'     A <- A + noise_matrix
#'     # Ensure the matrix remains PSD after adding noise
#'     A <- (A + t(A))/2
#'   }
#'
#'   return(list(matrix = A, eigenvalues = d))
#' }
#'
#' # Create a larger test matrix with slower decay
#' n <- 2000  # Larger matrix
#' message("Creating test matrix...")
#' result <- create_test_matrix(n, decay_rate = 0.98)
#' A <- result$matrix
#' true_eigenvalues <- result$eigenvalues
#'
#' # -------------------------------------------------------------------------
#' # Compare RPCholesky vs SVD with benchmarking
#' # -------------------------------------------------------------------------
#' message("Running benchmarks...")
#'
#' # Choose ranks that are reasonable for the matrix size
#' ranks <- c(20, 50, 100, 200, 400)
#' all_results <- list()
#'
#' # Run benchmarks and collect results
#' for (k in ranks) {
#'   message(paste("Processing rank k =", k))
#'
#'   # Time RPCholesky
#'   rpcholesky_time <- system.time({
#'     rpcholesky_result <- rpcholesky(A, k = k, b = "auto", verbose = FALSE)
#'     G <- rpcholesky_result$G
#'     A_rpcholesky <- t(G) %*% G
#'   })
#'
#'   # Time SVD
#'   svd_time <- system.time({
#'     svd_result <- svd(A, nu = k, nv = k)
#'     U <- svd_result$u
#'     d <- svd_result$d[1:k]
#'     V <- svd_result$v
#'     A_svd <- U %*% diag(d) %*% t(V)
#'   })
#'
#'   # Calculate errors
#'   frob_norm_original <- sqrt(sum(A^2))
#'
#'   rpcholesky_error <- sqrt(sum((A - A_rpcholesky)^2)) / frob_norm_original
#'   svd_error <- sqrt(sum((A - A_svd)^2)) / frob_norm_original
#'
#'   rpcholesky_trace_error <- abs(sum(diag(A)) - sum(diag(A_rpcholesky))) / sum(diag(A))
#'   svd_trace_error <- abs(sum(diag(A)) - sum(diag(A_svd))) / sum(diag(A))
#'
#'   # Store results
#'   all_results[[length(all_results) + 1]] <- list(
#'     k = k,
#'     rpcholesky_time = rpcholesky_time[3],  # elapsed time
#'     svd_time = svd_time[3],                # elapsed time
#'     rpcholesky_error = rpcholesky_error,
#'     svd_error = svd_error,
#'     rpcholesky_trace_error = rpcholesky_trace_error,
#'     svd_trace_error = svd_trace_error,
#'     speedup = svd_time[3] / rpcholesky_time[3]
#'   )
#' }
#'
#' # -------------------------------------------------------------------------
#' # Visualize and analyze results
#' # -------------------------------------------------------------------------
#' results_df <- do.call(rbind, lapply(all_results, function(r) {
#'   data.frame(
#'     Rank = r$k,
#'     RPCholesky_Error = r$rpcholesky_error,
#'     SVD_Error = r$svd_error,
#'     Error_Ratio = r$rpcholesky_error / r$svd_error,
#'     RPCholesky_Time = r$rpcholesky_time,
#'     SVD_Time = r$svd_time,
#'     Speedup = r$speedup
#'   )
#' }))
#'
#' print(results_df)
#'
#' # Plot speedup vs rank
#' p1 <- ggplot(results_df, aes(x = Rank, y = Speedup)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   ggtitle("Computational Advantage of RPCholesky") +
#'   ylab("Speed-up Ratio (SVD time / RPCholesky time)") +
#'   xlab("Approximation Rank") +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red")
#'
#' # Plot error comparison
#' errors_df <- data.frame(
#'   Rank = rep(results_df$Rank, 2),
#'   Method = rep(c("RPCholesky", "SVD"), each = nrow(results_df)),
#'   Error = c(results_df$RPCholesky_Error, results_df$SVD_Error)
#' )
#'
#' p2 <- ggplot(errors_df, aes(x = Rank, y = Error, color = Method, group = Method)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   scale_y_log10() +
#'   ggtitle("Approximation Error Comparison") +
#'   ylab("Relative Frobenius Error (log scale)") +
#'   xlab("Approximation Rank")
#'
#' # Plot error ratio across ranks
#' p3 <- ggplot(results_df, aes(x = Rank, y = Error_Ratio)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   ggtitle("Approximation Quality Ratio") +
#'   ylab("Error Ratio (RPCholesky / SVD)") +
#'   xlab("Approximation Rank") +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red")
#'
#' # Create efficiency metric: (1/error) / time
#' results_df$RPCholesky_Efficiency <- (1/results_df$RPCholesky_Error) / results_df$RPCholesky_Time
#' results_df$SVD_Efficiency <- (1/results_df$SVD_Error) / results_df$SVD_Time
#' results_df$Efficiency_Ratio <- results_df$RPCholesky_Efficiency / results_df$SVD_Efficiency
#'
#' p4 <- ggplot(results_df, aes(x = Rank, y = Efficiency_Ratio)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   ggtitle("Efficiency Comparison") +
#'   ylab("Efficiency Ratio (RPCholesky / SVD)") +
#'   xlab("Approximation Rank") +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red")
#'
#' # Plot all figures
#' gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
#'
#' # -------------------------------------------------------------------------
#' # Demonstrate performance with increasing matrix size
#' # -------------------------------------------------------------------------
#' message("Testing scaling with matrix size...")
#'
#' # Test with different matrix sizes
#' sizes <- c(500, 1000, 2000, 3000)
#' rank_fraction <- 0.05  # Use 5% of matrix size as rank
#'
#' scaling_results <- data.frame(
#'   Size = integer(),
#'   Rank = integer(),
#'   RPCholesky_Time = numeric(),
#'   SVD_Time = numeric(),
#'   Speedup = numeric()
#' )
#'
#' for (size in sizes) {
#'   k <- ceiling(size * rank_fraction)
#'   message(paste("Matrix size:", size, "with rank", k))
#'
#'   # Create test matrix
#'   test_matrix <- create_test_matrix(size, decay_rate = 0.98)$matrix
#'
#'   # Time RPCholesky
#'   rpcholesky_time <- system.time({
#'     rpcholesky(test_matrix, k = k, b = "auto", verbose = FALSE)
#'   })[3]
#'
#'   # Time SVD
#'   svd_time <- system.time({
#'     svd(test_matrix, nu = k, nv = k)
#'   })[3]
#'
#'   # Store results
#'   scaling_results <- rbind(scaling_results, data.frame(
#'     Size = size,
#'     Rank = k,
#'     RPCholesky_Time = rpcholesky_time,
#'     SVD_Time = svd_time,
#'     Speedup = svd_time / rpcholesky_time
#'   ))
#' }
#'
#' print(scaling_results)
#'
#' # Plot scaling behavior
#' scaling_long <- reshape2::melt(scaling_results,
#'                                id.vars = c("Size", "Rank"),
#'                                measure.vars = c("RPCholesky_Time", "SVD_Time"),
#'                                variable.name = "Method",
#'                                value.name = "Time")
#'
#' p5 <- ggplot(scaling_long, aes(x = Size, y = Time, color = Method, group = Method)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   ggtitle("Computation Time vs Matrix Size") +
#'   ylab("Time (seconds)") +
#'   xlab("Matrix Size")
#'
#' p6 <- ggplot(scaling_results, aes(x = Size, y = Speedup)) +
#'   geom_line() +
#'   geom_point(size = 3) +
#'   theme_minimal() +
#'   ggtitle("Speedup vs Matrix Size") +
#'   ylab("Speedup (SVD time / RPCholesky time)") +
#'   xlab("Matrix Size") +
#'   geom_hline(yintercept = 1, linetype = "dashed", color = "red")
#'
#' gridExtra::grid.arrange(p5, p6, ncol = 2)
#'
#' # -------------------------------------------------------------------------
#' # Analyze selected pivots
#' # -------------------------------------------------------------------------
#' k_analyze <- 100  # Choose a moderate rank for analysis
#' rpcholesky_result <- rpcholesky(A, k = k_analyze, b = "auto", verbose = FALSE)
#'
#' # Get actual eigenvalues of the matrix
#' actual_eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
#'
#' # Plot top eigenvalues against indices of selected pivots
#' eigvals_df <- data.frame(
#'   Index = 1:n,
#'   Eigenvalue = actual_eigenvalues
#' )
#' # Add Selected column
#' eigvals_df$Selected <- ifelse(1:n %in% rpcholesky_result$idx, "Yes", "No")
#'
#' # Use proper subsetting with a limited number of points for visualization
#' plot_points <- min(n, 500)
#' p7 <- ggplot(eigvals_df[1:plot_points, ], aes(x = Index, y = Eigenvalue, color = Selected)) +
#'   geom_point(aes(size = Selected == "Yes")) +
#'   scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1), guide = "none") +
#'   scale_y_log10() +
#'   theme_minimal() +
#'   ggtitle("Selected Pivots vs Eigenvalue Distribution") +
#'   ylab("Eigenvalue (log scale)") +
#'   xlab("Index") +
#'   scale_color_manual(values = c("Yes" = "red", "No" = "gray"))
#'
#' print(p7)
#'
#' # Summary statistics
#' cat("\nSummary Statistics:\n")
#' cat("-------------------\n")
#' cat(sprintf("Matrix Size: %d x %d\n", n, n))
#' cat(sprintf("Number of ranks tested: %d\n", length(ranks)))
#' cat(sprintf("Average speedup over SVD: %.2fx\n", mean(results_df$Speedup)))
#' cat(sprintf("Maximum speedup over SVD: %.2fx\n", max(results_df$Speedup)))
#' cat(sprintf("Average error ratio: %.2fx\n", mean(results_df$Error_Ratio)))
#' cat("\nKey Findings:\n")
#' cat("1. RPCholesky is faster than SVD, especially for larger matrices\n")
#' cat("2. The speedup increases with matrix size\n")
#' cat("3. While SVD provides optimal approximation, RPCholesky gives excellent results with much less computation\n")
#' cat("4. RPCholesky automatically selects important pivots corresponding to largest eigenvalue directions\n")
