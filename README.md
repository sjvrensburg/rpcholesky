# rpcholesky

## Overview

The `rpcholesky` package implements an accelerated randomly pivoted Cholesky decomposition for low‐rank approximation of positive semidefinite matrices. It provides functions for computing the decomposition, evaluating approximation error, selecting an optimal rank, comparing the approximation with singular value decomposition (SVD), and analyzing pivot selection.

## Installation

With the devtools package:

```r
devtools::install_local("path/to/rpcholesky")
```

The package requires Rcpp and RcppArmadillo.

## Usage

Load the package by running:

```r
library(rpcholesky)
```

## Functions

- **rpcholesky(A, k, b = "auto", stoptol = 1e-10, verbose = FALSE)**  
  Computes the accelerated randomly pivoted Cholesky decomposition of a positive semidefinite matrix `A` to produce a low‐rank factorization of target rank `k`.

- **rpcholesky_error(A, approx, factor = TRUE)**  
  Calculates various error metrics (Frobenius norm, spectral norm, trace error, and maximum entry error) between the original matrix `A` and its low‐rank approximation.

- **kernel_alignment(K, K_approx, factor = TRUE)**  
  Computes the kernel alignment between the original kernel matrix `K` and its approximated version.

- **find_optimal_rank(A, error_threshold = 0.01, error_type = "frobenius", max_rank = min(nrow(A), 100), step_size = 5, binary_search = TRUE, verbose = FALSE)**  
  Determines the minimum rank for the low‐rank approximation of `A` such that the relative error is below a specified threshold.

- **estimate_rank_eigenvalues(A, method = "elbow", threshold = 0.95, max_rank = min(nrow(A), 100), plot = TRUE)**  
  Estimates an optimal rank based on the decay of the eigenvalues of `A` using either the "elbow" or "threshold" method.

- **profile_rank_errors(A, ranks = seq(5, min(nrow(A), 100), by = 5), error_types = c("frobenius", "trace", "max"), verbose = TRUE, plot = TRUE)**  
  Profiles approximation error over a range of ranks for different error metrics.

- **compare_with_svd(A, ranks = seq(5, min(nrow(A), 50), by = 5), metrics = c("error", "time"), time_benchmark = TRUE, plot = TRUE)**  
  Compares the low‐rank approximations obtained from `rpcholesky` and SVD, reporting errors and computation times.

- **analyze_pivots(A, result, compute_eigendecomp = TRUE, plot = TRUE)**  
  Analyzes the pivot selection from the decomposition and compares properties of the selected pivots with those of the original matrix.

## Examples

### Example 1: Low‐Rank Approximation

```r
# Create a positive semidefinite matrix
set.seed(123)
n <- 100
A <- matrix(rnorm(n * n), n, n)
A <- A %*% t(A)

# Compute the low‐rank approximation with target rank 10
result <- rpcholesky(A, k = 10, b = "auto", stoptol = 1e-10, verbose = TRUE)
G <- result$G

# Reconstruct the approximated matrix
A_approx <- t(G) %*% G

# Compute error metrics
errors <- rpcholesky_error(A, result)
print(errors)
```

### Example 2: Optimal Rank Selection

```r
# Find the optimal rank such that the relative Frobenius error is below 0.05
optimal <- find_optimal_rank(A, error_threshold = 0.05, error_type = "frobenius", verbose = TRUE)
print(optimal)
```

### Example 3: Eigenvalue-Based Rank Estimation

```r
# Estimate the optimal rank using the eigenvalue decay (elbow method)
estimation <- estimate_rank_eigenvalues(A, method = "elbow", plot = FALSE)
print(estimation$optimal_rank)
```

## Package Structure

- **src/**: C++ source files and makefiles used for compiling the C++ components.
- **R/**: R scripts that define the package functions and interface.
- **inst/**: Includes header files and example scripts.
- **man/**: Documentation files for the functions.
- **DESCRIPTION**: Package metadata.
- **NAMESPACE**: Export and import directives.

## References and Acknowledgements

The implementation in this package is based on the following work:

- **Chen, Y., Epperly, E. N., Tropp, J. A., & Webber, R. J. (2023).** *Randomly pivoted Cholesky: Practical approximation of a kernel matrix with few entry evaluations.* [arXiv:2207.06503](https://arxiv.org/abs/2207.06503).

- **Epperly, E. N., Tropp, J. A., & Webber, R. J. (2024).** *Embrace rejection: Kernel matrix approximation by accelerated randomly pivoted Cholesky.* Manuscript in preparation.

This code is also based on the implementation available at [https://github.com/eepperly/Randomly-Pivoted-Cholesky](https://github.com/eepperly/Randomly-Pivoted-Cholesky).

## License

The `rpcholesky` package is licensed under the GNU General Public License version 3 (GPL-3).
