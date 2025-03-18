#ifndef ACCELERATED_RPCHOLESKY_H
#define ACCELERATED_RPCHOLESKY_H

// Forward declarations
#include <RcppArmadillo.h>

/**
 * Optimized weighted sampling with replacement
 *
 * @param n Size of the population to sample from
 * @param b Number of samples to draw
 * @param weights Vector of weights for sampling probabilities
 * @return Vector of indices sampled
 */
arma::uvec sample_with_replacement_optimized(const int n, const int b, const arma::vec &weights);

/**
 * Refined rejection Cholesky decomposition
 *
 * @param H Positive semidefinite matrix
 * @return List containing the Cholesky factor L and accepted pivot indices
 */
Rcpp::List rejection_cholesky(arma::mat H);

/**
 * Accelerated randomly pivoted Cholesky decomposition
 *
 * @param A Input positive semidefinite matrix
 * @param k Target rank for approximation
 * @param b_in Batch size (integer or "auto")
 * @param stoptol Early stopping tolerance relative to original trace
 * @param verbose Whether to print progress information
 * @return List containing low-rank factor G, selected pivot indices, and corresponding rows
 */
Rcpp::List accelerated_rpcholesky(const arma::mat &A, int k, SEXP b_in, double stoptol, bool verbose);

#endif // ACCELERATED_RPCHOLESKY_H
