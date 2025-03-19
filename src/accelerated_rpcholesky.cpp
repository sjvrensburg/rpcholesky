#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "../inst/include/rpcholesky.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rejection_cholesky(arma::mat H) {
  // Forward the call to the inline implementation in the namespace
  return rpcholesky::rejection_cholesky(H);
}

// [[Rcpp::export]]
List accelerated_rpcholesky(const arma::mat &A, int k, SEXP b_in, double stoptol, bool verbose) {
  // Forward the call to the inline implementation in the namespace
  return rpcholesky::accelerated_rpcholesky(A, k, b_in, stoptol, verbose);
}
