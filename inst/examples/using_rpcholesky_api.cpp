#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(rpcholesky)]]

// Include the combined header that has both inline implementations
// and exported functions
#include <rpcholesky.h>

using namespace Rcpp;
using namespace arma;

//' Low-rank matrix approximation using rpcholesky
 //'
 //' @param A A positive semidefinite matrix
 //' @param k Target rank for the approximation
 //' @param compute_error Whether to compute approximation error
 //'
 //' @return A list containing the approximation results
 //' @export
 // [[Rcpp::export]]
 List low_rank_approximation(arma::mat A, int k, bool compute_error = true) {
   // Create a low-rank approximation using rpcholesky
   List rpchol_result = rpcholesky::accelerated_rpcholesky(
     A, k, wrap("auto"), 1e-10, false);

   // Extract the low-rank factor G
   arma::mat G = rpchol_result["G"];
   arma::uvec idx = rpchol_result["idx"];

   // Construct the approximation A_approx = G^T * G
   arma::mat A_approx = G.t() * G;

   // Create result list
   List result = List::create(
     Named("approximation") = A_approx,
     Named("factor") = G,
     Named("selected_indices") = idx,
     Named("rank") = k
   );

   // Optionally compute approximation error
   if (compute_error) {
     double frob_norm_original = norm(A, "fro");
     double frob_error = norm(A - A_approx, "fro");
     double relative_error = frob_error / frob_norm_original;

     result["frobenius_error"] = frob_error;
     result["relative_error"] = relative_error;
   }

   return result;
 }

 //' Demonstrate using the helper functions from the rpcholesky namespace
 //'
 //' @param weights A vector of weights
 //' @param n_samples Number of samples to draw
 //'
 //' @return A vector of indices sampled according to weights
 //' @export
 // [[Rcpp::export]]
 arma::uvec weighted_sampling_demo(arma::vec weights, int n_samples) {
   // Use the inline implementation directly from the namespace
   return rpcholesky::sample_with_replacement_optimized(
     weights.n_elem, n_samples, weights);
 }
