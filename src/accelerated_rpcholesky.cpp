// [[Rcpp::depends(RcppArmadillo)]]
#include "accelerated_rpcholesky.h"
#include <chrono>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

//------------------------------------------------------------------------------
// Optimized Weighted Sampling with Replacement
//------------------------------------------------------------------------------
arma::uvec sample_with_replacement_optimized(const int n, const int b, const arma::vec &weights) {
  arma::vec cum_weights = arma::cumsum(weights);
  double total = cum_weights(cum_weights.n_elem - 1);
  arma::vec r = arma::randu(b) * total;
  arma::uvec samples(b);
  double* cum_ptr = cum_weights.memptr();
  for (int i = 0; i < b; i++) {
    // Use std::upper_bound on the raw pointer
    const double* pos = std::upper_bound(cum_ptr, cum_ptr + cum_weights.n_elem, r(i));
    samples(i) = pos - cum_ptr;
  }
  return samples;
}

//------------------------------------------------------------------------------
// Refined Rejection Cholesky Decomposition
//------------------------------------------------------------------------------
List rejection_cholesky(arma::mat H) {
  int b = H.n_rows;
  if (H.n_rows != H.n_cols)
    stop("rejection_cholesky requires a square matrix");
  if (arma::trace(H) <= 0)
    stop("rejection_cholesky requires a strictly positive trace");

  // Save the original diagonal of H.
  arma::vec u = H.diag();
  std::vector<int> accepted;
  arma::mat L = arma::zeros<arma::mat>(b, b);

  for (int j = 0; j < b; j++) {
    double r = arma::randu();
    // Accept pivot j if r * u[j] < H(j,j)
    if (r * u(j) < H(j,j)) {
      accepted.push_back(j);
      double pivot = H(j,j);
      double sqrt_pivot = std::sqrt(pivot);
      // Compute L[j:, j] = H[j:, j] / sqrt(H(j,j))
      L(arma::span(j, b-1), j) = H(arma::span(j, b-1), j) / sqrt_pivot;
      // Update trailing submatrix H[(j+1):, (j+1):] -= outer(L[(j+1):, j], L[(j+1):, j])
      if (j + 1 < b) {
        arma::vec L_col = L(arma::span(j+1, b-1), j);
        H(arma::span(j+1, b-1), arma::span(j+1, b-1)) -= L_col * L_col.t();
      }
    }
  }

  // Convert accepted indices to an arma::uvec using Armadillo's conversion.
  arma::uvec idx = arma::conv_to<arma::uvec>::from(accepted);
  // Extract the accepted submatrix from L: L(idx, idx)
  arma::mat L_accepted;
  if (idx.n_elem > 0) {
    L_accepted = L.submat(idx, idx);
  }

  return List::create(Named("L") = L_accepted,
                      Named("idx") = idx);
}

//------------------------------------------------------------------------------
// Refined Accelerated Randomized Pivoted Cholesky Decomposition
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List accelerated_rpcholesky(const arma::mat &A, int k, SEXP b_in, double stoptol, bool verbose) {
  if (A.n_rows != A.n_cols)
    stop("Matrix A must be square.");
  int n = A.n_rows;
  if (k <= 0 || k > n)
    stop("Invalid rank k.");

  // Initialize the diagonal residual and original trace.
  arma::vec diags = A.diag();
  double orig_trace = arma::sum(diags);
  if (orig_trace <= 0)
    stop("Matrix A must have a positive trace.");

  // Determine batch size b: if "auto", set b = ceil(k/10)
  int b;
  bool auto_b = false;
  if (TYPEOF(b_in) == STRSXP) {
    std::string b_str = as<std::string>(b_in);
    if (b_str == "auto") {
      b = std::ceil(k / 10.0);
      auto_b = true;
    } else {
      stop("Invalid b value. Use an integer or 'auto'.");
    }
  } else {
    b = as<int>(b_in);
  }

  // Initialize output matrices:
  // G (k x n) will store the low-rank factor rows,
  // rows (k x n) will store the corresponding rows of A,
  // arr_idx (length k) stores the global pivot indices.
  arma::mat G = arma::zeros<arma::mat>(k, n);
  arma::mat rows = arma::zeros<arma::mat>(k, n);
  arma::uvec arr_idx = arma::zeros<arma::uvec>(k);
  int counter = 0;

  // Timing variables for adaptive batching.
  double total_rej_time = 0.0;
  double total_proc_time = 0.0;
  std::chrono::high_resolution_clock::time_point start_time, end_time;

  while (counter < k) {
    // Sample candidate indices with replacement from 0 to n-1 using the optimized function.
    double sum_diags = arma::sum(diags);
    if (sum_diags <= 0)
      break;
    arma::vec probs = diags / sum_diags;
    arma::uvec cand = sample_with_replacement_optimized(n, b, probs);

    // Form submatrix H = A(cand, cand) minus the correction from previously accepted pivots.
    arma::mat A_cand = A.submat(cand, cand);
    arma::mat H;
    if (counter > 0) {
      arma::uvec row_indices = arma::regspace<arma::uvec>(0, counter - 1);
      arma::mat G_sub = G.submat(row_indices, cand);
      H = A_cand - G_sub.t() * G_sub;
    } else {
      H = A_cand;
    }

    if (auto_b) {
      start_time = std::chrono::high_resolution_clock::now();
    }

    // Run refined rejection sampling on H.
    List rej_res = rejection_cholesky(H);
    arma::mat L = as<arma::mat>(rej_res["L"]);
    arma::uvec accepted = as<arma::uvec>(rej_res["idx"]);

    if (auto_b) {
      end_time = std::chrono::high_resolution_clock::now();
    }
    double rej_time = auto_b ? std::chrono::duration<double>(end_time - start_time).count() : 0.0;
    total_rej_time += rej_time;

    int num_sel = accepted.n_elem;
    // Limit the number of accepted pivots if more than needed.
    if (num_sel > (k - counter)) {
      num_sel = k - counter;
      accepted = accepted.head(num_sel);
      if (L.n_rows >= (unsigned)num_sel && L.n_cols >= (unsigned)num_sel) {
        L = L.submat(0, 0, num_sel - 1, num_sel - 1);
      }
    }

    // Map accepted candidate indices back to global indices.
    arma::uvec accepted_cand = cand.elem(accepted);

    if (auto_b) {
      start_time = std::chrono::high_resolution_clock::now();
    }

    // Save the accepted global pivot indices.
    arr_idx.subvec(counter, counter + num_sel - 1) = accepted_cand;
    // Save the corresponding rows of A.
    rows.rows(counter, counter + num_sel - 1) = A.rows(accepted_cand);

    // Compute the update for G:
    arma::mat G_update;
    if (counter > 0) {
      arma::uvec row_indices_accepted = arma::regspace<arma::uvec>(0, counter - 1);
      arma::mat G_corr = G.submat(row_indices_accepted, accepted_cand).t()
        * G.rows(row_indices_accepted);
      G_update = A.rows(accepted_cand) - G_corr;
    } else {
      G_update = A.rows(accepted_cand);
    }

    // Solve L * X = G_update for X.
    arma::mat X;
    bool solved = arma::solve(X, L, G_update);
    if (!solved)
      stop("Linear system could not be solved.");
    // Store the result in G.
    G.rows(counter, counter + num_sel - 1) = X;

    // Update diags: subtract the contribution of the new G rows.
    arma::rowvec sq_sum = arma::sum(arma::square(X), 0);
    diags -= sq_sum.t();
    diags.elem(arma::find(diags < 0)).zeros();

    if (auto_b) {
      end_time = std::chrono::high_resolution_clock::now();
    }
    double proc_time = auto_b ? std::chrono::duration<double>(end_time - start_time).count() : 0.0;
    total_proc_time += proc_time;

    // Adaptive batch size update.
    if (auto_b && rej_time > 0) {
      double target = std::ceil(b * proc_time / (4.0 * rej_time));
      double candidate1 = target;
      double candidate2 = std::ceil(1.5 * b);
      double candidate3 = std::ceil(k / 3.0);
      double m = std::min({candidate1, candidate2, candidate3});
      double candidate4 = std::ceil(b / 3.0);
      b = std::max({10.0, m, candidate4});
      b = static_cast<int>(b);
    }

    counter += num_sel;

    if (stoptol > 0 && arma::sum(diags) <= stoptol * orig_trace) {
      G = G.rows(0, counter - 1);
      rows = rows.rows(0, counter - 1);
      arr_idx = arr_idx.subvec(0, counter - 1);
      break;
    }

    if (verbose) {
      Rcout << "Accepted " << num_sel << " / " << b << std::endl;
    }
  }

  return List::create(Named("G") = G,
                      Named("idx") = arr_idx,
                      Named("rows") = rows);
}
