#include "E_integral.h"
#include "rtruncatenormal.h"
#include "basic_functions.h"
#include "dmvtnormal.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo, mvtnorm)]]

using namespace Rcpp;
using namespace std::chrono;

arma::mat create_f_vector_matrix(const arma::vec& f_vector, int n_cols) {
  arma::mat f_vector_matrix = repmat(f_vector, 1, n_cols);

  return f_vector_matrix;
}


arma::vec col_means(const arma::mat& X) {
  arma::rowvec out = arma::mean(X, 0);
  return out.t();
}
//
// // [[Rcpp::export]]
// double MC_int_xy_cpp(const arma::vec &s_i, const arma::vec &y_i, double mu_tte, double sd_tte,
//                      const arma::vec &mu_y, double sigma_y_sq, const arma::vec &mu_r, const arma::mat &sigma_r,
//                      int sample_J, double x_min, bool obs = false)
// {
//   // Function dmvnorm("dmvnorm", Environment::namespace_env("mvtnorm"));
//
//   int p = s_i.n_elem;
//
//   arma::vec t_i_j_vector(sample_J);
//   if (obs)
//   {
//     t_i_j_vector.fill(x_min);
//   }
//   else
//   {
//     arma::vec lower_bounds(sample_J);
//     lower_bounds.fill(x_min);
//     arma::vec upper_bounds(sample_J);
//     upper_bounds.fill(INFINITY);
//     arma::vec means(sample_J);
//     means.fill(mu_tte);
//     arma::vec sds(sample_J);
//     sds.fill(sd_tte);
//
//     t_i_j_vector = rtruncnorm(sample_J, lower_bounds, upper_bounds, means, sds);
//   }
//
//   double mu_omega = mu_r[0];
//   arma::vec mu_b = mu_r.subvec(1, mu_r.n_elem - 1);
//   double sigma_omega_sq = sigma_r(0, 0);
//   arma::vec sigma_bomega = sigma_r.submat(1, 0, mu_r.n_elem - 1, 0);
//   arma::mat sigma_b = sigma_r.submat(1, 1, mu_r.n_elem - 1, mu_r.n_elem - 1);
//
//   arma::vec f_vector(sample_J);
//
//   arma::vec lower_bounds(sample_J); // Vector of zeros
//   lower_bounds.fill(0.0);
//   arma::vec upper_bounds(sample_J);
//   for (int j = 0; j < sample_J; j++)
//   {
//     double t_i_j = t_i_j_vector(j);
//     upper_bounds.fill(exp(t_i_j)); // Vector filled with exp(t_i_j)
//   }
//   arma::vec means(sample_J); // Vector filled with mu_omega
//   means.fill(mu_omega);
//   arma::vec sds(sample_J);
//   sds.fill(sqrt(sigma_omega_sq)); // Vector filled with sqrt(sigma_omega_sq)
//
//   arma::vec omega_i_j_vector = rtruncnorm(sample_J, lower_bounds, upper_bounds, means, sds);
//
//   for (int j = 0; j < sample_J; j++)
//   {
//
//     double y = omega_i_j_vector(j);
//     arma::mat Z_i_value = Z_i_cpp(y, s_i);
//
//     arma::vec mean_vec = mu_y + Z_i_value * (mu_b + sigma_bomega / sigma_omega_sq * (y - mu_omega));
//
//     arma::mat sigma_matrix = sigma_y_sq * arma::eye(p, p) + Z_i_value * (sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) * Z_i_value.t();
//
//     arma::rowvec y_i_matrix = y_i.t();
//     arma::rowvec mean_vec_row = mean_vec.t();
//
//     arma::mat density = dmvnorm_arma(y_i_matrix, mean_vec_row, sigma_matrix, false);
//     f_vector(j) = density[0];
//   }
//
//   return mean(f_vector);
// }


// [[Rcpp::export]]
List MC_int_xy_all_cpp(const arma::vec& s_i, const arma::vec& y_i, double mu_tte, double sd_tte,
                       const arma::vec& mu_y, double sigma_y_sq, const arma::vec& mu_r, const arma::mat& sigma_r,
                       int sample_J, double x_min, bool obs) {
  // record time
  // Time point variables
  auto start = steady_clock::now();
  auto end = steady_clock::now();

  start = steady_clock::now();

  int p = s_i.n_elem;
  double mu_omega = mu_r[0];
  arma::vec mu_b = mu_r.subvec(1, mu_r.n_elem - 1);
  double sigma_omega_sq = sigma_r(0, 0);
  arma::vec sigma_bomega = sigma_r.submat(1, 0, mu_r.n_elem - 1, 0);
  arma::mat sigma_b = sigma_r.submat(1, 1, mu_r.n_elem - 1, mu_r.n_elem - 1);

  arma::vec t_i_j_vector(sample_J);
  if (obs) {
    t_i_j_vector.fill(x_min);
  } else {
    arma::vec lower_bounds(sample_J);
    lower_bounds.fill(x_min);
    arma::vec upper_bounds(sample_J);
    upper_bounds.fill(INFINITY);
    arma::vec means(sample_J);
    means.fill(mu_tte);
    arma::vec sds(sample_J);
    sds.fill(sd_tte);

    t_i_j_vector = rtruncnorm_cpp(sample_J, lower_bounds, upper_bounds, means, sds);
  }

  arma::vec f_vector(sample_J);

  arma::vec lower_bounds(sample_J); // Vector of zeros
  lower_bounds.fill(0.0);
  arma::vec upper_bounds(sample_J);
  for (int j = 0; j < sample_J; j++)
  {
    double t_i_j = t_i_j_vector(j);
    upper_bounds.fill(exp(t_i_j)); // Vector filled with exp(t_i_j)
  }
  arma::vec means(sample_J); // Vector filled with mu_omega
  means.fill(mu_omega);
  arma::vec sds(sample_J);
  sds.fill(sqrt(sigma_omega_sq)); // Vector filled with sqrt(sigma_omega_sq)

  arma::vec omega_i_j_vector = rtruncnorm_cpp(sample_J, lower_bounds, upper_bounds, means, sds);

  // record time
  end = steady_clock::now();
  // Rcout << "Sampling time: " << duration<double>(end - start).count() << " seconds\n";

  start = steady_clock::now();


  arma::mat E_b_i_matrix(sample_J, 3, arma::fill::zeros);
  arma::mat E_b_i_b_i_T_matrix(sample_J, 9, arma::fill::zeros);
  arma::mat E_Z_i_b_i_matrix(sample_J, p, arma::fill::zeros);
  arma::mat E_Z_i_T_Z_i_b_i_b_i_T_matrix(sample_J, 9, arma::fill::zeros);

  for (int j = 0; j < sample_J; ++j) {
    double y = omega_i_j_vector(j);
    arma::mat Z_i_result = Z_i_cpp(y, s_i);

    arma::vec mean_vec = mu_y + Z_i_result * (mu_b + sigma_bomega / sigma_omega_sq * (y - mu_omega));
    arma::mat sigma_matrix = sigma_y_sq * arma::eye(p, p) + Z_i_result * (sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) * Z_i_result.t();
    arma::mat inv_sigma_matrix = inv(sigma_matrix);

    arma::rowvec y_i_matrix = y_i.t();
    arma::rowvec mean_vec_row = mean_vec.t();

    arma::mat density = dmvnorm_arma(y_i_matrix, mean_vec_row, sigma_matrix, false);
    f_vector(j) = density[0];

    arma::vec E_b_i = mu_b + sigma_bomega / sigma_omega_sq * (y - mu_omega) +
      (sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) * Z_i_result.t() *
      inv_sigma_matrix * (y_i - mean_vec);

    arma::mat Var_b_i = (sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) -
      (sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) * Z_i_result.t() *
      inv_sigma_matrix * trans((sigma_b - 1 / sigma_omega_sq * sigma_bomega * sigma_bomega.t()) * Z_i_result.t());
    arma::mat E_b_i_b_i_T = E_b_i * E_b_i.t() + Var_b_i;

    E_b_i_matrix.row(j) = E_b_i.t();
    E_b_i_b_i_T_matrix.row(j) = vectorise(E_b_i_b_i_T).t();
    E_Z_i_b_i_matrix.row(j) = (Z_i_result * E_b_i).t();
    E_Z_i_T_Z_i_b_i_b_i_T_matrix.row(j) = vectorise(Z_i_result.t() * Z_i_result * E_b_i_b_i_T).t();
  }
  // record time
  end = steady_clock::now();
  // Rcout << "Integral time: " << duration<double>(end - start).count() << " seconds\n";
  start = steady_clock::now();


  double E_denominator = mean(f_vector);
  double E_ti = mean(t_i_j_vector % f_vector);
  double E_ti_sq = mean(square(t_i_j_vector) % f_vector);
  double E_g0_ti = mean(g_0_cpp(t_i_j_vector, mu_omega, sqrt(sigma_omega_sq)) % f_vector);
  double E_g1_ti = mean(g_1_cpp(t_i_j_vector, mu_omega, sqrt(sigma_omega_sq)) % f_vector);
  double E_g2_ti = mean(g_2_cpp(t_i_j_vector, mu_omega, sqrt(sigma_omega_sq)) % f_vector);

  double E_omegai = mean(omega_i_j_vector % f_vector);
  double E_omegai_sq = mean(square(omega_i_j_vector) % f_vector);

  arma::mat f_vector_matrix = create_f_vector_matrix(f_vector, E_b_i_matrix.n_cols);
  arma::mat f_vector_matrix_9 = create_f_vector_matrix(f_vector, 9);
  arma::mat omega_i_j_matrix = create_f_vector_matrix(omega_i_j_vector, E_b_i_matrix.n_cols);
  arma::mat f_vector_matrix_p = create_f_vector_matrix(f_vector, E_Z_i_b_i_matrix.n_cols);

  arma::vec E_bi = col_means(E_b_i_matrix % f_vector_matrix);
  arma::vec E_ri(4);
  E_ri(0) = E_omegai;
  E_ri.subvec(1, 3) = E_bi;

  arma::mat E_bi_bi_T = reshape(col_means(E_b_i_b_i_T_matrix % f_vector_matrix_9), 3, 3);

  arma::vec E_omegai_bi = col_means(E_b_i_matrix % f_vector_matrix % omega_i_j_matrix);

  arma::mat E_ri_riT(4, 4, arma::fill::zeros);
  E_ri_riT(0, 0) = E_omegai_sq;
  E_ri_riT.submat(0, 1, 0, 3) = E_omegai_bi.t();
  E_ri_riT.submat(1, 0, 3, 0) = E_omegai_bi;
  E_ri_riT.submat(1, 1, 3, 3) = E_bi_bi_T;

  arma::vec E_Z_i_b_i = col_means(E_Z_i_b_i_matrix % f_vector_matrix_p);
  arma::mat E_Z_i_T_Z_i_b_i_b_i_T = reshape(col_means(E_Z_i_T_Z_i_b_i_b_i_T_matrix % f_vector_matrix_9), 3, 3);

  double E_Z_i_b_i_sq = trace(E_Z_i_T_Z_i_b_i_b_i_T);
  // record time
  end = steady_clock::now();
  // Rcout << "Final computations time: " << duration<double>(end - start).count() << " seconds\n";

  return List::create(Named("E_denominator") = E_denominator,
                      Named("E_ti") = E_ti, Named("E_ti_sq") = E_ti_sq,
                      Named("E_g0_ti") = E_g0_ti, Named("E_g1_ti") = E_g1_ti,
                      Named("E_g2_ti") = E_g2_ti, Named("E_ri") = E_ri,
                      Named("E_ri_riT") = E_ri_riT, Named("E_Z_i_b_i") = E_Z_i_b_i,
                      Named("E_Z_i_b_i_sq") = E_Z_i_b_i_sq);
}

