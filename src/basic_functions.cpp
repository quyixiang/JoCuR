#define ARMA_WARN_LEVEL 1
#include "basic_functions.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Function to compute g_0
// // [[Rcpp::export]]
// double g_0(double t_i, double mu_omega, double sigma_omega) {
//   double exp_t_i = std::exp(t_i);
//   double z1 = (exp_t_i - mu_omega) / sigma_omega;
//   double z2 = -mu_omega / sigma_omega;
//   return std::log(R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0));
// }

// [[Rcpp::export]]
arma::vec g_0_cpp(const arma::vec& t, double mu_omega, double sigma_omega) {
  arma::vec result(t.n_elem);
  for (size_t i = 0; i < t.n_elem; i++) {
    double exp_t_i = std::exp(t(i));
    double z1 = (exp_t_i - mu_omega) / sigma_omega;
    double z2 = -mu_omega / sigma_omega;
    double pnorm_diff = R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0);
    result(i) = std::log(pnorm_diff);
  }
  return result;
}

// // Function to compute g_1
// // [[Rcpp::export]]
// double g_1(double t_i, double mu_omega, double sigma_omega) {
//   double exp_t_i = std::exp(t_i);
//   double z1 = (exp_t_i - mu_omega) / sigma_omega;
//   double z2 = -mu_omega / sigma_omega;
//   double pnorm_diff = R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0);
//   double dnorm_diff = R::dnorm(z1, 0.0, 1.0, 0) - R::dnorm(z2, 0.0, 1.0, 0);
//   return dnorm_diff / (pnorm_diff * sigma_omega);
// }

// [[Rcpp::export]]
arma::vec g_1_cpp(const arma::vec& t, double mu_omega, double sigma_omega) {
  arma::vec result(t.n_elem);
  for (size_t i = 0; i < t.n_elem; i++) {
    double exp_t_i = std::exp(t(i));
    double z1 = (exp_t_i - mu_omega) / sigma_omega;
    double z2 = -mu_omega / sigma_omega;
    double pnorm_diff = R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0);
    double dnorm_diff = R::dnorm(z1, 0.0, 1.0, 0) - R::dnorm(z2, 0.0, 1.0, 0);
    result(i) = dnorm_diff / (pnorm_diff * sigma_omega);
  }
  return result;
}

// Function to compute g_2
// // [[Rcpp::export]]
// double g_2(double t_i, double mu_omega, double sigma_omega) {
//   double exp_t_i = std::exp(t_i);
//   double z1 = (exp_t_i - mu_omega) / sigma_omega;
//   double z2 = -mu_omega / sigma_omega;
//   double pnorm_diff = R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0);
//   double dnorm_z1 = R::dnorm(z1, 0.0, 1.0, 0);
//   double dnorm_z2 = R::dnorm(z2, 0.0, 1.0, 0);
//   return (dnorm_z1 * -z1 - dnorm_z2 * z2) / pnorm_diff;
// }

// [[Rcpp::export]]
arma::vec g_2_cpp(const arma::vec& t, double mu_omega, double sigma_omega) {
  arma::vec result(t.n_elem);
  for (size_t i = 0; i < t.n_elem; i++) {
    double exp_t_i = std::exp(t(i));
    double z1 = (exp_t_i - mu_omega) / sigma_omega;
    double z2 = -mu_omega / sigma_omega;
    double pnorm_diff = R::pnorm(z1, 0.0, 1.0, 1, 0) - R::pnorm(z2, 0.0, 1.0, 1, 0);
    double dnorm_z1 = R::dnorm(z1, 0.0, 1.0, 0);
    double dnorm_z2 = R::dnorm(z2, 0.0, 1.0, 0);
    result(i) = (dnorm_z1 * (-z1) - dnorm_z2 * (-z2)) / pnorm_diff;
  }
  return result;
}

// Function to compute Z_i
// [[Rcpp::export]]
arma::mat Z_i_cpp(double y, const arma::vec &visittime_vec_i) {
  int n = visittime_vec_i.n_elem;
  arma::vec Delta_i = visittime_vec_i - y;

  arma::vec sec_col = Delta_i % (Delta_i <= 0);
  arma::vec third_col = Delta_i % (Delta_i > 0);

  arma::mat result(n, 3);
  result.col(0).ones();
  result.col(1) = sec_col;
  result.col(2) = third_col;

  return result;
}


