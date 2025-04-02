#ifndef E_INTEGRAL_H
#define E_INTEGRAL_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtruncatenormal.h"
#include "basic_functions.h"
#include "dmvtnormal.h"
#include <chrono>
using namespace Rcpp;

arma::mat create_f_vector_matrix(const arma::vec& f_vector, int n_cols);
arma::vec col_means(const arma::mat& X);

List MC_int_xy_all_cpp(const arma::vec& s_i, const arma::vec& y_i, double mu_tte, double sd_tte,
                       const arma::vec& mu_y, double sigma_y_sq, const arma::vec& mu_r, const arma::mat& sigma_r,
                       int sample_J, double x_min, bool obs = false);
#endif // E_INTEGRAL_H
