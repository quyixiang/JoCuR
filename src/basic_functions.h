#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include <RcppArmadillo.h>

// Function prototypes
arma::vec g_0_cpp(const arma::vec& t, double mu_omega, double sigma_omega);
arma::vec g_1_cpp(const arma::vec& t, double mu_omega, double sigma_omega);
arma::vec g_2_cpp(const arma::vec& t, double mu_omega, double sigma_omega);
arma::mat Z_i_cpp(double y, const arma::vec &visittime_vec_i);

#endif // BASIC_FUNCTIONS_H
