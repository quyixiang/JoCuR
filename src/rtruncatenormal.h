#ifndef RTRUNCATENORMAL_H
#define RTRUNCATENORMAL_H

#include <RcppArmadillo.h>

// Declare the function
arma::vec rtruncnorm_cpp(int n, arma::vec a, arma::vec b, arma::vec mean, arma::vec sd);
Rcpp::NumericVector rtruncnorm_function(const int n, Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu, Rcpp::NumericVector sigma);
#endif // RTRUNCATENORMAL_H
