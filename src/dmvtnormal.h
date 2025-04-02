#ifndef DMVTNORMAL_H
#define DMVTNORMAL_H

#include <RcppArmadillo.h>

// Declare the Mahalanobis distance function
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov);

// Declare the multivariate normal density function
arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false);

#endif // DMVTNORMAL_H
