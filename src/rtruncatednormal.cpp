#define ARMA_WARN_LEVEL 1

#include "rtruncatenormal.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
// #include <RcppDist.h>
//// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

const double t1 = 0.15;
const double t2 = 2.18;
const double t3 = 0.725;
// const double t4 = 0.45;

// Helper function to perform Exponential Rejection Sampling (a, inf)
double ers_a_inf(double a) {
  double x, rho;
  do {
    x = R::rexp(1.0 / a) + a;
    rho = exp(-0.5 * (x - a) * (x - a));
  } while (R::runif(0, 1) > rho);
  return x;
}

// Helper function to perform Exponential Rejection Sampling (a, b)
double ers_a_b(double a, double b) {
  double x, rho;
  do {
    x = R::rexp(1.0 / a) + a; // R::rexp generates exponential random numbers
    rho = exp(-0.5 * (x - a) * (x - a)); // Calculate the acceptance ratio
  } while (R::runif(0, 1) > rho || x > b); // Repeat until a valid sample is generated within (a, b)
  return x;
}

// Helper function to perform Normal Rejection Sampling (a, inf)
double nrs_a_inf(double a) {
  double x = -DBL_MAX;
  while (x < a) {
    x = R::rnorm(0, 1);
  }
  return x;
}

// Helper function to perform Normal Rejection Sampling (a, b)
double nrs_a_b(double a, double b) {
  double x = -DBL_MAX;
  while (x < a || x > b) {
    x = R::rnorm(0, 1); // R::rnorm generates normal random numbers with mean 0 and standard deviation 1
  }
  return x;
}


// Helper function to perform Half-Normal Rejection Sampling (a, b)
double hnrs_a_b(double a, double b) {
  double x = a - 1.0;
  while (x < a || x > b) {
    x = R::rnorm(0, 1); // Generate normal random numbers with mean 0 and standard deviation 1
    x = std::fabs(x); // Take the absolute value
  }
  return x;
}

// Helper function to perform Uniform Rejection Sampling (a, b)
double urs_a_b(double a, double b) {
  const double phi_a = R::dnorm(a, 0.0, 1.0, false); // Calculate the normal density at 'a'
  double x = 0.0;//, u = 0.0;

  // Upper bound of normal density on [a, b]
  const double ub = (a < 0 && b > 0) ? M_1_SQRT_2PI : phi_a;
  do {
    x = R::runif(a, b); // Generate uniform random numbers between 'a' and 'b'
  } while (R::runif(0, 1) * ub > R::dnorm(x, 0, 1, false)); // Acceptance condition
  return x;
}



// Main function for left truncated normal distribution
double r_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double t4 = 0.45;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha);
  } else {
    return mean + sd * ers_a_inf(alpha);
  }
}

// Helper function to perform Right-Truncated Normal Sampling
double r_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

// Helper function to perform Truncated Normal Sampling
double r_truncnorm_my(double a, double b, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, false);
  const double phi_b = R::dnorm(beta, 0.0, 1.0, false);
  if (beta <= alpha) {
    return NA_REAL;
  } else if (alpha <= 0 && 0 <= beta) { /* 2 */
  if (phi_a <= t1 || phi_b <= t1) {   /* 2 (a) */
  return mean + sd * nrs_a_b(alpha, beta);
  } else { /* 2 (b) */
  return mean + sd * urs_a_b(alpha, beta);
  }
  } else if (alpha > 0) {      /* 3 */
  if (phi_a / phi_b <= t2) { /* 3 (a) */
  return mean + sd * urs_a_b(alpha, beta);
  } else {
    if (alpha < t3) { /* 3 (b) */
  return mean + sd * hnrs_a_b(alpha, beta);
    } else { /* 3 (c) */
  return mean + sd * ers_a_b(alpha, beta);
    }
  }
  } else {                     /* 3s */
  if (phi_b / phi_a <= t2) { /* 3s (a) */
  return mean - sd * urs_a_b(-beta, -alpha);
  } else {
    if (beta > -t3) { /* 3s (b) */
  return mean - sd * hnrs_a_b(-beta, -alpha);
    } else { /* 3s (c) */
  return mean - sd * ers_a_b(-beta, -alpha);
    }
  }
  }
}

// [[Rcpp::export]]
arma::vec rtruncnorm_cpp(int n, arma::vec a, arma::vec b, arma::vec mean, arma::vec sd) {
  arma::vec ret(n);

  for (int i = 0; i < n; ++i) {
    double ca = a[i % a.size()];
    double cb = b[i % b.size()];
    double cmean = mean[i % mean.size()];
    double csd = sd[i % sd.size()];

    if (R_IsNA(ca) || R_IsNA(cb) || R_IsNA(cmean) || R_IsNA(csd)) {
      ret[i] = NA_REAL;
    } else if (ca == R_NegInf && cb == R_PosInf) {
      ret[i] = R::rnorm(cmean, csd);
    } else if (ca == R_NegInf) {
      ret[i] = r_righttruncnorm(cb, cmean, csd);
    } else if (cb == R_PosInf) {
      ret[i] = r_lefttruncnorm(ca, cmean, csd);
    } else {
      ret[i] = r_truncnorm_my(ca, cb, cmean, csd);
    }
  }

  return ret;
}


// // [[Rcpp::export]]
// Rcpp::NumericVector rtruncnorm_function(const int n, Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu, Rcpp::NumericVector sigma) {
//   Rcpp::NumericVector x(n);
//   for (int i = 0; i < n; i++) {
//     NumericVector result = rtruncnorm(1, mu[i], sigma[i], a[i], b[i]);
//     x[i] = result[0];
//     }
//   return x;
// }
//
// Rcpp::NumericVector rtruncnorm(const int n, const double mu,
// const double sigma, const double a, const double b)
