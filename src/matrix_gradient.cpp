#include <RcppArmadillo.h>
using namespace Rcpp;
#include "matrix_gradient.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Function to calculate dSigmadsigmaij
// [[Rcpp::export]]
arma::mat dSigmadsigmaij_cpp(int i, int j, arma::mat Sigma) {
  int dimn = Sigma.n_rows;
  arma::sp_mat J_ij(dimn, dimn);
  J_ij(i, j) = 1;
  arma::mat result = arma::mat(J_ij) + arma::trans(arma::mat(J_ij)) - arma::mat(J_ij) * arma::mat(J_ij);
  return result;
}

// Function to calculate dldsigmaij
// [[Rcpp::export]]
double dldsigmaij_cpp(int i, int j, arma::mat Sigma, arma::mat E_ri_riT) {
  arma::mat partial = dSigmadsigmaij_cpp(i, j, Sigma);
  arma::mat Sigma_inv = arma::inv(Sigma);
  double a = arma::trace(Sigma_inv * partial);
  double b = arma::trace(-Sigma_inv * partial * Sigma_inv * E_ri_riT);
  return -0.5 * (a + b);
}


// Helper function to compute the vech (half-vectorization) of a matrix
//[[Rcpp::export]]
arma::vec vech_cpp(const arma::mat& mat) {
  arma::uvec indices = arma::trimatl_ind(arma::size(mat));
  return mat(indices);
}

// Function to compute dldvechSigma
// [[Rcpp::export]]
arma::vec dldvechSigma_cpp(arma::mat Sigma, arma::mat E_ri_riT) {
  arma::mat Sigma_inv = arma::inv(Sigma);
  arma::mat tmp = Sigma_inv - Sigma_inv * E_ri_riT * Sigma_inv;

  // Symmetrize tmp
  tmp = arma::trans(tmp) + tmp - arma::diagmat(tmp);

  arma::vec result = vech_cpp(tmp);
  return -0.5 * result;
}

// [[Rcpp::export]]
arma::mat duplication_matrix_cpp(const int &n) {
  arma::mat out((n*(n+1))/2, n*n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      arma::vec u((n*(n+1))/2, arma::fill::zeros);
      u(j*n+i-((j+1)*j)/2) = 1.0;
      arma::mat T(n,n, arma::fill::zeros);
      T(i,j) = 1.0;
      T(j,i) = 1.0;
      out += u * arma::trans(arma::vectorise(T));
    }
  }
  return out.t();
}

// [[Rcpp::export]]
arma::mat elimination_matrix_cpp(const int &n) {
  arma::mat out((n*(n+1))/2, n*n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    arma::rowvec e_j(n, arma::fill::zeros);
    e_j(j) = 1.0;
    for (int i = j; i < n; ++i) {
      arma::vec u((n*(n+1))/2, arma::fill::zeros);
      u(j*n+i-((j+1)*j)/2) = 1.0;
      arma::rowvec e_i(n, arma::fill::zeros);
      e_i(i) = 1.0;
      out += arma::kron(u, arma::kron(e_j, e_i));
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::mat commutation_matrix_cpp(int m, int n) {
  arma::mat K(m * n, m * n, arma::fill::zeros);

  int count = 0;
  for (int k = 0; k < m; ++k) {
    for (int l = 0; l < n; ++l) {
      K(count, m * l + k) = 1;
      count++;
    }
  }

  return K;
}

arma::mat dvechSigmadvechQ_cpp(const arma::mat& Q) {
  int k = Q.n_rows;
  arma::mat L_k = elimination_matrix_cpp(k);
  arma::mat D_k0 = duplication_matrix_cpp(k);
  arma::mat D_k = arma::diagmat(D_k0 * L_k) * D_k0;
  arma::mat K_k = commutation_matrix_cpp(k, k);
  arma::mat I_k = arma::eye(k, k);
  arma::mat result = arma::trans(L_k * (arma::kron(Q, I_k) + arma::kron(I_k, Q) * K_k) * D_k);
  return result;
}

// [[Rcpp::export]]
arma::vec vechdQdP_cpp(const arma::mat& Q) {
  arma::mat A = arma::ones<arma::mat>(Q.n_rows, Q.n_cols);
  A.diag() = Q.diag();
  arma::vec result = vech_cpp(A);
  return result;
}


// [[Rcpp::export]]
arma::vec dldP_cpp(const arma::mat& Q, const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::uvec& new_id) {
  arma::vec tmp_vec(10, arma::fill::zeros);

  for (arma::uword i = 1; i <= arma::max(new_id); ++i) {
    arma::mat tmp_E_ri_riT = E_r_rT.submat(4 * (i - 1), 0, 4 * i - 1, 3);
    tmp_vec += dldvechSigma_cpp(Sigma, tmp_E_ri_riT);
  }

  arma::mat dldvechQ = dvechSigmadvechQ_cpp(Q) * tmp_vec;
  arma::vec dldvechP = dldvechQ % vechdQdP_cpp(Q);  // element-wise multiplication
  return dldvechP;
}

// Function to compute dldP_cen
// [[Rcpp::export]]
arma::vec dldP_cen_cpp(const arma::mat& Q, const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::uvec& new_id, const arma::vec& E_Delta_cen) {
  arma::vec tmp_vec(10, arma::fill::zeros);

  for (arma::uword i = 1; i <= arma::max(new_id); ++i) {
    arma::mat tmp_E_ri_riT = E_r_rT.submat(4 * (i - 1), 0, 4 * i - 1, 3);
    double E_1_Delta = 1.0 - E_Delta_cen(i - 1);
    tmp_vec += dldvechSigma_cpp(Sigma, tmp_E_ri_riT) * E_1_Delta;
  }

  arma::mat dldvechQ = dvechSigmadvechQ_cpp(Q) * tmp_vec;
  arma::vec dldvechP = dldvechQ % vechdQdP_cpp(Q);  // element-wise multiplication
  return dldvechP;
}


// [[Rcpp::export]]
arma::vec dQdvechP_obs_cpp(const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::vec& mu_r, const arma::vec& t, const arma::uvec& new_id) {
  int n = arma::max(new_id);
  arma::mat Qr = trans(chol(Sigma));

  arma::vec first = dldP_cpp(Qr, Sigma, E_r_rT, new_id);

  double sigma_omega = std::sqrt(Sigma(0, 0));
  double mu_omega = mu_r(0);
  double second = 0;

  for (int i = 0; i < n; ++i) {
    double upper = (std::exp(t(i)) - mu_omega) / sigma_omega;
    double lower = -mu_omega / sigma_omega;
    second += - (R::dnorm(upper, 0, 1, 0) * (-upper) - R::dnorm(lower, 0, 1, 0) * (-lower)) /
      (R::pnorm(upper, 0, 1, 1, 0) - R::pnorm(lower, 0, 1, 1, 0));
  }

  arma::vec result = first;
  result(0) += second;

  return result;
}

// [[Rcpp::export]]
arma::vec dQdvechP_cen_cpp(const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::vec& mu_r, const arma::vec& E_g2_ti, const arma::uvec& new_id, const arma::vec& E_Delta_cen) {
  arma::mat Qr = trans(chol(Sigma));

  // Assuming dldP_cen_cpp is already defined and works appropriately
  arma::vec first = dldP_cen_cpp(Qr, Sigma, E_r_rT, new_id, E_Delta_cen);

  double second = -arma::accu(E_g2_ti % (1 - E_Delta_cen));

  arma::vec result = first;
  result(0) += second;

  return result;
}

// [[Rcpp::export]]
arma::vec dQdmur_obs_cpp(const arma::mat& Sigma, const arma::vec& E_r, const arma::vec& mu_r, const arma::vec& t, const arma::uvec& new_id) {
  int n = arma::max(new_id);
  double sigma_omega = std::sqrt(Sigma(0, 0));
  double mu_omega = mu_r(0);
  arma::vec first = arma::zeros<arma::vec>(4);
  double second = 0;

  arma::mat Sigma_inv = inv(Sigma);

  for (int i = 0; i < n; ++i) {
    double upper = (std::exp(t(i)) - mu_omega) / sigma_omega;
    double lower = -mu_omega / sigma_omega;
    first += Sigma_inv * (E_r.subvec(4 * i, 4 * i + 3) - mu_r);
    second += (R::dnorm(upper, 0, 1, false) - R::dnorm(lower, 0, 1, false)) /
      (R::pnorm(upper, 0, 1, true, false) - R::pnorm(lower, 0, 1, true, false)) / sigma_omega;
  }

  arma::vec result = first;
  result(0) += second;
  return result;
}

// [[Rcpp::export]]
arma::vec dQdmur_cen_cpp(const arma::mat& Sigma, const arma::vec& E_r, const arma::vec& mu_r, const arma::vec& E_g1_ti, const arma::uvec& new_id, const arma::vec& E_Delta_cen) {
  int n = arma::max(new_id);
  arma::vec first = arma::zeros<arma::vec>(4);
  double second = arma::accu(E_g1_ti % (1 - E_Delta_cen));

  arma::mat Sigma_inv = inv(Sigma);

  for (int i = 0; i < n; ++i) {
    first += Sigma_inv * (E_r.subvec(4 * i, 4 * i + 3) - mu_r) * (1 - E_Delta_cen(i));
  }

  arma::vec result = first;
  result(0) += second;
  return result;
}


// [[Rcpp::export]]
double Estep_1i_cpp(int i, const arma::vec& E_ri_all, const arma::mat& E_ri_riT_all, const arma::vec& mu_r, const arma::mat& Sigma_r) {
  // Adjusting the index for C++ 0-based indexing
  int idx = 4 * i;

  // Extracting subvector and submatrix corresponding to the ith subject
  arma::vec E_ri = E_ri_all.subvec(idx, idx + 3);
  arma::mat E_ri_riT = E_ri_riT_all.submat(idx, 0, idx + 3, 3);

  // Calculating E_r_rT as specified in the R function
  arma::mat E_r_rT = E_ri_riT - E_ri * trans(mu_r) - mu_r * trans(E_ri) + mu_r * trans(mu_r);

  // Computing and returning the trace of the product of inverse Sigma_r and E_r_rT
  return trace(inv(Sigma_r) * E_r_rT);
}


// [[Rcpp::export]]
double Q_function_cpp(int Nobs, int nobs, const arma::vec& E_r_obs, const arma::mat& E_r_rT_obs, const arma::vec& Estep_2_obs,
                  const arma::vec& yobs, const arma::mat& Xobs, const arma::mat& Xtte_obs, const arma::vec& visittime_obs, const arma::vec& tobs, const arma::uvec& new_id_obs,
                  int Ncen, int ncen, const arma::vec& E_r_cen, const arma::mat& E_r_rT_cen, const arma::vec& Estep_2_cen,
                  const arma::vec& ycen, const arma::mat& Xcen, const arma::mat& Xtte_cen, const arma::vec& visittime_cen, const arma::vec& E_g0_ti, const arma::vec& E_ti, const arma::vec& E_ti_sq, const arma::uvec& new_id_cen,
                  const arma::vec& beta_tte, double sigma_tte_sq, const arma::vec& mu_r, const arma::mat& Sigma_r, const arma::vec& beta_y, double sigma_y_sq,
                  const arma::vec& E_Delta_cen) {

  double first = 0, second = 0, third = 0, fourth = 0;
  double sigma_omega = std::sqrt(Sigma_r(0, 0));
  double mu_omega = mu_r(0);

  for (int i = 0; i < nobs; ++i) {
    first -= 0.5 * (std::log(det(Sigma_r)) + Estep_1i_cpp(i, E_r_obs, E_r_rT_obs, mu_r, Sigma_r));
    double upper = (std::exp(tobs(i)) - mu_omega) / sigma_omega;
    double lower = -mu_omega / sigma_omega;
    second -= std::log(R::pnorm(upper, 0, 1, true, false) - R::pnorm(lower, 0, 1, true, false));
    third -= Estep_2_obs(i) / (2 * sigma_y_sq);
    fourth -= 0.5 * (std::log(sigma_tte_sq) + std::pow(tobs(i) - as_scalar(Xtte_obs.row(i) * beta_tte), 2) / sigma_tte_sq);
  }

  for (int i = 0; i < ncen; ++i) {
    double E_1_Delta = 1 - E_Delta_cen(i);
    first -= E_1_Delta * 0.5 * (std::log(det(Sigma_r)) + Estep_1i_cpp(i, E_r_cen, E_r_rT_cen, mu_r, Sigma_r));
    second -= E_1_Delta * E_g0_ti(i);
    third -= E_1_Delta * Estep_2_cen(i) / (2 * sigma_y_sq);
    double tmp_mu_tte = as_scalar(Xtte_cen.row(i) * beta_tte);
    fourth -= E_1_Delta * 0.5 * (std::log(sigma_tte_sq) + (E_ti_sq(i) - 2 * E_ti(i) * tmp_mu_tte + std::pow(tmp_mu_tte, 2)) / sigma_tte_sq);
  }

  return -(Nobs + Ncen - arma::accu(E_Delta_cen)) * std::log(sigma_y_sq) / 2 + first + second + third + fourth;
}


