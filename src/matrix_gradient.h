#ifndef MATRIX_GRADIENT_H
#define MATRIX_GRADIENT_H

#include <RcppArmadillo.h>

arma::mat dSigmadsigmaij_cpp(int i, int j, arma::mat Sigma);
double dldsigmaij_cpp(int i, int j, arma::mat Sigma, arma::mat E_ri_riT);
arma::vec vech_cpp(const arma::mat& mat);
arma::vec dldvechSigma_cpp(arma::mat Sigma, arma::mat E_ri_riT);
arma::mat duplication_matrix_cpp(const int &n);
arma::mat elimination_matrix_cpp(const int &n);
arma::mat commutation_matrix_cpp(int m, int n);
arma::mat dvechSigmadvechQ_cpp(const arma::mat& Q);
arma::vec vechdQdP_cpp(const arma::mat& Q);
arma::vec dldP_cpp(const arma::mat& Q, const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::uvec& new_id);
arma::vec dldP_cen_cpp(const arma::mat& Q, const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::uvec& new_id, const arma::vec& E_Delta_cen);
arma::vec dQdvechP_obs_cpp(const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::vec& mu_r, const arma::vec& t, const arma::uvec& new_id);
arma::vec dQdvechP_cen_cpp(const arma::mat& Sigma, const arma::mat& E_r_rT, const arma::vec& mu_r, const arma::vec& E_g2_ti, const arma::uvec& new_id, const arma::vec& E_Delta_cen);
arma::vec dQdmur_obs_cpp(const arma::mat& Sigma, const arma::mat& E_r, const arma::vec& mu_r, const arma::vec& t, const arma::uvec& new_id);
arma::vec dQdmur_cen_cpp(const arma::mat& Sigma, const arma::mat& E_r, const arma::vec& mu_r, const arma::vec& E_g1_ti, const arma::uvec& new_id, const arma::vec& E_Delta_cen);
double Estep_1i_cpp(int i, const arma::vec& E_ri_all, const arma::mat& E_ri_riT_all, const arma::vec& mu_r, const arma::mat& Sigma_r);
double Q_function_cpp(int Nobs, int nobs, const arma::vec& E_r_obs, const arma::mat& E_r_rT_obs, const arma::vec& Estep_2_obs,
                      const arma::vec& yobs, const arma::mat& Xobs, const arma::mat& Xtte_obs, const arma::vec& visittime_obs, const arma::vec& tobs, const arma::uvec& new_id_obs,
                      int Ncen, int ncen, const arma::vec& E_r_cen, const arma::mat& E_r_rT_cen, const arma::vec& Estep_2_cen,
                      const arma::vec& ycen, const arma::mat& Xcen, const arma::mat& Xtte_cen, const arma::vec& visittime_cen, const arma::vec& E_g0_ti, const arma::vec& E_ti, const arma::vec& E_ti_sq, const arma::uvec& new_id_cen,
                      const arma::vec& beta_tte, double sigma_tte_sq, const arma::vec& mu_r, const arma::mat& Sigma_r, const arma::vec& beta_y, double sigma_y_sq,
                      const arma::vec& E_Delta_cen);

#endif // MATRIX_GRADIENT_H
