#include <RcppArmadillo.h>
#include "E_integral.h"
#include "dmvtnormal.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat Z_ci_cpp(const arma::vec& visittime_vec_i) {
  int n = visittime_vec_i.n_elem;
  arma::mat result(n, 2, fill::ones);  // Initialize a matrix with n rows and 2 columns filled with ones
  result.col(1) = visittime_vec_i;     // Set the second column to visittime_vec_i

  return result;
}

arma::vec E_b_cure_i_cpp(const arma::vec& mu_r_cure, const arma::mat& Sigma_r_cure,
                     const arma::vec& beta_y_cure, double sigma_y_cure_sq,
                     const arma::vec& y_i, const arma::mat& x_i,
                     const arma::vec& visittime_vec_i) {
  int p = visittime_vec_i.n_elem;
  arma::mat Z_ci_value = Z_ci_cpp(visittime_vec_i); // Function must return arma::mat

  arma::mat second = Sigma_r_cure * Z_ci_value.t() *
    inv(diagmat(arma::vec(p, fill::ones) * sigma_y_cure_sq) +
    Z_ci_value * Sigma_r_cure * Z_ci_value.t());
  arma::vec third = y_i - x_i * beta_y_cure - Z_ci_value * mu_r_cure;
  return mu_r_cure + second * third;
}


arma::mat E_b_b_T_mu_cure_i_cpp(const arma::vec& mu_r_cure, const arma::mat& Sigma_r_cure,
                            const arma::vec& beta_y_cure, double sigma_y_cure_sq,
                            const arma::vec& y_i, const arma::mat& x_i,
                            const arma::vec& visittime_vec_i) {
  int p = visittime_vec_i.n_elem;
  arma::mat Z_ci_value = Z_ci_cpp(visittime_vec_i); // Function must return arma::mat
  arma::vec E_b_mu_cure_i = E_b_cure_i_cpp(mu_r_cure, Sigma_r_cure, beta_y_cure,
                                       sigma_y_cure_sq, y_i, x_i, visittime_vec_i) - mu_r_cure;
  arma::mat first = E_b_mu_cure_i * E_b_mu_cure_i.t();
  arma::mat third = Sigma_r_cure * Z_ci_value.t() *
    inv(diagmat(arma::vec(p, fill::ones) * sigma_y_cure_sq) +
    Z_ci_value * Sigma_r_cure * Z_ci_value.t()) *
    Z_ci_value * Sigma_r_cure;
  return first + Sigma_r_cure - third;
}

double E_Delta1_cpp(const arma::vec& mu_r_cure, const arma::mat& Sigma_r_cure,
                const arma::vec& beta_y_cure, double sigma_y_cure_sq,
                const arma::vec& tmp_y_i, const arma::mat& tmp_X_i,
                const arma::vec& tmp_visittime) {
  arma::mat Z_i_cure = Z_ci_cpp(tmp_visittime);
  arma::mat tmp_var = sigma_y_cure_sq * arma::eye(size(Z_i_cure.n_rows, Z_i_cure.n_rows)) +
    Z_i_cure * Sigma_r_cure * Z_i_cure.t();
  arma::vec tmp_mean = tmp_X_i * beta_y_cure + Z_i_cure * mu_r_cure;

  arma::mat density = dmvnorm_arma(tmp_y_i.t(), tmp_mean.t(), tmp_var, false);
  return density[0];
}

// [[Rcpp::export]]
void E_test(const arma::vec& visittime_cen,
            const arma::vec& new_id_cen,
            const arma::mat& Xcen){
  // Rcout << "visittime_cen: " << visittime_cen << std::endl;

}

// [[Rcpp::export]]
List E_loop_obs_cpp(const arma::vec& visittime_obs,
                    const arma::vec& new_id_obs,
                    const arma::mat& Xobs,
                    const arma::vec& yobs,
                    const arma::mat& Xtte_obs,
                    const arma::vec& tobs,
                    const arma::vec& beta_y,
                    const arma::vec& beta_tte,
                    double sigma_tte_sq,
                    double sigma_y_sq,
                    const arma::vec& mu_r,
                    const arma::mat& Sigma_r,
                    int nobs,
                    int Nobs,
                    int Sample_J) {
  arma::vec E_denominator_obs(nobs, fill::zeros);
  arma::vec E_r_obs(4 * nobs, fill::zeros);
  arma::mat E_r_rT_obs(4 * nobs, 4, fill::zeros);
  arma::mat E_r_rT_mu_obs(4 * nobs, 4, fill::zeros);
  arma::vec E_Z_b_obs(Nobs, fill::zeros);
  arma::vec Estep_2_obs(nobs, fill::zeros);

  for (int i = 0; i < nobs; ++i) {
    arma::uvec indices = find(new_id_obs == i + 1);
    arma::vec tmp_visittime = visittime_obs.elem(indices);
    arma::mat tmp_X_i = Xobs.rows(indices);
    arma::vec tmp_y_i = yobs.elem(indices);
    arma::vec tmp_Xtte_i = Xtte_obs.row(i).t();
    arma::vec tmp_m_i = tmp_y_i - tmp_X_i * beta_y;

    List tmp_MC_int = MC_int_xy_all_cpp(tmp_visittime, tmp_y_i,
                                        dot(tmp_Xtte_i, beta_tte), sqrt(sigma_tte_sq),
                                        tmp_X_i * beta_y, sigma_y_sq,
                                        mu_r, Sigma_r, Sample_J, tobs(i), true);

    double tmp_denominator = std::max(as<double>(tmp_MC_int["E_denominator"]), 1e-200);

    arma::vec E_ri_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_ri"]) / tmp_denominator;
    arma::mat E_ri_riT_value_tmp_MC = as<arma::mat>(tmp_MC_int["E_ri_riT"]) / tmp_denominator;
    arma::mat E_r_rT_mu_value_tmp_MC = E_ri_riT_value_tmp_MC - 2 * E_ri_value_tmp_MC * mu_r.t() + mu_r * mu_r.t();

    arma::vec E_Zb_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_Z_i_b_i"]) / tmp_denominator;
    arma::vec E_Zb_sq_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_Z_i_b_i_sq"]) / tmp_denominator;
    double E_Zb_sq_value_tmp_sum_MC = sum(square(tmp_m_i)) - 2 * dot(tmp_m_i, E_Zb_value_tmp_MC) + sum(E_Zb_sq_value_tmp_MC);

    E_denominator_obs(i) = tmp_denominator;
    E_r_obs.subvec(4 * i, 4 * i + 3) = E_ri_value_tmp_MC;
    E_r_rT_obs.submat(4 * i, 0, 4 * i + 3, 3) = E_ri_riT_value_tmp_MC;
    E_r_rT_mu_obs.submat(4 * i, 0, 4 * i + 3, 3) = E_r_rT_mu_value_tmp_MC;
    E_Z_b_obs.elem(indices) = E_Zb_value_tmp_MC;
    Estep_2_obs(i) = E_Zb_sq_value_tmp_sum_MC;
  }

  return List::create(Named("E_denominator_obs") = E_denominator_obs,
                      Named("E_r_obs") = E_r_obs,
                      Named("E_r_rT_obs") = E_r_rT_obs,
                      Named("E_r_rT_mu_obs") = E_r_rT_mu_obs,
                      Named("E_Z_b_obs") = E_Z_b_obs,
                      Named("Estep_2_obs") = Estep_2_obs);
}

// [[Rcpp::export]]
List E_loop_cpp(const arma::vec& visittime_cen,
                             const arma::vec& new_id_cen,
                             const arma::mat& Xcen,
                             const arma::vec& ycen,
                             const arma::mat& Xtte_cen,
                             const arma::vec& tcen,
                             const arma::vec& beta_y,
                             const arma::vec& beta_tte,
                             double sigma_tte_sq,
                             double sigma_y_sq,
                             const arma::vec& mu_r,
                             const arma::mat& Sigma_r,
                             const arma::vec mu_r_cure,
                             const arma::mat Sigma_r_cure,
                             const arma::vec beta_cure,
                             double sigma_y_cure_sq,
                             int ncen,
                             int Ncen,
                             int Sample_J) {
// Rcout << "start: " << ncen << std::endl;
  arma::vec E_r_cen(4 * ncen, fill::zeros);
  arma::mat E_r_rT_cen(4 * ncen, 4, fill::zeros);
  arma::mat E_r_rT_mu_cen(4 * ncen, 4, fill::zeros);
  arma::vec E_Z_b_cen(Ncen, fill::zeros);
  arma::vec Estep_2_cen(ncen, fill::zeros);
  arma::vec E_t_cen(ncen, fill::zeros);
  arma::vec E_t_sq_cen(ncen, fill::zeros);
  arma::vec E_g0_t_cen(ncen, fill::zeros);
  arma::vec E_g1_t_cen(ncen, fill::zeros);
  arma::vec E_g2_t_cen(ncen, fill::zeros);
  arma::vec E_Delta1(ncen, fill::zeros);
  arma::vec E_b_cure(2 * ncen, fill::zeros);
  arma::mat E_b_b_T_mu_cure(2 * ncen, 2, fill::zeros);
  arma::vec E_denominator_cen(ncen, fill::zeros);

  for (int i = 0; i < ncen; ++i) {
    // Rcout << "i: " << i << std::endl;
    arma::uvec indices = find(new_id_cen == i + 1);
    arma::vec tmp_visittime = visittime_cen.elem(indices);
    arma::mat tmp_X_i = Xcen.rows(indices);
    // Rcout << "tmp_X_i: " << tmp_X_i << std::endl;
    arma::vec tmp_y_i = ycen.elem(indices);
    arma::vec tmp_Xtte_i = Xtte_cen.row(i).t();
    // Rcout << "tmp_Xtte_i: " << tmp_Xtte_i << std::endl;

    arma::vec tmp_m_i = tmp_y_i - tmp_X_i * beta_y;

    double mu_tte = dot(tmp_Xtte_i, beta_tte);
    double sd_tte = sqrt(sigma_tte_sq);

    List tmp_MC_int = MC_int_xy_all_cpp(tmp_visittime, tmp_y_i,
                                        mu_tte, sd_tte, tmp_X_i * beta_y, sigma_y_sq,
                                        mu_r, Sigma_r, Sample_J, tcen(i),false);

    double tmp_denominator = std::max(as<double>(tmp_MC_int["E_denominator"]), 1e-200);
// Rcout << "tmp_denominator: " << tmp_denominator << std::endl;
    E_t_cen(i) = as<double>(tmp_MC_int["E_ti"]) / tmp_denominator;
    E_t_sq_cen(i) = as<double>(tmp_MC_int["E_ti_sq"]) / tmp_denominator;
    E_g0_t_cen(i) = as<double>(tmp_MC_int["E_g0_ti"]) / tmp_denominator;
    E_g1_t_cen(i) = as<double>(tmp_MC_int["E_g1_ti"]) / tmp_denominator;
    E_g2_t_cen(i) = as<double>(tmp_MC_int["E_g2_ti"]) / tmp_denominator;
    arma::vec E_ri_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_ri"]) / tmp_denominator;
    arma::mat E_ri_riT_value_tmp_MC = as<arma::mat>(tmp_MC_int["E_ri_riT"]) / tmp_denominator;
// Rcout << "E_ri_riT_value_tmp_MC: " << E_ri_riT_value_tmp_MC << std::endl;
    // Calculate E[(r_i-mu_r)*(r_i-mu_r)^T]
    arma::mat E_r_rT_mu_value_tmp_MC = E_ri_riT_value_tmp_MC - E_ri_value_tmp_MC * mu_r.t() - mu_r * E_ri_value_tmp_MC.t() + mu_r * mu_r.t();
// Rcout << "E_r_rT_mu_value_tmp_MC: " << E_r_rT_mu_value_tmp_MC << std::endl;
    arma::vec E_Zb_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_Z_i_b_i"]) / tmp_denominator;
    arma::vec E_Zb_sq_value_tmp_MC = as<arma::vec>(tmp_MC_int["E_Z_i_b_i_sq"]) / tmp_denominator;
    double E_Zb_sq_value_tmp_sum_MC = sum(square(tmp_m_i)) - 2 * dot(tmp_m_i, E_Zb_value_tmp_MC) + sum(E_Zb_sq_value_tmp_MC);
    E_denominator_cen(i) = tmp_denominator;

    E_r_cen.subvec(4 * i, 4 * i + 3) = E_ri_value_tmp_MC;
// Rcout << "A" << std::endl;
    E_r_rT_cen.submat(4 * i, 0, 4 * i + 3, 3) = E_ri_riT_value_tmp_MC;
// Rcout << "B" << std::endl;
    E_r_rT_mu_cen.submat(4 * i, 0, 4 * i + 3, 3) = E_r_rT_mu_value_tmp_MC;
// Rcout << "C" << std::endl;
// Rcout <<"E_Zb_value_tmp_MC: " << E_Zb_value_tmp_MC << std::endl;
    E_Z_b_cen.elem(indices) = E_Zb_value_tmp_MC;  // Assuming new_id_cen is defined elsewhere
// Rcout << "D" << std::endl;
    Estep_2_cen(i) = E_Zb_sq_value_tmp_sum_MC;

// Rcout <<"start censor" << std::endl;
    E_Delta1(i) = E_Delta1_cpp(mu_r_cure, Sigma_r_cure, beta_cure, sigma_y_cure_sq,
                                 tmp_y_i, tmp_X_i, tmp_visittime);
// Rcout <<"A"<< std::endl;

    E_b_cure.subvec(2 * i, 2 * i + 1) = E_b_cure_i_cpp(mu_r_cure, Sigma_r_cure, beta_cure, sigma_y_cure_sq,
                    tmp_y_i, tmp_X_i, tmp_visittime);
// Rcout<<"B"<< std::endl;
    E_b_b_T_mu_cure.submat(2 * i, 0, 2 * i + 1, 1) = E_b_b_T_mu_cure_i_cpp(mu_r_cure, Sigma_r_cure, beta_cure, sigma_y_cure_sq,
                           tmp_y_i, tmp_X_i, tmp_visittime);
// Rcout << "E_b_b_T_mu_cure: " << E_b_b_T_mu_cure << std::endl;
  }

  return List::create(Named("E_r_cen") = E_r_cen,
                      Named("E_r_rT_cen") = E_r_rT_cen,
                      Named("E_r_rT_mu_cen") = E_r_rT_mu_cen,
                      Named("E_Z_b_cen") = E_Z_b_cen,
                      Named("Estep_2_cen") = Estep_2_cen,
                      Named("E_t_cen") = E_t_cen,
                      Named("E_t_sq_cen") = E_t_sq_cen,
                      Named("E_g0_t_cen") = E_g0_t_cen,
                      Named("E_g1_t_cen") = E_g1_t_cen,
                      Named("E_g2_t_cen") = E_g2_t_cen,
                      Named("E_denominator_cen") = E_denominator_cen,
                      Named("E_Delta1") = E_Delta1,
                      Named("E_b_cure") = E_b_cure,
                      Named("E_b_b_T_mu_cure") = E_b_b_T_mu_cure);
}
