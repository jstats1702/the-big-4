#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

inline double rtruncnorm_one(double mean, double a, double b) {
     const double lower = std::isfinite(a) ? R::pnorm5(a, mean, 1.0, 1, 0) : 0.0;
     const double upper = std::isfinite(b) ? R::pnorm5(b, mean, 1.0, 1, 0) : 1.0;
     
     double lo = std::min(lower, upper);
     double hi = std::max(lower, upper);
     
     // margen de seguridad en (0,1)
     const double eps = 1e-10;
     lo = std::max(lo, eps);
     hi = std::min(hi, 1.0 - eps);
     
     if (hi <= lo) {
          // fallback estable cerca del bound válido
          return std::isfinite(a) ? std::max(mean, a) :
          std::isfinite(b) ? std::min(mean, b) : mean;
     }
     
     const double u = R::runif(lo, hi);
     return R::qnorm5(u, mean, 1.0, 1, 0);
}

arma::cube sample_z_cpp(const arma::cube& y,
                        const arma::vec&  mu,
                        const arma::mat&  delta,
                        const arma::mat&  U,
                        const arma::vec&  lambda,
                        const double      zeta,
                        arma::cube        z,
                        const int         K,
                        const int         n) {
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = delta(i, k);
               const arma::rowvec ui   = U.row(i);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = delta(j, k);
                    const arma::rowvec uj   = U.row(j);
                    const double       bi   = -exp_lam_k * arma::norm(ui - uj, 2);
                    
                    const double mean_ijk = mu_k + zeta + d_ik + d_jk + bi;
                    
                    const bool edge_ijk = (y(i, j, k) > 0.5);
                    double z_ijk;
                    if (edge_ijk) {
                         z_ijk = rtruncnorm_one(mean_ijk, 0.0, R_PosInf);
                    } else {
                         z_ijk = rtruncnorm_one(mean_ijk, R_NegInf, 0.0);
                    } 
                    
                    z(i, j, k) = z_ijk;
                    z(j, i, k) = z_ijk;
               } 
          }
     }
     
     return z;
}

double sample_zeta_cpp(const arma::cube& z,
                       const arma::vec&  mu,
                       const arma::mat&  delta,
                       const arma::mat&  U,
                       const arma::vec&  lambda,
                       const double      omega2,
                       const int         K,
                       const int         n) {
     
     const int    M       = n * (n - 1) / 2;
     const double v2_zeta = 1.0 / (1.0 / omega2 + static_cast<double>(K * M));
     
     double sum_res = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const arma::rowvec ui   = U.row(i);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const arma::rowvec uj   = U.row(j);
                    const double       bi   = -exp_lam_k * arma::norm(ui - uj, 2);
                    
                    sum_res += ( z(i, j, k) - mu_k - d_ik - d_jk - bi);
               }
          }
     }
     
     return R::rnorm(v2_zeta * sum_res, std::sqrt(v2_zeta));
}

arma::vec sample_mu_cpp(const arma::cube& z,
                        arma::vec         mu,
                        const arma::mat&  delta,
                        const arma::mat&  U,
                        const arma::vec&  lambda,
                        const double      zeta,
                        const double      sigma2,
                        const int         K,
                        const int         n) {
     
     const int    M     = n * (n - 1) / 2;
     const double v2_mu = 1.0 / (1.0 / sigma2 + double(M));
     
     for (int k = 0; k < K; ++k) {
          double          sum_res   = 0.0;
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const arma::rowvec ui   = U.row(i);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const arma::rowvec uj   = U.row(j);
                    const double       bi   = -exp_lam_k * arma::norm(ui - uj, 2);
                    
                    sum_res += (z(i, j, k) - zeta - d_ik - d_jk - bi);
               } 
          }
          
          const double m_mu = v2_mu * sum_res;
          mu[k] = R::rnorm(m_mu, std::sqrt(v2_mu));
     } 
     
     return mu;
}

arma::mat sample_delta_cpp(const arma::cube& z,
                           const arma::vec&  mu,
                           const double      tau2,
                           arma::mat         delta,
                           const arma::vec&  vartheta,
                           const arma::mat&  U,
                           const arma::vec&  lambda,
                           const double      zeta,
                           const int         K,
                           const int         n) {
     
     const double v2_delta = 1.0 / (1.0 / tau2 + double(n - 1));
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n; ++i) {
               double             sum_res = 0.0;
               const arma::rowvec ui      = U.row(i);
               
               for (int j = 0; j < n; ++j) {
                    if (j != i) {
                         const double       d_jk = del_k[j];
                         const arma::rowvec uj   = U.row(j);
                         const double       bi   = -exp_lam_k * arma::norm(ui - uj, 2);
                         
                         sum_res += (z(i, j, k) - mu_k - zeta - d_jk - bi);
                    }
               }  
               
               const double m_delta = v2_delta * (vartheta[i] / tau2 + sum_res);
               delta(i, k) = R::rnorm(m_delta, std::sqrt(v2_delta));
          }  
     }
     return delta;
}  

arma::vec sample_vartheta_cpp(const arma::mat& delta,
                              arma::vec        vartheta,
                              const double     kappa2,
                              const double     tau2,
                              const int        n,
                              const int        K) {
     
     const double v2_vartheta = 1.0 / (1.0 / kappa2 + static_cast<double>(K) / tau2);
     
     for (int i = 0; i < n; ++i) {
          const double m_vartheta  = v2_vartheta * (arma::accu(delta.row(i)) / tau2);
          
          vartheta[i] = R::rnorm(m_vartheta, std::sqrt(v2_vartheta));
     }  
     
     return vartheta;
}  

inline double logpost_lambda(const double lam,
                             const double S1,
                             const double S2,
                             const double upsilon2) {
     
     const double el = std::exp(lam);
     double out = - el * S1 - 0.5 * el * el * S2;
     out += -0.5 * (lam * lam) / upsilon2;
     return out;
}

Rcpp::List sample_lambda_cpp(const arma::cube& z,
                             const arma::vec&  mu,
                             const arma::mat&  delta,
                             const arma::mat&  U,
                             arma::vec         lambda,
                             const double      zeta,
                             const double      upsilon2,
                             arma::vec         log_step,
                             const int         K,
                             const int         n,
                             const int         iter,
                             const int         burn_in,
                             const double      target,
                             const double      eta0) {
     
     int n_acc = 0;
     
     // Polynomially decaying learning rate (Robbins–Monro)
     const double gamma0  = eta0;
     const double gamma_t = gamma0 / std::sqrt(1.0 + static_cast<double>(iter));
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k = mu[k];
          
          // S1_k = sum_{i<j} r_ij,k * ||u_i - u_j||
          // S2   = sum_{i<j} ||u_i - u_j||^2
          double S1 = 0.0;
          double S2 = 0.0;
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = delta(i, k);
               const arma::rowvec u_i  = U.row(i);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = delta(j, k);
                    const arma::rowvec u_j  = U.row(j);
                    const double       d_ij = arma::norm(u_i - u_j, 2);
                    const double       r_ij = z(i, j, k) - mu_k - zeta - d_ik - d_jk;
                    
                    S1 += r_ij * d_ij;
                    S2 += d_ij * d_ij;
               }
          }
          
          const double lam_curr = lambda[k];
          
          // Per-layer log-step
          double log_step_k = log_step[k];
          double step_k     = std::exp(log_step_k);
          
          // Gaussian RW on lambda_k
          const double lam_prop   = lam_curr + step_k * R::rnorm(0.0, 1.0);
          
          const double lpost_curr = logpost_lambda(lam_curr, S1, S2, upsilon2);
          const double lpost_prop = logpost_lambda(lam_prop, S1, S2, upsilon2);
          const double log_r      = lpost_prop - lpost_curr;
          const double log_u      = std::log(R::runif(0.0, 1.0));
          
          int accepted = 0;
          if (log_u < log_r) {
               lambda[k] = lam_prop;
               ++n_acc;
               accepted = 1;
          }
          
          // Per-layer adaptation of log_step_k (Robbins–Monro)
          if (burn_in > 0 && iter <= burn_in) {
               log_step_k += gamma_t * (static_cast<double>(accepted) - target);
          }
          
          // Clamp per-layer log-step (lambda can tolerate a bit wider range)
          log_step_k  = std::min(std::max(log_step_k, std::log(1e-6)), std::log(5.0));
          log_step[k] = log_step_k;
     }
     
     const double acc_rate = static_cast<double>(n_acc) / static_cast<double>(K);
     
     return Rcpp::List::create(
          Rcpp::Named("lambda")     = lambda,
          Rcpp::Named("log_step")   = log_step,
          Rcpp::Named("acc_rate")   = acc_rate,
          Rcpp::Named("n_accepted") = n_acc
     );
}

inline double logpost_ui(const arma::rowvec& ui,
                         const int           i,
                         const arma::cube&   z,
                         const arma::vec&    mu,
                         const arma::mat&    delta,
                         const arma::mat&    U,
                         const arma::vec&    lambda,
                         const double        zeta,
                         const int           K,
                         const int           n) {
     
     double S1_total = 0.0;  // sum_k e^{λ_k}  * sum_{j≠i} r_{ij,k} d_ij
     double S2_total = 0.0;  // sum_k e^{2λ_k} * sum_{j≠i} d_ij^2
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k   = mu[k];
          const double    el     = std::exp(lambda[k]);
          const double    el2    = el * el;
          
          double S1_k = 0.0;
          double S2_k = 0.0;
          
          for (int j = 0; j < n; ++j) {
               if (j != i) {
                    const arma::rowvec uj   = U.row(j);
                    const double       rij  = z(i, j, k) - mu_k - zeta - delta(i, k) - delta(j, k);
                    const double       dij  = arma::norm(ui - uj, 2);
                    
                    S1_k += rij * dij;
                    S2_k += dij * dij;
               }
          }
          
          S1_total += el  * S1_k;
          S2_total += el2 * S2_k;
     }
     
     // Log-posterior kernel (up to additive constant)
     return - S1_total - 0.5 * S2_total - 0.5 * arma::dot(ui, ui);
}

Rcpp::List sample_U_cpp(const arma::cube& z,
                        const arma::vec&  mu,
                        const arma::mat&  delta,
                        arma::mat         U,
                        const arma::vec&  lambda,
                        const double      zeta,
                        arma::vec         log_step,
                        const int         K,
                        const int         n,
                        const int         d,
                        const int         iter,
                        const int         burn_in,
                        const double      target,
                        const double      eta0) {
     
     int n_acc = 0;
     
     // Polynomially decaying learning rate (Robbins–Monro)
     const double gamma0  = eta0;
     const double gamma_t = gamma0 / std::sqrt(1.0 + static_cast<double>(iter));
     
     for (int i = 0; i < n; ++i) {
          const arma::rowvec ui_curr = U.row(i);
          
          // Gaussian noise
          arma::rowvec eps(d);
          for (int j = 0; j < d; ++j)
               eps[j] = R::rnorm(0.0, 1.0);
          
          // Per-node log-step and dimension-invariant scaling
          double log_step_i = log_step[i];
          double step_i     = std::exp(log_step_i);
          double scale_i    = step_i / std::sqrt(static_cast<double>(d));
          
          const arma::rowvec ui_prop = ui_curr + scale_i * eps;
          
          const double lp_curr = logpost_ui(ui_curr, i, z, mu, delta, U, lambda, zeta, K, n);
          const double lp_prop = logpost_ui(ui_prop, i, z, mu, delta, U, lambda, zeta, K, n);
          const double log_r   = lp_prop - lp_curr;
          const double log_u   = std::log(R::runif(0.0, 1.0));
          
          int accepted = 0;
          if (log_u < log_r) {
               U.row(i) = ui_prop;
               ++n_acc;
               accepted = 1;
          }
          
          // Per-node adaptation of log_step_i (Robbins–Monro)
          if (burn_in > 0 && iter <= burn_in) {
               // use the Bernoulli(accepted) as a local "acc_rate_i"
               log_step_i += gamma_t * (static_cast<double>(accepted) - target);
          }
          
          // Clamp per-node log-step (U often needs smaller steps)
          log_step_i  = std::min(std::max(log_step_i, std::log(1e-4)), std::log(1.0));
          log_step[i] = log_step_i;
     }
     
     const double acc_rate = static_cast<double>(n_acc) / static_cast<double>(n);
     
     return Rcpp::List::create(
          Rcpp::Named("U")          = U,
          Rcpp::Named("log_step")   = log_step,
          Rcpp::Named("acc_rate")   = acc_rate,
          Rcpp::Named("n_accepted") = n_acc
     );
}

double sample_omega2_cpp(const double zeta,
                         const double a_omega,
                         const double b_omega) {
     
     const double A = a_omega + 0.5;
     const double B = b_omega + 0.5 * zeta * zeta;
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

double sample_sigma2_cpp(const arma::vec& mu,
                         const double a_sigma,
                         const double b_sigma,
                         const int    K) {
     
     const double A = a_sigma + 0.5 * static_cast<double>(K);
     const double B = b_sigma + 0.5 * arma::accu(arma::square(mu));
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

double sample_tau2_cpp(const arma::mat& delta, 
                       const arma::vec& vartheta, 
                       const double     a_tau,
                       const double     b_tau,
                       const int        n,
                       const int        K) {
     
     double ss = 0.0;
     
     for (int i = 0; i < n; ++i) {
          const double vi = vartheta[i];
          
          for (int k = 0; k < K; ++k) {
               const double diff = delta(i, k) - vi;
               
               ss += diff * diff;
          }
     }
     
     const double A = a_tau + 0.5 * static_cast<double>(n * K);
     const double B = b_tau + 0.5 * ss;
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

double sample_kappa2_cpp(const arma::vec& vartheta,
                         const double a_kappa,
                         const double b_kappa,
                         const int    n) {
     
     const double A = a_kappa + 0.5 * static_cast<double>(n);
     const double B = b_kappa + 0.5 * arma::accu(arma::square(vartheta));
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

double sample_upsilon2_cpp(const arma::vec& lambda,
                           const double     a_u,
                           const double     b_u,
                           const int        K) {
     
     const double A = a_u + 0.5 * static_cast<double>(K);
     const double B = b_u + 0.5 * arma::dot(lambda, lambda);
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

// [[Rcpp::export]]
Rcpp::List gibbs_sampler_multilayer(const arma::cube& y,
                                    const int         d,
                                    const int n_iter, const int n_burn, const int n_thin,
                                    const double a_omega,    const double b_omega,
                                    const double a_sigma,    const double b_sigma,
                                    const double a_tau,      const double b_tau,
                                    const double a_kappa,    const double b_kappa,
                                    const double a_upsilon,  const double b_upsilon) {
     
     // Dimensions
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     // Initialize
     double omega2    = 1.0 / R::rgamma(a_omega,    1.0 / b_omega);
     double sigma2    = 1.0 / R::rgamma(a_sigma,    1.0 / b_sigma);
     double tau2      = 1.0 / R::rgamma(a_tau,      1.0 / b_tau);
     double kappa2    = 1.0 / R::rgamma(a_kappa,    1.0 / b_kappa);
     double upsilon2  = 1.0 / R::rgamma(a_upsilon,  1.0 / b_upsilon);
     
     double zeta = R::rnorm(0.0, std::sqrt(omega2));

     arma::vec mu(K);
     for (int k = 0; k < K; ++k) mu[k] = R::rnorm(0.0, std::sqrt(sigma2));
     
     arma::vec vartheta(n);
     for (int i = 0; i < n; ++i) vartheta[i] = R::rnorm(0.0, std::sqrt(kappa2));
     
     arma::mat  delta(n, K, arma::fill::zeros);
     arma::cube z(n, n, K, arma::fill::zeros);
     
     arma::vec lambda(K);
     for (int k = 0; k < K; ++k) lambda[k] = R::rnorm(0.0, std::sqrt(upsilon2));
     
     arma::mat U(n, d);
     for (int i = 0; i < n; ++i)
          for (int t = 0; t < d; ++t)
               U(i, t) = R::rnorm(0.0, 1.0);
     
     // RW–MH tuning: per-parameter log-step vectors
     arma::vec log_step_U(n);
     log_step_U.fill(std::log(0.1));
     
     arma::vec log_step_lambda(K);
     log_step_lambda.fill(std::log(0.1));
     
     // typical targets
     const double target_U       = (d > 1 ? 0.234 : 0.44);
     const double target_lambda  = 0.44;
     
     // initial learning rate
     const double eta0 = 0.05;
     
     // Cumulative acceptance tracking
     long long cum_acc_U = 0, cum_trials_U = 0;
     long long cum_acc_L = 0, cum_trials_L = 0;
     
     // Storage
     const int n_samples = (n_iter - n_burn) / n_thin;
     
     arma::cube  store_delta     (n_samples, n, K, arma::fill::zeros);
     arma::cube  store_U         (n_samples, n, d, arma::fill::zeros);
     
     arma::mat   store_mu        (n_samples, K, arma::fill::zeros);
     arma::mat   store_vartheta  (n_samples, n, arma::fill::zeros);
     arma::mat   store_lambda    (n_samples, K, arma::fill::zeros);
     
     arma::vec   store_zeta      (n_samples, arma::fill::zeros);
     arma::vec   store_omega2    (n_samples, arma::fill::zeros);
     arma::vec   store_sigma2    (n_samples, arma::fill::zeros);
     arma::vec   store_tau2      (n_samples, arma::fill::zeros);
     arma::vec   store_kappa2    (n_samples, arma::fill::zeros);
     arma::vec   store_upsilon2  (n_samples, arma::fill::zeros);
     
     // Sampling
     Rcpp::Rcout << "Initializing Gibbs sampler...\n";
     const int step = std::max(1, n_iter / 20);
     
     for (int iter = 1; iter <= n_iter; ++iter) {

          // RW–MH blocks with shared, adaptive step sizes
          {
               // U update
               Rcpp::List resU = sample_U_cpp(z, mu, delta, U, lambda, zeta, log_step_U, K, n, d, iter, n_burn, target_U, eta0);
               U            = Rcpp::as<arma::mat>(resU["U"]);
               log_step_U   = Rcpp::as<arma::vec>(resU["log_step"]);
               cum_acc_U   += Rcpp::as<int>(resU["n_accepted"]);
               cum_trials_U += n;
          } 
          {
               // lambda update
               Rcpp::List resL = sample_lambda_cpp(z, mu, delta, U, lambda, zeta, upsilon2, log_step_lambda, K, n, iter, n_burn, target_lambda, eta0);
               lambda          = Rcpp::as<arma::vec>(resL["lambda"]);
               log_step_lambda = Rcpp::as<arma::vec>(resL["log_step"]);
               cum_acc_L      += Rcpp::as<int>(resL["n_accepted"]);
               cum_trials_L   += K;
          }
          
          // Update model parameters
          z         = sample_z_cpp(y, mu, delta, U, lambda, zeta, z, K, n);
          zeta      = sample_zeta_cpp(z, mu, delta, U, lambda, omega2, K, n);
          mu        = sample_mu_cpp(z, mu, delta, U, lambda, zeta, sigma2, K, n);
          delta     = sample_delta_cpp(z, mu, tau2, delta, vartheta, U, lambda, zeta, K, n);
          vartheta  = sample_vartheta_cpp(delta, vartheta, kappa2, tau2, n, K);
          omega2    = sample_omega2_cpp(zeta, a_omega, b_omega);
          sigma2    = sample_sigma2_cpp(mu, a_sigma, b_sigma, K);
          tau2      = sample_tau2_cpp(delta, vartheta, a_tau, b_tau, n, K);
          kappa2    = sample_kappa2_cpp(vartheta, a_kappa, b_kappa, n);
          upsilon2  = sample_upsilon2_cpp(lambda, a_upsilon, b_upsilon, K);
          
          // Store
          if ((iter > n_burn) && ((iter - n_burn) % n_thin == 0)) {
               const int pos = (iter - n_burn) / n_thin - 1;
               
               for (int k = 0; k < K; ++k)
                    for (int i = 0; i < n; ++i)
                         store_delta(pos, i, k) = delta(i, k);
               
               for (int k = 0; k < K; ++k)
                    store_lambda(pos, k) = lambda[k];
               
               for (int i = 0; i < n; ++i)
                    for (int t = 0; t < d; ++t)
                         store_U(pos, i, t) = U(i, t);
               
               store_mu.row(pos)        = mu.t();
               store_vartheta.row(pos)  = vartheta.t();
               
               store_zeta     [pos]     = zeta;
               store_omega2   [pos]     = omega2;
               store_sigma2   [pos]     = sigma2;
               store_tau2     [pos]     = tau2;
               store_kappa2   [pos]     = kappa2;
               store_upsilon2 [pos]     = upsilon2;
          }
          
          // Progress: print cumulative acceptance so far
          if ((iter % step) == 0) {
               int pct = static_cast<int>(std::round(100.0 * iter / n_iter));
               double accU_cum = (cum_trials_U > 0) ? (static_cast<double>(cum_acc_U) / cum_trials_U) : 0.0;
               double accL_cum = (cum_trials_L > 0) ? (static_cast<double>(cum_acc_L) / cum_trials_L) : 0.0;
               int accU_pct = static_cast<int>(std::round(100.0 * accU_cum));
               int accL_pct = static_cast<int>(std::round(100.0 * accL_cum));
               Rcpp::Rcout << "Progress: " << pct << "% completed"
                           << " | cumulative acc_U ≈ " << accU_pct << "%, "
                           << "acc_lambda ≈ " << accL_pct << "%\n";
          } 
     } 
     Rcpp::Rcout << "Sampler completed.\n";
     
     return Rcpp::List::create(
          Rcpp::Named("zeta")       = store_zeta,
          Rcpp::Named("mu")         = store_mu,
          Rcpp::Named("delta")      = store_delta,
          Rcpp::Named("vartheta")   = store_vartheta,
          Rcpp::Named("lambda")     = store_lambda,
          Rcpp::Named("U")          = store_U,
          Rcpp::Named("omega2")     = store_omega2,
          Rcpp::Named("sigma2")     = store_sigma2,
          Rcpp::Named("tau2")       = store_tau2,
          Rcpp::Named("kappa2")     = store_kappa2,
          Rcpp::Named("upsilon2")   = store_upsilon2
     );
}

// [[Rcpp::export]]
arma::vec log_likelihood_multilayer_cpp(const arma::cube& y,
                                        const Rcpp::List& samples) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     arma::vec  zeta   = samples["zeta"];
     arma::mat  mu     = samples["mu"];
     arma::cube delta  = samples["delta"];
     arma::mat  lambda = samples["lambda"];
     arma::cube Ucube  = samples["U"];
     
     const int d         = Ucube.n_slices;
     const int n_samples = mu.n_rows;
     arma::vec log_lik_samples(n_samples, arma::fill::zeros);
     
     const double eps = 1e-10;
     
     for (int s = 0; s < n_samples; ++s) {
          double log_lik = 0.0;
          const double zeta_s  = zeta[s];
          
          arma::mat U_s(n, d);
          for (int t = 0; t < d; ++t) {
               U_s.col(t) = Ucube.slice(t).row(s).t();
          }
          
          for (int k = 0; k < K; ++k) {
               const double       mu_k      = mu(s, k);
               const arma::rowvec del_k     = delta.slice(k).row(s);
               const double       exp_lam_k = std::exp(lambda(s, k));
               
               for (int i = 0; i < n - 1; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                         const double    bi   = - exp_lam_k * arma::norm(U_s.row(i) - U_s.row(j), 2);
                         
                         const double eta_ijk = mu_k + zeta_s + del_k[i] + del_k[j] + bi;
                         const double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                         const double y_ijk   = y(i, j, k);
                         
                         log_lik += y_ijk * std::log(p_ijk + eps) + (1.0 - y_ijk) * std::log(1.0 - p_ijk + eps);
                    } 
               }
          }
          log_lik_samples[s] = log_lik;
     }
     
     return log_lik_samples;
} 

// [[Rcpp::export]]
double log_likelihood_iter_cpp(const arma::cube& y,
                               const arma::vec&  mu,
                               const arma::mat&  delta,
                               const arma::mat&  S,
                               const arma::vec&  lambda,
                               const double      zeta,
                               const double      eps = 1e-10) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     double log_lik = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    // Distance via Gram: ||u_i - u_j|| = sqrt(Sii + Sjj - 2 Sij)
                    const double Sii  = S(i, i);
                    const double Sjj  = S(j, j);
                    const double Sij  = S(i, j);
                    const double dij2 = std::max(0.0, Sii + Sjj - 2.0 * Sij);
                    const double dij  = std::sqrt(dij2);
                    
                    const double    d_jk = del_k[j];
                    const double    bi   = - exp_lam_k * dij;
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + bi;
                    const double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                    const double y_ijk   = y(i, j, k);
                    
                    log_lik += y_ijk * std::log(p_ijk + eps) + (1.0 - y_ijk) * std::log(1.0 - p_ijk + eps);
               }
          }
     }
     
     return log_lik;
}

// [[Rcpp::export]]
arma::vec log_likelihood_pointwise_cpp(const arma::cube& y,
                                       const arma::vec&  mu,
                                       const arma::mat&  delta,
                                       const arma::mat&  U, 
                                       const arma::vec&  lambda,
                                       const double      zeta,
                                       const double      eps = 1e-10) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     const int M = n * (n - 1) / 2;
     arma::vec out(K * M, arma::fill::none);
     
     int idx = 0;
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const arma::rowvec ui   = U.row(i);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const arma::rowvec uj   = U.row(j);
                    const double       bi   = - exp_lam_k * arma::norm(ui - uj, 2);
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + bi;
                    const double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                    const double y_ijk   = y(i, j, k);
                    
                    out[idx++] = y_ijk * std::log(p_ijk + eps) + (1.0 - y_ijk) * std::log(1.0 - p_ijk + eps);
               }
          }
     }
     
     return out;
}

// [[Rcpp::export]]
arma::cube interaction_prob_cpp(const arma::vec&  mu,
                                const arma::mat&  delta,
                                const arma::mat&  U,
                                const arma::vec&  lambda,
                                const double      zeta) {
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     // Precompute pairwise distances D (shared across layers)
     arma::mat D(n, n, arma::fill::zeros);
     for (int i = 0; i < n - 1; ++i) {
          const arma::rowvec ui = U.row(i);
          for (int j = i + 1; j < n; ++j) {
               const arma::rowvec uj = U.row(j);
               const double dij = arma::norm(ui - uj, 2);
               D(i, j) = dij;
               D(j, i) = dij;
          } 
     } 
     
     arma::cube P(n, n, K, arma::fill::zeros);
     
     // Build probabilities (upper triangle; mirror to ensure symmetry)
     for (int k = 0; k < K; ++k) {
          const double    mu_k      = mu[k];
          const arma::vec del_k     = delta.col(k);
          const double    exp_lam_k = std::exp(lambda[k]);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk = del_k[j];
                    const double bi   = - exp_lam_k * D(i, j);
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + bi;
                    const double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                    
                    P(i, j, k) = p_ijk;
                    P(j, i, k) = p_ijk;
               } 
          }
          
          for (int i = 0; i < n; ++i) P(i, i, k) = 0.0;
     } 
     
     return P;
} 

// [[Rcpp::export]]
arma::cube simulate_multilayer_network_cpp(const arma::vec&  mu,
                                           const arma::mat&  delta,
                                           const arma::mat&  U,
                                           const arma::vec&  lambda,
                                           const double      zeta) {
     
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     arma::cube Y(n, n, K, arma::fill::zeros);
     
     arma::cube P = interaction_prob_cpp(mu, delta, U, lambda, zeta);
     
     for (int k = 0; k < K; ++k) {
          for (int i = 0; i < n - 1; ++i) {
               for (int j = i + 1; j < n; ++j) {
                    const double p_ijk = P(i, j, k);
                    const double y_ijk = (R::unif_rand() < p_ijk) ? 1.0 : 0.0;
                    
                    Y(i, j, k) = y_ijk;
                    Y(j, i, k) = y_ijk;
               }
               Y(i, i, k) = 0.0;
          } 
          Y(n - 1, n - 1, k) = 0.0;
     } 
     
     return Y;
} 