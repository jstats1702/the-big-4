#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

inline double rtruncnorm_one(double mean, 
                             double a, 
                             double b) {
     
     const double lower = std::isfinite(a) ? R::pnorm5(a, mean, 1.0, 1, 0) : 0.0;
     const double upper = std::isfinite(b) ? R::pnorm5(b, mean, 1.0, 1, 0) : 1.0;
     
     double lo = std::min(lower, upper);
     double hi = std::max(lower, upper);
     
     const double eps = 1e-10;
     lo = std::max(lo, eps);
     hi = std::min(hi, 1.0 - eps);
     
     if (hi <= lo) {
          return std::isfinite(a) ? std::max(mean, a) :
          std::isfinite(b) ? std::min(mean, b) : mean;
     }
     
     const double u = R::runif(lo, hi);
     return R::qnorm5(u, mean, 1.0, 1, 0);
}

arma::cube sample_z_cpp(const arma::cube& y,
                        const arma::vec&  mu,
                        const arma::mat&  delta,
                        const double      zeta,
                        arma::cube        z,
                        const int         K,
                        const int         n) {
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk     = del_k[j];
                    const double mean_ijk = mu_k + zeta + d_ik + d_jk;
                    const bool   edge_ijk = (y(i, j, k) > 0.5);
                    
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
                       const double      omega2,
                       const int         K,
                       const int         n) {
     
     const int    M       = n * (n - 1) / 2;
     const double v2_zeta = 1.0 / ( (1.0 / omega2) + static_cast<double>(K * M) );
     
     double sum_res = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk  = del_k[j];
                    
                    sum_res += (z(i, j, k) - mu_k - d_ik - d_jk);
               }
          }
     }
     
     return R::rnorm(v2_zeta * sum_res, std::sqrt(v2_zeta));
}

arma::vec sample_mu_cpp(const arma::cube& z,
                        arma::vec         mu,
                        const arma::mat&  delta,
                        const double      zeta,
                        const double      sigma2,
                        const int         K,
                        const int         n) {

     const int M = n * (n - 1) / 2;
     const double v2_mu = 1.0 / (1.0 / sigma2 + static_cast<double>(M));
     
     for (int k = 0; k < K; ++k) {
          double sum_res = 0.0;
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk = del_k[j];
                    sum_res += (z(i, j, k) - zeta - d_ik - d_jk);
               }
          }
          
          const double m_mu = v2_mu * sum_res;
          mu[k] = R::rnorm(m_mu, std::sqrt(v2_mu));
     }
     
     return mu;
}

arma::mat sample_delta_cpp(const arma::cube& z,
                           const arma::vec&  mu,
                           const double      zeta,
                           const double      tau2,
                           arma::mat         delta,
                           const arma::vec&  vartheta,
                           const int         K,
                           const int         n) {
     
     const double v2_delta = 1.0 / (1.0 / tau2 + static_cast<double>(n - 1));
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n; ++i) {
               double sum_res = 0.0;
               
               for (int j = 0; j < n; ++j) {
                    if (j != i) {
                         const double d_jk = del_k[j];
                         
                         sum_res += (z(i, j, k) - mu_k - zeta - d_jk);
                    }
               }
               
               const double m_delta = v2_delta * ((vartheta[i] / tau2) + sum_res);
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
          const double m_vartheta = v2_vartheta * (arma::accu(delta.row(i)) / tau2);
          
          vartheta[i] = R::rnorm(m_vartheta, std::sqrt(v2_vartheta));
     }
     
     return vartheta;
}

double sample_omega2_cpp(const double zeta,
                         const double a_omega,
                         const double b_omega) {
     
     const double A = a_omega + 0.5;
     const double B = b_omega + 0.5 * zeta * zeta;
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

double sample_sigma2_cpp(const arma::vec& mu,
                         const double     a_sigma,
                         const double     b_sigma,
                         const int        K) {
     
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
                         const double     a_kappa,
                         const double     b_kappa,
                         const int        n) {
     
     const double A = a_kappa + 0.5 * static_cast<double>(n);
     const double B = b_kappa + 0.5 * arma::accu(arma::square(vartheta));
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

// [[Rcpp::export]]
Rcpp::List gibbs_sampler_multilayer(const arma::cube& y, 
                                    const int n_iter, const int n_burn, const int n_thin,
                                    const double a_omega,   const double b_omega,
                                    const double a_sigma,   const double b_sigma, 
                                    const double a_tau,     const double b_tau,
                                    const double a_kappa,   const double b_kappa) {
     // Dimensions
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     // Initialize
     double omega2 = 1.0 / R::rgamma(a_omega, 1.0 / b_omega);
     double sigma2 = 1.0 / R::rgamma(a_sigma, 1.0 / b_sigma);
     double tau2   = 1.0 / R::rgamma(a_tau,   1.0 / b_tau);
     double kappa2 = 1.0 / R::rgamma(a_kappa, 1.0 / b_kappa);
     
     double zeta = R::rnorm(0.0, std::sqrt(omega2));
     
     arma::vec mu(K);
     for (int k = 0; k < K; ++k) mu[k] = R::rnorm(0.0, std::sqrt(sigma2));
     
     arma::vec vartheta(n);
     for (int i = 0; i < n; ++i) vartheta[i] = R::rnorm(0.0, std::sqrt(kappa2));
     
     arma::mat  delta(n, K, arma::fill::zeros);
     arma::cube z(n, n, K, arma::fill::zeros);

     // Storage
     const int n_samples = (n_iter - n_burn) / n_thin;

     arma::cube store_delta    (n_samples, n, K, arma::fill::zeros);
     
     arma::mat  store_mu       (n_samples, K, arma::fill::zeros);
     arma::mat  store_vartheta (n_samples, n, arma::fill::zeros);
     
     arma::vec  store_zeta     (n_samples, arma::fill::zeros);
     arma::vec  store_omega2   (n_samples, arma::fill::zeros);
     arma::vec  store_sigma2   (n_samples, arma::fill::zeros);
     arma::vec  store_tau2     (n_samples, arma::fill::zeros);
     arma::vec  store_kappa2   (n_samples, arma::fill::zeros);
     
     // Sampling
     Rcpp::Rcout << "Initializing Gibbs sampler...\n";
     const int step = n_iter / 20;
     
     for (int iter = 1; iter <= n_iter; ++iter) {
          
          // Update model parameters 
          z        = sample_z_cpp(y, mu, delta, zeta, z, K, n);
          zeta     = sample_zeta_cpp(z, mu, delta, omega2, K, n);
          mu       = sample_mu_cpp(z, mu, delta, zeta, sigma2, K, n);
          delta    = sample_delta_cpp(z, mu, zeta, tau2, delta, vartheta, K, n);
          vartheta = sample_vartheta_cpp(delta, vartheta, kappa2, tau2, n, K);
          omega2   = sample_omega2_cpp(zeta, a_omega, b_omega);
          sigma2   = sample_sigma2_cpp(mu, a_sigma, b_sigma, K);
          tau2     = sample_tau2_cpp(delta, vartheta, a_tau, b_tau, n, K);
          kappa2   = sample_kappa2_cpp(vartheta, a_kappa, b_kappa, n);
          
          // Store
          if ((iter > n_burn) && ((iter - n_burn) % n_thin == 0)) {
               const int pos = (iter - n_burn) / n_thin - 1;
               
               for (int k = 0; k < K; ++k)
                    for (int i = 0; i < n; ++i)
                         store_delta(pos, i, k) = delta(i, k);
               
               store_mu.row(pos)        = mu.t();
               store_vartheta.row(pos)  = vartheta.t();
               
               store_zeta  [pos]        = zeta;
               store_omega2[pos]        = omega2;
               store_sigma2[pos]        = sigma2;
               store_tau2  [pos]        = tau2;
               store_kappa2[pos]        = kappa2;
          }
          
          // Progress 
          if (step > 0 && (iter % step) == 0) {
               int pct = static_cast<int>(std::round(100.0 * iter / n_iter));
               Rcpp::Rcout << "Progress: " << pct << "% completed\n";
          }
     }
     Rcpp::Rcout << "Sampler completed.\n";
     
     return Rcpp::List::create(
          Rcpp::Named("zeta")      = store_zeta,
          Rcpp::Named("mu")        = store_mu,
          Rcpp::Named("delta")     = store_delta,
          Rcpp::Named("vartheta")  = store_vartheta,
          Rcpp::Named("omega2")    = store_omega2,
          Rcpp::Named("sigma2")    = store_sigma2,
          Rcpp::Named("tau2")      = store_tau2,
          Rcpp::Named("kappa2")    = store_kappa2
     );
}

// [[Rcpp::export]]
arma::vec log_likelihood_multilayer_cpp(const arma::cube& y, 
                                        const Rcpp::List& samples) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     arma::vec  zeta  = samples["zeta"];
     arma::mat  mu    = samples["mu"];
     arma::cube delta = samples["delta"];
     
     int n_samples = mu.n_rows;
     arma::vec log_lik_samples(n_samples, arma::fill::zeros);
     
     const double eps = 1e-10;
     
     for (int s = 0; s < n_samples; ++s) {
          double       log_lik = 0.0;
          const double zeta_s  = zeta[s];
          
          for (int k = 0; k < K; ++k) { 
               double        mu_k   = mu(s, k);
               arma::rowvec delta_k = delta.slice(k).row(s);
               
               for (int i = 0; i < n - 1; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                         double eta_ijk = mu_k + zeta_s + delta_k[i] + delta_k[j];
                         double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                         double y_ijk   = y(i, j, k);
                         
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
                               const double      zeta,
                               const double      eps = 1e-10) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     double log_lik = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk = del_k[j];
                    const double eta  = mu_k + zeta + d_ik + d_jk;
                    
                    const double p_ijk = R::pnorm5(eta, 0.0, 1.0, 1, 0);
                    const double y_ijk = y(i, j, k);
                    
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
                                       const double      zeta,
                                       const double      eps = 1e-10) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     const int M = n * (n - 1) / 2;
     arma::vec out(K * M, arma::fill::none);
     
     int idx = 0;
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk    = del_k[j];
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk;
                    
                    const double p_ijk   = R::pnorm5(eta_ijk, 0.0, 1.0, 1, 0);
                    const double y_ijk   = y(i, j, k);
                    
                    out[idx++] = y_ijk * std::log(p_ijk + eps) + (1.0 - y_ijk) * std::log(1.0 - p_ijk + eps);
               }
          }
     }
     
     return out;
}

// [[Rcpp::export]]
arma::cube interaction_prob_cpp(const arma::vec& mu,     
                                const arma::mat& delta,
                                const double     zeta) {
     
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     arma::cube P(n, n, K, arma::fill::zeros);
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const double d_jk    = del_k[j];
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk;
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
arma::cube simulate_multilayer_network_cpp(const arma::vec& mu,
                                           const arma::mat& delta,
                                           const double     zeta) {
     
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     arma::cube P = interaction_prob_cpp(mu, delta, zeta);
     
     arma::cube Y(n, n, K, arma::fill::zeros);
     
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