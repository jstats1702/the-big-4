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
                        const arma::cube& Theta,
                        const arma::umat& xi,
                        const double      zeta,
                        arma::cube        z,
                        const int         K,
                        const int         n) {
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k   = mu[k];
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = delta(i, k);
               const unsigned int ci   = xi(i, k);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = delta(j, k);
                    const unsigned int cj   = xi(j, k);
                    const double       tabk = Theta(ci, cj, k);
                    
                    const double mean_ijk = zeta + mu_k + d_ik + d_jk + tabk;
                    
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
                       const arma::cube& Theta,
                       const arma::umat& xi,
                       const double      omega2,
                       const int         K,
                       const int         n) {
     
     const int    M       = n * (n - 1) / 2;
     const double v2_zeta = 1.0 / (1.0 / omega2 + static_cast<double>(K * M));
     
     double sum_res = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k   = mu[k];
          const arma::vec del_k  = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const unsigned int ci   = xi(i, k);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const unsigned int cj   = xi(j, k); 
                    const double       tabk = Theta(ci, cj, k);
                    
                    sum_res += (z(i, j, k) - mu_k - d_ik - d_jk - tabk);
               }
          }
     }
     
     return R::rnorm(v2_zeta * sum_res, std::sqrt(v2_zeta));
}

arma::vec sample_mu_cpp(const arma::cube& z,
                        arma::vec         mu,
                        const arma::mat&  delta, 
                        const arma::cube& Theta,
                        const arma::umat& xi,
                        const double      zeta,
                        const double      sigma2,
                        const int         K,
                        const int         n) {
     
     const int    M     = n * (n - 1) / 2;
     const double v2_mu = 1.0 / (1.0 / sigma2 + double(M));
     
     for (int k = 0; k < K; ++k) {
          double          sum_res = 0.0;
          const arma::vec del_k   = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const unsigned int ci   = xi(i, k);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const unsigned int cj   = xi(j, k); 
                    const double       tabk = Theta(ci, cj, k);
                    
                    sum_res += (z(i, j, k) - zeta - d_ik - d_jk - tabk);
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
                           const arma::cube& Theta,
                           const arma::umat& xi,
                           const double      zeta,
                           const int         K,
                           const int         n) {
     
     const double v2_delta = 1.0 / (1.0 / tau2 + double(n - 1));
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k   = mu[k];
          const arma::vec del_k  = delta.col(k);
          
          for (int i = 0; i < n; ++i) {
               double             sum_res = 0.0;
               const unsigned int ci      = xi(i, k);
               
               for (int j = 0; j < n; ++j) {
                    if (j != i) {
                         const double       d_jk = del_k[j];
                         const unsigned int cj   = xi(j, k);
                         const double       tabk = Theta(ci, cj, k);
                         
                         sum_res += (z(i, j, k) - mu_k - zeta - d_jk - tabk);
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

arma::cube sample_Theta_cpp(const arma::cube& z,
                            const arma::vec&  mu,
                            const arma::mat&  delta, 
                            const arma::umat& xi,
                            arma::cube        Theta, 
                            const double      rho2,
                            const double      zeta,
                            const int         K,
                            const int         n,
                            const int         C) {
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          // Accumulate sufficient stats per (a,b) with a<=b
          arma::mat sum_res(C, C, arma::fill::zeros);
          arma::mat counts (C, C, arma::fill::zeros);
          
          for (int i = 0; i < n - 1; ++i) {
               const unsigned int ai   = xi(i, k);
               const double       d_ik = del_k[i];
               
               for (int j = i + 1; j < n; ++j) {
                    const unsigned int bj    = xi(j, k);
                    const double       d_jk  = del_k[j];
                    const double       r_ijk = z(i, j, k) - zeta - mu_k - d_ik - d_jk;
                    
                    unsigned int a = ai;
                    unsigned int b = bj;
                    if (a > b) std::swap(a, b);
                    
                    sum_res(a, b) += r_ijk;
                    counts (a, b) += 1.0;
               }
          }
         
          for (int a = 0; a < C; ++a) {
               for (int b = a; b < C; ++b) {
                    const double N_ab = counts(a, b);  // number of dyads in block (a,b) at layer k
                    const double v2   = 1.0 / (1.0 / rho2 + N_ab);
                    const double m    = v2 * sum_res(a, b);
                    
                    const double th = R::rnorm(m, std::sqrt(v2));
                    Theta(a, b, k)  = th;
                    Theta(b, a, k)  = th;
               }
          }
     }
    
     return Theta;
}

arma::umat sample_xi_cpp(const arma::cube& z,
                         const arma::vec&  mu, 
                         const arma::mat&  delta,
                         arma::umat        xi,
                         const arma::cube& Theta, 
                         const double      zeta,
                         const double      alpha,
                         const int         K,
                         const int         n,
                         const int         C) {
     
     // Scratch objects
     arma::vec logp(C);
     arma::vec probs(C);
     
     for (int k = 0; k < K; ++k) {
          // Current cluster counts for layer k
          arma::ivec counts(C, arma::fill::zeros);
          for (int i = 0; i < n; ++i) counts[ xi(i, k) ]++;
          
          const double    mu_k   = mu[k];
          
          for (int i = 0; i < n; ++i) {
               // Remove i from its current cluster (layer k)
               counts[ xi(i, k) ]--;
               
               // For each candidate cluster c, compute log posterior up to a constant (layer k only)
               for (int c = 0; c < C; ++c) {
                    double ll = std::log( static_cast<double>(counts[c]) + alpha / static_cast<double>(C) ); 
                    
                    const double base_i_k = zeta + mu_k + delta(i, k);
                    for (int j = 0; j < n; ++j) {
                         if (j == i) continue;
                         
                         const double    base_j_k = delta(j, k);
                         const unsigned  cj       = xi(j, k);
                         const double    theta    = Theta(c, cj, k);
                         
                         const double mean_ijk = base_i_k + base_j_k + theta;
                         const double zij      = (i < j) ? z(i, j, k) : z(j, i, k);
                         
                         ll += -0.5 * (zij - mean_ijk) * (zij - mean_ijk);
                    }
                    logp[c] = ll;
               }
               
               // Normalize with log-sum-exp
               const double mlog = logp.max();
               double sumexp = 0.0;
               for (int c = 0; c < C; ++c) {
                    probs[c] = std::exp(logp[c] - mlog);
                    sumexp  += probs[c];
               }
               for (int c = 0; c < C; ++c) probs[c] /= sumexp;
               
               // Sample new cluster for i at layer k
               const double u = R::runif(0.0, 1.0);
               double cum = 0.0;
               unsigned int new_c = 0u;
               for (int c = 0; c < C; ++c) {
                    cum += probs[c];
                    if (u <= cum) { 
                         new_c = static_cast<unsigned int>(c); 
                         break; 
                    }
               }
               
               xi(i, k) = new_c;
               counts[new_c]++;  // put i back
          }
     }
    
     return xi;
}

arma::mat sample_omega_cpp(const arma::umat& xi,
                           const double      alpha,
                           const int         K,
                           const int         n,
                           const int         C) {
     
     arma::mat omega(K, C, arma::fill::zeros);
     
     const double base = alpha / static_cast<double>(C);
     arma::ivec counts(C, arma::fill::zeros);
     arma::vec  shape(C);
     arma::vec  g(C);
     
     for (int k = 0; k < K; ++k) {
          // Counts for layer k
          counts.zeros();
          for (int i = 0; i < n; ++i) counts[ xi(i, k) ]++;
          
          // Dirichlet parameters for layer k
          for (int c = 0; c < C; ++c) shape[c] = base + static_cast<double>(counts[c]);
          
          // Sample via normalized Gammas
          double sumg = 0.0;
          for (int c = 0; c < C; ++c) {
               g[c] = R::rgamma(shape[c], 1.0);
               sumg += g[c];
          }
          for (int c = 0; c < C; ++c) omega(k, c) = g[c] / sumg;
     }

     return omega;
}

double sample_alpha_cpp(const arma::umat& xi,
                        const double      a_alpha,
                        const double      b_alpha,
                        const double      alpha,
                        const int         K,
                        const int         n) {
     
     // Escobar–West (multi-group) auxiliary scheme:
     // For each layer k: η_k ~ Beta(α+1, n), s_k ~ Bernoulli(n/(n+α))
     // Then: α ~ Gamma(a_alpha + Σ_k (m_k - s_k),  b_alpha - Σ_k log η_k)
     
     double sum_log_eta = 0.0;
     int    sum_m       = 0;
     int    sum_s       = 0;
     
     for (int k = 0; k < K; ++k) {
          // Occupied clusters m_k in layer k
          arma::uword Ck_infer = xi.col(k).max() + 1u;
          arma::ivec counts(Ck_infer, arma::fill::zeros);
          for (int i = 0; i < n; ++i) counts[ xi(i, k) ]++;
          int m_k = static_cast<int>(arma::sum(counts > 0));
          
          // Auxiliary variables
          const double eta_k = R::rbeta(alpha + 1.0, static_cast<double>(n));
          const double p_k   = static_cast<double>(n) / (static_cast<double>(n) + alpha);
          const int    s_k   = (R::runif(0.0, 1.0) < p_k) ? 1 : 0;
          
          sum_log_eta += std::log(eta_k);
          sum_m       += m_k;
          sum_s       += s_k;
     }
     
     const double A = a_alpha + static_cast<double>(sum_m - sum_s);
     const double B  = b_alpha - sum_log_eta;
     
     return R::rgamma(A, 1.0 / B);
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

double sample_rho2_cpp(const arma::cube& Theta,
                       const double      a_rho,
                       const double      b_rho,
                       const int         K,
                       const int         C) {
     
     const int M  = K * C * (C + 1) / 2;
     double    ss = 0.0;

     for (int k = 0; k < K; ++k) {
          for (int a = 0; a < C; ++a) {
               for (int b = a; b < C; ++b) {
                    const double th = Theta(a, b, k);
                    ss += th * th;
               }
          }
     }
     
     const double A = a_rho + 0.5 * static_cast<double>(M);
     const double B = b_rho + 0.5 * ss;
     
     return 1.0 / R::rgamma(A, 1.0 / B);
}

// [[Rcpp::export]]
Rcpp::List gibbs_sampler_multilayer(const arma::cube& y,
                                    const int         C,
                                    arma::umat        xi,
                                    const int n_iter, const int n_burn, const int n_thin,
                                    const double a_omega,    const double b_omega,
                                    const double a_sigma,    const double b_sigma,
                                    const double a_tau,      const double b_tau,
                                    const double a_kappa,    const double b_kappa,
                                    const double a_rho,      const double b_rho,
                                    const double a_alpha,    const double b_alpha) {
     // Dimensions
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     // Initialize
     double omega2    = 1.0 / R::rgamma(a_omega,    1.0 / b_omega);
     double sigma2    = 1.0 / R::rgamma(a_sigma,    1.0 / b_sigma);
     double tau2      = 1.0 / R::rgamma(a_tau,      1.0 / b_tau);
     double kappa2    = 1.0 / R::rgamma(a_kappa,    1.0 / b_kappa);
     double rho2      = 1.0 / R::rgamma(a_rho,      1.0 / b_rho);
     
     double zeta  = R::rnorm(0.0, std::sqrt(omega2));
     double alpha = R::rgamma(a_alpha, 1.0 / b_alpha);
     
     arma::vec mu(K);
     for (int k = 0; k < K; ++k) mu[k] = R::rnorm(0.0, std::sqrt(sigma2));
     
     arma::vec vartheta(n);
     for (int i = 0; i < n; ++i) vartheta[i] = R::rnorm(0.0, std::sqrt(kappa2));

     arma::mat  delta(n, K, arma::fill::zeros);
     arma::cube z(n, n, K, arma::fill::zeros);
     
     arma::mat omega(K, C, arma::fill::ones);
     for (int k = 0; k < K; ++k) omega.row(k) /= arma::accu(omega.row(k));
     
     arma::cube Theta(C, C, K, arma::fill::zeros);
     for (int k = 0; k < K; ++k) {
          for (int a = 0; a < C; ++a) {
               for (int b = a; b < C; ++b) {
                    const double th = R::rnorm(0.0, std::sqrt(rho2));
                    Theta(a, b, k) = th;
                    Theta(b, a, k) = th;
               }
          }
     }
     
     // Storage
     const int n_samples = (n_iter - n_burn) / n_thin;
     
     arma::cube  store_delta     (n_samples, n, K, arma::fill::zeros);
     arma::cube  store_Theta     (C, C, K * n_samples, arma::fill::zeros);
     
     arma::mat   store_mu        (n_samples, K, arma::fill::zeros);
     arma::mat   store_vartheta  (n_samples, n, arma::fill::zeros);
     arma::cube  store_omega     (n_samples, K, C, arma::fill::zeros); 
     arma::ucube store_xi        (n_samples, n, K, arma::fill::zeros);
     
     arma::vec   store_zeta      (n_samples, arma::fill::zeros);
     arma::vec   store_alpha     (n_samples, arma::fill::zeros);
     arma::vec   store_omega2    (n_samples, arma::fill::zeros);
     arma::vec   store_sigma2    (n_samples, arma::fill::zeros);
     arma::vec   store_tau2      (n_samples, arma::fill::zeros);
     arma::vec   store_kappa2    (n_samples, arma::fill::zeros);
     arma::vec   store_rho2      (n_samples, arma::fill::zeros);
     
     // Sampling
     Rcpp::Rcout << "Initializing Gibbs sampler...\n";
     const int step = std::max(1, n_iter / 20);
     
     for (int iter = 1; iter <= n_iter; ++iter) {
          
          // Update model parameters (UPDATED signatures use layer-specific xi)
          z         = sample_z_cpp(y, mu, delta, Theta, xi, zeta, z, K, n);
          zeta      = sample_zeta_cpp(z, mu, delta, Theta, xi, omega2, K, n);
          mu        = sample_mu_cpp(z, mu, delta, Theta, xi, zeta, sigma2, K, n);
          delta     = sample_delta_cpp(z, mu, tau2, delta, vartheta, Theta, xi, zeta, K, n);
          vartheta  = sample_vartheta_cpp(delta, vartheta, kappa2, tau2, n, K);
          Theta     = sample_Theta_cpp(z, mu, delta, xi, Theta, rho2, zeta, K, n, C);
          xi        = sample_xi_cpp(z, mu, delta, xi, Theta, zeta, alpha, K, n, C);
          omega     = sample_omega_cpp(xi, alpha, K, n, C);
          alpha     = sample_alpha_cpp(xi, a_alpha, b_alpha, alpha, K, n);
          omega2    = sample_omega2_cpp(zeta, a_omega, b_omega);
          sigma2    = sample_sigma2_cpp(mu, a_sigma, b_sigma, K);
          tau2      = sample_tau2_cpp(delta, vartheta, a_tau, b_tau, n, K);
          kappa2    = sample_kappa2_cpp(vartheta, a_kappa, b_kappa, n);
          rho2      = sample_rho2_cpp(Theta, a_rho, b_rho, K, C);
          
          // Store
          if ((iter > n_burn) && ((iter - n_burn) % n_thin == 0)) {
               const int pos = (iter - n_burn) / n_thin - 1;
               
               // delta
               for (int k = 0; k < K; ++k)
                    for (int i = 0; i < n; ++i)
                         store_delta(pos, i, k) = delta(i, k);
               
               // Theta: stack K slices per sample
               for (int k = 0; k < K; ++k) {
                    store_Theta.slice(pos * K + k) = Theta.slice(k);
               }
               
               // xi (n x K per sample)
               for (int k = 0; k < K; ++k)
                    for (int i = 0; i < n; ++i)
                         store_xi(pos, i, k) = xi(i, k);
               
               // omega (K x C per sample)
               for (int k = 0; k < K; ++k)
                    for (int c = 0; c < C; ++c)
                         store_omega(pos, k, c) = omega(k, c);
          
               // vectors/scalars
               store_mu.row(pos)       = mu.t();
               store_vartheta.row(pos) = vartheta.t();
          
               store_zeta     [pos]    = zeta;
               store_alpha    [pos]    = alpha;
               store_omega2   [pos]    = omega2;
               store_sigma2   [pos]    = sigma2;
               store_tau2     [pos]    = tau2;
               store_kappa2   [pos]    = kappa2;
               store_rho2     [pos]    = rho2;
          }
     
          // Progress (print m_k for each layer)
          if ((iter % step) == 0) {
               int pct = static_cast<int>(std::round(100.0 * iter / n_iter));
               
               arma::ivec mK(K, arma::fill::zeros);
               for (int k = 0; k < K; ++k) {
                    arma::uword C_infer = xi.col(k).max() + 1u;
                    arma::ivec counts_k(C_infer, arma::fill::zeros);
                    for (int i = 0; i < n; ++i) counts_k[ xi(i, k) ]++;
                    mK[k] = static_cast<int>(arma::sum(counts_k > 0));
               }
          
               Rcpp::Rcout << "Progress: " << pct << "% completed"
                           << " | occupied clusters per layer m = ";
               for (int k = 0; k < K; ++k) {
                    Rcpp::Rcout << mK[k] << (k+1<K ? ' ' : '\n');
               }
          }
     }
     
     Rcpp::Rcout << "Sampler completed.\n";
     
     return Rcpp::List::create(
          Rcpp::Named("zeta")       = store_zeta,
          Rcpp::Named("mu")         = store_mu,
          Rcpp::Named("delta")      = store_delta,
          Rcpp::Named("vartheta")   = store_vartheta,
          Rcpp::Named("Theta")      = store_Theta, 
          Rcpp::Named("xi")         = store_xi,
          Rcpp::Named("omega")      = store_omega,
          Rcpp::Named("alpha")      = store_alpha, 
          Rcpp::Named("omega2")     = store_omega2,
          Rcpp::Named("sigma2")     = store_sigma2,
          Rcpp::Named("tau2")       = store_tau2,
          Rcpp::Named("kappa2")     = store_kappa2,
          Rcpp::Named("rho2")       = store_rho2
     );
}

// [[Rcpp::export]]
arma::vec log_likelihood_multilayer_cpp(const arma::cube& y,
                                        const Rcpp::List& samples) {
     
     const int n = y.n_rows; 
     const int K = y.n_slices;
     
     arma::vec   zeta   = samples["zeta"];
     arma::mat   mu     = samples["mu"];
     arma::cube  delta  = samples["delta"];
     arma::cube  Theta  = samples["Theta"]; 
     arma::ucube xi_all = samples["xi"];
     
     const int n_samples = mu.n_rows;
     arma::vec log_lik_samples(n_samples, arma::fill::zeros);
     
     const double eps = 1e-10;
     
     for (int s = 0; s < n_samples; ++s) {
          double       log_lik = 0.0;
          const double zeta_s  = zeta[s];
          
          for (int k = 0; k < K; ++k) {
               const double       mu_k  = mu(s, k);
               const arma::rowvec del_k = delta.slice(k).row(s);
               
               for (int i = 0; i < n - 1; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                         const unsigned int ci = xi_all(s, i, k); 
                         const unsigned int cj = xi_all(s, j, k);
                         const double       tabk = Theta(ci, cj, s * K + k);
                         
                         const double eta_ijk = mu_k + zeta_s + del_k[i] + del_k[j] + tabk;
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
                               const arma::cube& Theta, 
                               const arma::umat& xi,
                               const double      zeta,
                               const double      eps = 1e-10) {
     
     const int n = y.n_rows;
     const int K = y.n_slices;
     
     double log_lik = 0.0;
     
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const unsigned int ci   = xi(i, k); 
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const unsigned int cj   = xi(j, k); 
                    const double       tabk = Theta(ci, cj, k);
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + tabk;
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
                                       const arma::cube& Theta,
                                       const arma::umat& xi,
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
               const double       d_ik = del_k[i];
               const unsigned int ci   = xi(i, k);
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const unsigned int cj   = xi(j, k); 
                    const double       tabk = Theta(ci, cj, k); 
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + tabk;
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
                                const arma::cube& Theta,
                                const arma::umat& xi,
                                const double      zeta) {
     
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     arma::cube P(n, n, K, arma::fill::zeros);
     
     // Build probabilities (upper triangle; mirror to ensure symmetry)
     for (int k = 0; k < K; ++k) {
          const double    mu_k  = mu[k];
          const arma::vec del_k = delta.col(k);
          
          for (int i = 0; i < n - 1; ++i) {
               const double       d_ik = del_k[i];
               const unsigned int ci   = xi(i, k); 
               
               for (int j = i + 1; j < n; ++j) {
                    const double       d_jk = del_k[j];
                    const unsigned int cj   = xi(j, k);
                    const double       tabk = Theta(ci, cj, k);
                    
                    const double eta_ijk = mu_k + zeta + d_ik + d_jk + tabk;
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
                                           const arma::cube& Theta,
                                           const arma::umat& xi,
                                           const double      zeta) {
     
     const int n = delta.n_rows;
     const int K = mu.n_elem;
     
     arma::cube Y(n, n, K, arma::fill::zeros);
     
     arma::cube P = interaction_prob_cpp(mu, delta, Theta, xi, zeta);
     
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