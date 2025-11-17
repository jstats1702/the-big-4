prior_sim_with_covariates <- function(X, 
                                      n_sim = 1e6,
                                      a_sigma    = 3, b_sigma    = 2 * 1.5,
                                      a_kappa    = 3, b_kappa    = 2 * 1.5,
                                      a_tau      = 3, b_tau      = 2 * 1.5,
                                      a_varsigma = 3, b_varsigma = 2 * 100,
                                      a_omega    = 3, b_omega    = 2 * 1.5) {
     
     stopifnot(length(dim(X)) == 3)
     n <- dim(X)[1]; p <- dim(X)[3]
     up <- upper.tri(matrix(0, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Collect dyad covariates into an M x p matrix (rows are dyads i<j)
     X_dyads <- matrix(NA_real_, nrow = M, ncol = p)
     for (s in seq_len(p)) X_dyads[, s] <- X[,,s][up]
     
     # ---- Draw hyper-variances ~ IG(a,b)  ===  1 / Gamma(shape=a, rate=1/b) ----
     sigma2    <- 1 / rgamma(n_sim, shape = a_sigma,    rate = 1 / b_sigma)
     kappa2    <- 1 / rgamma(n_sim, shape = a_kappa,    rate = 1 / b_kappa)
     tau2      <- 1 / rgamma(n_sim, shape = a_tau,      rate = 1 / b_tau)
     varsigma2 <- 1 / rgamma(n_sim, shape = a_varsigma, rate = 1 / b_varsigma)
     omega2    <- 1 / rgamma(n_sim, shape = a_omega,    rate = 1 / b_omega)
     
     # ---- Baseline components ----
     mu         <- rnorm(n_sim, mean = 0, sd = sqrt(sigma2))
     zeta       <- rnorm(n_sim, mean = 0, sd = sqrt(omega2))   # NEW global random effect
     vartheta_i <- rnorm(n_sim, mean = 0, sd = sqrt(kappa2))
     vartheta_j <- rnorm(n_sim, mean = 0, sd = sqrt(kappa2))
     delta_i    <- rnorm(n_sim, mean = vartheta_i, sd = sqrt(tau2))
     delta_j    <- rnorm(n_sim, mean = vartheta_j, sd = sqrt(tau2))
     
     # ---- Covariate contribution x^T beta ----
     # Sample a random dyad per iteration (uniform over observed dyads)
     idx   <- sample.int(M, size = n_sim, replace = TRUE)
     Xpick <- X_dyads[idx, , drop = FALSE]           # n_sim x p
     
     # Draw beta_k ~ N(0, varsigma2 I_p) for each iteration (independent across p)
     beta_mat <- matrix(rnorm(n_sim * p), nrow = n_sim, ncol = p) * sqrt(varsigma2)
     xbeta    <- rowSums(beta_mat * Xpick)           # x_{ij}^T beta for the picked dyad
     
     # ---- Prior-predictive probability ----
     eta <- mu + zeta + delta_i + delta_j + xbeta    # NEW: + zeta
     p   <- pnorm(eta)
     
     # ---- Summary ----
     cat("Summary of p (prior with covariates + zeta):\n")
     print(c(
          mean   = mean(p),
          sd     = sd(p),
          q2.5   = unname(quantile(p, 0.025)),
          median = unname(quantile(p, 0.5)),
          q97.5  = unname(quantile(p, 0.975))
     ))
     cat(sprintf("\nMass at extremes: P(p < 0.05) = %.3f, P(p > 0.95) = %.3f\n",
                 mean(p < 0.05), mean(p > 0.95)))
     
     # ---- Optional diagnostics ----
     # Average squared norm of dyads (should be ~ p if X is standardized)
     qbar <- mean(rowSums(X_dyads^2))
     cat(sprintf("\nqbar (mean ||x_ij||^2 over dyads): %.3f;  E[varsigma^2] = b/(a-1) = %.5f;  E[omega^2] = %.5f\n",
                 qbar,
                 b_varsigma / (a_varsigma - 1),
                 b_omega  / (a_omega  - 1)))
     
     hist(p, freq = FALSE,
          main = "",
          xlab = expression(theta),
          ylab = "Density",
          border = "grey85", col = "grey85")
     
     invisible(list(p = p, qbar = qbar))
}

load("data_metallica.RData")

X <- cov_metallica$X

set.seed(123)

a_sigma    <- 3; b_sigma    <- 2 * 1.5
a_kappa    <- 3; b_kappa    <- 2 * 1.5
a_tau      <- 3; b_tau      <- 2 * 1.5
a_varsigma <- 3; b_varsigma <- 2 * 100
a_omega    <- 3; b_omega    <- 2 * 1.5


pdf(file = paste0("prior_simulation_model_2.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

res2 <- prior_sim_with_covariates(
     X,
     n_sim = 1e6,
     a_sigma    = a_sigma,    b_sigma    = b_sigma,
     a_kappa    = a_kappa,    b_kappa    = b_kappa,
     a_tau      = a_tau,      b_tau      = b_tau,
     a_varsigma = a_varsigma, b_varsigma = b_varsigma,
     a_omega    = a_omega,    b_omega    = b_omega
)

dev.off()