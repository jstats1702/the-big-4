prior_sim_model3 <- function(X, 
                             n_sim = 1e6, 
                             d = 2,
                             a_sigma,    b_sigma,
                             a_kappa,    b_kappa, 
                             a_tau,      b_tau,
                             a_varsigma, b_varsigma, 
                             a_upsilon,  b_upsilon,
                             a_omega,    b_omega) {
     stopifnot(length(dim(X)) == 3)
     n <- dim(X)[1]; p <- dim(X)[3]
     up <- upper.tri(matrix(0, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Collect dyad covariates into an M x p matrix (rows are dyads i<j)
     X_dyads <- matrix(NA_real_, nrow = M, ncol = p)
     for (s in seq_len(p)) X_dyads[, s] <- X[ , , s][up]
     
     # ---- Hyper-variances ~ IG(a,b)  (1 / Gamma(shape=a, rate=1/b)) ----
     sigma2    <- 1 / rgamma(n_sim, shape = a_sigma,    rate = 1 / b_sigma)
     kappa2    <- 1 / rgamma(n_sim, shape = a_kappa,    rate = 1 / b_kappa)
     tau2      <- 1 / rgamma(n_sim, shape = a_tau,      rate = 1 / b_tau)
     varsigma2 <- 1 / rgamma(n_sim, shape = a_varsigma, rate = 1 / b_varsigma)
     upsilon2  <- 1 / rgamma(n_sim, shape = a_upsilon,  rate = 1 / b_upsilon)
     omega2    <- 1 / rgamma(n_sim, shape = a_omega,    rate = 1 / b_omega)
     
     # ---- Baseline components ----
     mu         <- rnorm(n_sim, 0, sqrt(sigma2))
     vartheta_i <- rnorm(n_sim, 0, sqrt(kappa2))
     vartheta_j <- rnorm(n_sim, 0, sqrt(kappa2))
     delta_i    <- rnorm(n_sim, vartheta_i, sqrt(tau2))
     delta_j    <- rnorm(n_sim, vartheta_j, sqrt(tau2))
     zeta       <- rnorm(n_sim, 0, sqrt(omega2))  # global random shift
     
     # ---- Pick random dyads (uniform over observed dyads) ----
     idx   <- sample.int(M, size = n_sim, replace = TRUE)
     Xpick <- X_dyads[idx, , drop = FALSE]      # n_sim x p
     
     # ---- beta_k ~ N(0, varsigma2 I_p) per iteration ----
     beta_mat <- matrix(rnorm(n_sim * p), nrow = n_sim, ncol = p) * sqrt(varsigma2)
     xbeta    <- rowSums(beta_mat * Xpick)      # x_{ij}^T beta
     
     # ---- Bilinear term: lambda * (u_i^T u_j) ----
     # U_i, U_j ~ N_d(0, I), independent; lambda ~ N(0, upsilon2)
     Ui  <- matrix(rnorm(n_sim * d), nrow = n_sim, ncol = d)
     Uj  <- matrix(rnorm(n_sim * d), nrow = n_sim, ncol = d)
     sij <- rowSums(Ui * Uj)                    # u_i^T u_j per iteration
     lambda   <- rnorm(n_sim, 0, sqrt(upsilon2))
     bilinear <- lambda * sij
     
     # ---- Prior-predictive probability ----
     eta <- zeta + mu + delta_i + delta_j + xbeta + bilinear
     p   <- pnorm(eta)
     
     # ---- Summary ----
     cat("Summary of p (Model 3 prior with covariates + bilinear + zeta):\n")
     print(c(
          mean   = mean(p),
          sd     = sd(p),
          q2.5   = unname(quantile(p, 0.025)),
          median = unname(quantile(p, 0.5)),
          q97.5  = unname(quantile(p, 0.975))
     ))
     cat(sprintf("\nMass at extremes: P(p < 0.05) = %.3f, P(p > 0.95) = %.3f\n",
                 mean(p < 0.05), mean(p > 0.95)))
     
     # ---- Diagnostics ----
     qbar <- mean(rowSums(X_dyads^2))           # ~ p if X standardized
     cat(sprintf("\nqbar (mean ||x_ij||^2 over dyads): %.3f;  E[varsigma^2] = %.5f\n",
                 qbar, b_varsigma / (a_varsigma - 1)))
     cat(sprintf("E[upsilon^2] = b_u/(a_u-1) = %.5f;  mean(s_ij^2) ≈ %.3f (theory: d)\n",
                 b_upsilon / (a_upsilon - 1), mean(sij^2)))
     cat(sprintf("E[omega^2] = b_omega/(a_omega-1) = %.5f;  sd(zeta) ≈ %.3f\n",
                 b_omega / (a_omega - 1), sd(zeta)))
     cat(sprintf("sd(lambda) ≈ %.3f;  sd(bilinear) ≈ %.3f\n",
                 sd(lambda), sd(bilinear)))
     
     hist(p, freq = FALSE,
          main = "",
          xlab = expression(theta),
          ylab = "Density",
          border = "grey85", col = "grey85")
     
     invisible(list(p = p,
                    eta = eta,
                    zeta = zeta,
                    xbeta = xbeta,
                    bilinear = bilinear,
                    sij = sij,
                    lambda = lambda,
                    qbar = qbar,
                    sigma2 = sigma2,
                    kappa2 = kappa2,
                    tau2 = tau2,
                    varsigma2 = varsigma2,
                    upsilon2 = upsilon2,
                    omega2 = omega2))
}

load("data_metallica.RData")

X <- cov_metallica$X

set.seed(123)

pdf(file = paste0("prior_simulation_model_3.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

res3 <- prior_sim_model3(
     X,
     n_sim = 1e6,
     d = 5,
     a_sigma    = 3, b_sigma    = 2*1.5,
     a_kappa    = 3, b_kappa    = 2*1.5,
     a_tau      = 3, b_tau      = 2*1.5,
     a_varsigma = 3, b_varsigma = 2*100,
     a_upsilon  = 3, b_upsilon  = 2*100,
     a_omega    = 3, b_omega    = 2*1.5
)

dev.off()