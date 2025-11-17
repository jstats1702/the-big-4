prior_sim_model4 <- function(X, 
                             n_sim = 1e6, 
                             d = 2,
                             a_sigma,    b_sigma,
                             a_kappa,    b_kappa, 
                             a_tau,      b_tau,
                             a_varsigma, b_varsigma, 
                             a_upsilon,  b_upsilon) {
     stopifnot(length(dim(X)) == 3)
     n <- dim(X)[1]; p <- dim(X)[3]
     up <- upper.tri(matrix(0, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Collect dyad covariates into an M x p matrix (rows are dyads i<j)
     X_dyads <- matrix(NA_real_, nrow = M, ncol = p)
     for (s in seq_len(p)) X_dyads[, s] <- X[,,s][up]
     
     # ---- Hyper-variances (draw ONCE per prior replicate; here: one replicate) ----
     sigma2    <- 1 / rgamma(1, shape = a_sigma,    rate = 1 / b_sigma)
     kappa2    <- 1 / rgamma(1, shape = a_kappa,    rate = 1 / b_kappa)
     tau2      <- 1 / rgamma(1, shape = a_tau,      rate = 1 / b_tau)
     varsigma2 <- 1 / rgamma(1, shape = a_varsigma, rate = 1 / b_varsigma)
     upsilon2  <- 1 / rgamma(1, shape = a_upsilon,  rate = 1 / b_upsilon)
     
     # ---- Baseline components ----
     mu         <- rnorm(n_sim, 0, sqrt(sigma2))
     vartheta_i <- rnorm(n_sim, 0, sqrt(kappa2))
     vartheta_j <- rnorm(n_sim, 0, sqrt(kappa2))
     delta_i    <- rnorm(n_sim, vartheta_i, sqrt(tau2))
     delta_j    <- rnorm(n_sim, vartheta_j, sqrt(tau2))
     
     # ---- Pick random dyads ----
     idx   <- sample.int(M, size = n_sim, replace = TRUE)
     Xpick <- X_dyads[idx, , drop = FALSE]      # n_sim x p
     
     # ---- beta ~ N(0, varsigma2 I_p) per draw ----
     beta_mat <- matrix(rnorm(n_sim * p), nrow = n_sim, ncol = p) * sqrt(varsigma2)
     xbeta    <- rowSums(beta_mat * Xpick)
     
     # ---- Distance term: -exp(lambda) * ||u_i - u_j|| (Model 4) ----
     Ui  <- matrix(rnorm(n_sim * d), nrow = n_sim, ncol = d)
     Uj  <- matrix(rnorm(n_sim * d), nrow = n_sim, ncol = d)
     dij <- sqrt(rowSums((Ui - Uj)^2))
     lambda <- rnorm(n_sim, 0, sqrt(upsilon2))     # upsilon2 fixed across draws
     coeff  <- -exp(lambda)
     dist_term <- coeff * dij
     
     # ---- Prior-predictive probability ----
     eta    <- mu + delta_i + delta_j + xbeta + dist_term
     p_edge <- pnorm(eta)
     
     # ---- Summary ----
     cat("Summary of p (Model 4 prior with covariates + distance decay):\n")
     print(c(
          mean   = mean(p_edge),
          sd     = sd(p_edge),
          q2.5   = unname(quantile(p_edge, 0.025)),
          median = unname(quantile(p_edge, 0.5)),
          q97.5  = unname(quantile(p_edge, 0.975))
     ))
     cat(sprintf("\nMass at extremes: P(p < 0.05) = %.3f, P(p > 0.95) = %.3f\n",
                 mean(p_edge < 0.05), mean(p_edge > 0.95)))
     
     # ---- Diagnostics ----
     qbar <- mean(rowSums(X_dyads^2))           # ~ p if X standardized
     cat(sprintf("\nqbar (mean ||x_ij||^2 over dyads): %.3f;  E[varsigma^2] ≈ %.5f\n",
                 qbar, b_varsigma / (a_varsigma - 1)))
     cat(sprintf("Drawn upsilon2 = %.6g (E ≈ %.5f);  mean(d_ij^2) ≈ %.3f (theory: 2*d)\n",
                 upsilon2, b_upsilon / (a_upsilon - 1), mean(dij^2)))
     cat(sprintf("sd(lambda) ≈ %.3f;  sd(distance_term) ≈ %.3f\n",
                 sd(lambda), sd(dist_term)))
     
     hist(p_edge, freq = FALSE,
          main = "",
          xlab = expression(theta),
          ylab = "Density",
          border = "grey85", col = "grey85")
     
     invisible(list(p = p_edge,
                    eta = eta,
                    xbeta = xbeta,
                    dist_term = dist_term,
                    dij = dij,
                    lambda = lambda,
                    upsilon2 = upsilon2,
                    qbar = qbar))
}

load("data_metallica.RData")

X <- cov_metallica$X

set.seed(123)

pdf(file = paste0("prior_simulation_model_4.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

res4 <- prior_sim_model4(
     X,
     n_sim = 1e6,
     d = 3,
     a_sigma    = 3, b_sigma    = 2*1.5,
     a_kappa    = 3, b_kappa    = 2*1.5,
     a_tau      = 3, b_tau      = 2*1.5,
     a_varsigma = 3, b_varsigma = 2*100,
     a_upsilon  = 3, b_upsilon  = 2*1
)

dev.off()