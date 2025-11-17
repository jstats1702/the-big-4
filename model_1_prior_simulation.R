set.seed(123)

# Hyperparameters (Inverse-Gamma with scale parameter b)
a_sigma <- 3; b_sigma <- 2 * 1.5   # mu
a_kappa <- 3; b_kappa <- 2 * 1.5   # vartheta
a_tau   <- 3; b_tau   <- 2 * 1.5   # delta
a_omega <- 3; b_omega <- 2 * 1.5   # zeta

n_sim <- 1000000L

# Draw variances ~ Inverse-Gamma(a, b)  ===  1 / Gamma(shape = a, rate = 1/b)
sigma2 <- 1 / rgamma(n_sim, shape = a_sigma, rate = 1 / b_sigma)
kappa2 <- 1 / rgamma(n_sim, shape = a_kappa, rate = 1 / b_kappa)
tau2   <- 1 / rgamma(n_sim, shape = a_tau,   rate = 1 / b_tau)
omega2 <- 1 / rgamma(n_sim, shape = a_omega, rate = 1 / b_omega)

# Global layer effects
mu <- rnorm(n_sim, mean = 0, sd = sqrt(sigma2))

# Band random effect
zeta <- rnorm(n_sim, mean = 0, sd = sqrt(omega2))

# Node-specific baselines
vartheta_i <- rnorm(n_sim, mean = 0, sd = sqrt(kappa2))
vartheta_j <- rnorm(n_sim, mean = 0, sd = sqrt(kappa2))

# Sociability effects per layer (conditional on vartheta and tau2)
delta_i <- rnorm(n_sim, mean = vartheta_i, sd = sqrt(tau2))
delta_j <- rnorm(n_sim, mean = vartheta_j, sd = sqrt(tau2))

# Prior-predictive connection probability
eta <- mu + zeta + delta_i + delta_j
p   <- pnorm(eta)

# Numerical summary
cat("Summary of p (prior):\n")
print(c(
     mean   = mean(p),
     sd     = sd(p),
     q2.5   = unname(quantile(p, 0.025)),
     median = unname(quantile(p, 0.5)),
     q97.5  = unname(quantile(p, 0.975))
))
cat(sprintf("\nMass at extremes: P(p < 0.05) = %.3f, P(p > 0.95) = %.3f\n",
            mean(p < 0.05), mean(p > 0.95)))

# Histogram
pdf(file = paste0("prior_simulation_model_1.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

hist(p, freq = FALSE,
     main = "",
     xlab = expression(theta),
     ylab = "Density",
     border = "grey85", col = "grey85")

dev.off()
