# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

source("helper_functions.R")

Rcpp::sourceCpp("model_3_mcmc_algorithm.cpp")

# Data =========================================================================

load("data_anthrax.RData")

# Build 3D adjacency array y (n x n x K)
metrics <- c("logRMS", "logSC", "logitSFM", "logFlux")
K <- length(metrics)
n <- nrow(res_anthrax$logRMS$adjacency)

y <- array(NA_real_, dim = c(n, n, K))
y[, , 1] <- res_anthrax$logRMS$adjacency
y[, , 2] <- res_anthrax$logSC$adjacency
y[, , 3] <- res_anthrax$logitSFM$adjacency
y[, , 4] <- res_anthrax$logFlux$adjacency

X <- cov_anthrax$X

# Hyperparameters ==============================================================

d          <- 3
a_omega    <- 3; b_omega    <- 2 * 1.5  # zeta
a_sigma    <- 3; b_sigma    <- 2 * 1.5  # mu
a_tau      <- 3; b_tau      <- 2 * 1.5  # delta
a_kappa    <- 3; b_kappa    <- 2 * 1.5  # vartheta
a_varsigma <- 3; b_varsigma <- 2 * 100  # beta
a_upsilon  <- 3; b_upsilon  <- 2 * 100  # lambda

# MCMC settings ================================================================

n_iter <- 1000000 + 200000
n_burn <- 200000
n_thin <- 20
n_keep <- (n_iter - n_burn) / n_thin

# Run sampler ==================================================================

set.seed(42)
samples <- gibbs_sampler_multilayer(
     y,
     X,
     d,
     n_iter, n_burn, n_thin,
     a_omega,    b_omega,
     a_sigma,    b_sigma,
     a_tau,      b_tau,
     a_kappa,    b_kappa,
     a_varsigma, b_varsigma,
     a_upsilon,  b_upsilon
)

# Compute log-likelihood chain
log_lik <- log_likelihood_multilayer_cpp(y, X, samples)

# Save
save(y, X, samples, log_lik, file = "samples_model_3_anthrax.RData")

# End ==========================================================================