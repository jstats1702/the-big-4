# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

source("helper_functions.R")

Rcpp::sourceCpp("model_1_mcmc_algorithm_b.cpp")

# Data =========================================================================

load("data_micro.RData")

n <- I
K <- J
y <- Ycube

rm(I, J, Ycube)

# Hyperparameters ==============================================================

a_omega <- 3; b_omega <- 2 * 1.5  # zeta
a_tau   <- 3; b_tau   <- 2 * 1.5  # delta
a_kappa <- 3; b_kappa <- 2 * 1.5  # vartheta

# MCMC settings ================================================================

n_iter <- 1000000 + 200000
n_burn <- 200000
n_thin <- 20
n_keep <- (n_iter - n_burn) / n_thin

# Run sampler ==================================================================

set.seed(42)
samples <- gibbs_sampler_multilayer(
     y, 
     n_iter, n_burn, n_thin,
     a_omega, b_omega,
     a_tau,   b_tau,
     a_kappa, b_kappa
)

# Compute log-likelihood chain
log_lik <- log_likelihood_multilayer_cpp(y, samples)

# Save
save(y, samples, log_lik, file = "samples_model_1_micro.RData")

# End ==========================================================================