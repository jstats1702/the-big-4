# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

source("helper_functions_b.R")

Rcpp::sourceCpp("model_4_mcmc_algorithm_b.cpp")

# Data =========================================================================

load("data_aarhus.RData")

n <- I
K <- J
y <- Ycube

rm(I, J, Ycube)

# Hyperparameters ==============================================================

d          <- 3
a_omega    <- 3; b_omega    <- 2 * 1.5  # zeta
a_sigma    <- 3; b_sigma    <- 2 * 1.5  # mu
a_tau      <- 3; b_tau      <- 2 * 1.5  # delta
a_kappa    <- 3; b_kappa    <- 2 * 1.5  # vartheta
a_upsilon  <- 3; b_upsilon  <- 2 * 1    # lambda

# MCMC settings ================================================================

n_iter <- 1000000 + 200000
n_burn <- 200000
n_thin <- 20
n_keep <- (n_iter - n_burn) / n_thin

# Run sampler ==================================================================

set.seed(42)
samples <- gibbs_sampler_multilayer(
     y,
     d,
     n_iter, n_burn, n_thin,
     a_omega,    b_omega,
     a_sigma,    b_sigma,
     a_tau,      b_tau,
     a_kappa,    b_kappa,
     a_upsilon,  b_upsilon
)

# Compute log-likelihood chain
log_lik <- log_likelihood_multilayer_cpp(y, samples)

# Save
save(y, samples, log_lik, file = "samples_model_4_aarhus.RData")

# End ==========================================================================