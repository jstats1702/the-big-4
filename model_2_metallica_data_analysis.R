# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

Rcpp::sourceCpp("model_2_mcmc_algorithm.cpp")

source("helper_functions.R")

library(igraph)
library(ggplot2)
library(gridExtra)
library(coda)

# Data and posterior samples ===================================================

load("data_metallica.RData")

load("samples_model_2_metallica.RData")

catalog <- res_metallica$logRMS$catalog

# Log-likelihood ===============================================================

x11()
plot_log_likelihood (log_lik)

# Mean effective sample sizes ==================================================

MCE <- c(
     mce(samples$zeta),
     mean(apply(samples$mu,       2,       mce)),
     mean(apply(samples$delta,    c(2, 3), mce)),
     mean(apply(samples$vartheta, 2,       mce)),
     mean(apply(samples$beta,     c(2,3),  mce)),
     mce(samples$omega2),
     mce(samples$sigma2),
     mce(samples$tau2),
     mce(samples$kappa2),
     mce(samples$varsigma2)
)
names(MCE) <- c("zeta", "mu", "delta", "vartheta", "beta", 
                "omega2", "sigma2", "tau2", "kappa2", "varsigma2")

round(MCE, 5)

# Inference on band random effect ==============================================

tab <- summ_vec(samples$zeta)
names(tab) <- c("Mean", "SD", "Q2.5", "Q97.5")

round(tab, 3)

# Inference on variance components =============================================

tab <- rbind(
     summ_vec(sqrt(samples$omega2)),
     summ_vec(sqrt(samples$sigma2)),
     summ_vec(sqrt(samples$tau2)),
     summ_vec(sqrt(samples$kappa2)),
     summ_vec(sqrt(samples$varsigma2))
)
rownames(tab) <- c("omega", "sigma", "tau", "kappa", "varsigma")
colnames(tab) <- c("Mean", "SD", "Q2.5", "Q97.5")

round(tab, 3)

# Inference on intercepts parameters ===========================================

x11()
plot_mu_by_layer(samples$mu)

tab <- cbind(
     apply(samples$mu, 2, mean),
     apply(samples$mu, 2, sd),
     t(apply(samples$mu, 2, quantile, probs = c(0.025, 0.975)))
)
rownames(tab) <- c("Loudness", "Brightness", "Tonality", "Rhythm")
colnames(tab) <- c("Mean", "SD", "Q2.5", "Q97.5")

round(tab, 3)

# Inference on sociality parameters ============================================

x11()
plot_delta_by_layer(samples$delta)

# Inference on mean sociality parameters =======================================

x11()
plot_vartheta(samples$vartheta)

# Inference on covariates coefficients =========================================

x11()
plot_beta_by_layer(samples$beta)

# Songs with significant effects ===============================================

summarize_vartheta(catalog, samples)

# Correlation between layers ===================================================

layers_correlation(samples$delta)

# Interaction probabilities ====================================================

metrics <- prob_metrics_per_layer_model_2(samples, X, y)

round(metrics$auc$mean,     3)
round(metrics$auprc$mean,   3)
round(metrics$brier$mean,   3)
round(metrics$logloss$mean, 3)

# Information criteria =========================================================

dic_model_2(log_lik, samples, y, X)

waic_model_2(samples, y, X)

# Test statistics ==============================================================

x11()
plot_test_stats_model_2(samples, y, X, layer_titles = c(1,2,3,4))

# Originality score ============================================================

originality_by_layer(samples)

# End ==========================================================================