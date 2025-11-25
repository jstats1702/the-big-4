setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

source("helper_functions.R")

library(igraph)
library(ggplot2)
library(gridExtra)
library(coda)

compute_model_metrics_table <- function(dataset) {
     
     # ------------------------------------------------------------------
     # Load data
     # Expects files named: data_<dataset>.RData
     # ------------------------------------------------------------------
     data_file <- paste0("data_", dataset, ".RData")
     if (!file.exists(data_file)) {
          stop("Data file not found: ", data_file)
     }
     
     load(data_file)
     
     # ------------------------------------------------------------------
     # Model 1: SMN
     # Expects: samples_model_1_<dataset>.RData
     # and C++ / helper functions already available or in same folder
     # ------------------------------------------------------------------
     Rcpp::sourceCpp("model_1_mcmc_algorithm.cpp")
     load(paste0("samples_model_1_", dataset, ".RData"))
     
     metrics <- prob_metrics_per_layer_model_1(samples, y)
     mean_metrics <- c(
          mean(metrics$auc$mean),
          mean(metrics$brier$mean),
          mean(metrics$logloss$mean)
     )
     IC <- c(
          dic_model_1(log_lik, samples, y)$DIC,
          waic_model_1(samples, y)$WAIC
     )
     metrics_summary_1 <- c(mean_metrics, IC)
     
     # ------------------------------------------------------------------
     # Model 2: SMN-C
     # ------------------------------------------------------------------
     Rcpp::sourceCpp("model_2_mcmc_algorithm.cpp")
     load(paste0("samples_model_2_", dataset, ".RData"))
     
     metrics <- prob_metrics_per_layer_model_2(samples, X, y)
     mean_metrics <- c(
          mean(metrics$auc$mean),
          mean(metrics$brier$mean),
          mean(metrics$logloss$mean)
     )
     IC <- c(
          dic_model_2(log_lik, samples, y, X)$DIC,
          waic_model_2(samples, y, X)$WAIC
     )
     metrics_summary_2 <- c(mean_metrics, IC)
     
     # ------------------------------------------------------------------
     # Model 3: SMN-C-BG
     # ------------------------------------------------------------------
     Rcpp::sourceCpp("model_3_mcmc_algorithm.cpp")
     load(paste0("samples_model_3_", dataset, ".RData"))
     
     metrics <- prob_metrics_per_layer_model_3(samples, X, y)
     mean_metrics <- c(
          mean(metrics$auc$mean),
          mean(metrics$brier$mean),
          mean(metrics$logloss$mean)
     )
     IC <- c(
          dic_model_3(log_lik, samples, y, X)$DIC,
          waic_model_3(samples, y, X)$WAIC
     )
     metrics_summary_3 <- c(mean_metrics, IC)
     
     # ------------------------------------------------------------------
     # Model 4: SMN-C-LD
     # ------------------------------------------------------------------
     Rcpp::sourceCpp("model_4_mcmc_algorithm.cpp")
     load(paste0("samples_model_4_", dataset, ".RData"))
     
     metrics <- prob_metrics_per_layer_model_4(samples, X, y)
     mean_metrics <- c(
          mean(metrics$auc$mean),
          mean(metrics$brier$mean),
          mean(metrics$logloss$mean)
     )
     IC <- c(
          dic_model_4(log_lik, samples, y, X)$DIC,
          waic_model_4(samples, y, X)$WAIC
     )
     metrics_summary_4 <- c(mean_metrics, IC)
     
     # ------------------------------------------------------------------
     # Model 5: SMN-C-SB
     # ------------------------------------------------------------------
     Rcpp::sourceCpp("model_5_mcmc_algorithm.cpp")
     load(paste0("samples_model_5_", dataset, ".RData"))
     
     metrics <- prob_metrics_per_layer_model_5(samples, X, y)
     mean_metrics <- c(
          mean(metrics$auc$mean),
          mean(metrics$brier$mean),
          mean(metrics$logloss$mean)
     )
     IC <- c(
          dic_model_5(log_lik, samples, y, X)$DIC,
          waic_model_5(samples, y, X)$WAIC
     )
     metrics_summary_5 <- c(mean_metrics, IC)
     
     # ------------------------------------------------------------------
     # Bind and format output
     # ------------------------------------------------------------------
     out <- rbind(
          metrics_summary_1,
          metrics_summary_2,
          metrics_summary_3,
          metrics_summary_4,
          metrics_summary_5
     )
     
     colnames(out) <- c("auc", "brier", "logloss", "dic", "waic")
     rownames(out) <- c("SMN", "SMN-C", "SMN-C-BG", "SMN-C-LD", "SMN-C-SB")
     
     round(out, 3)
}

# Results ======================================================================

compute_model_metrics_table("metallica")
compute_model_metrics_table("slayer")
compute_model_metrics_table("megadeth")
compute_model_metrics_table("anthrax")

# End ==========================================================================