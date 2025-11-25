# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

source("helper_functions.R")

library(igraph)
library(ggplot2)
library(gridExtra)
library(coda)

# Helper functions =============================================================

ppc_rmse <- function(stat_name, test) {
     
     # Extract observed and simulated values for the chosen statistic
     obs_vec <- test$observed[[stat_name]]
     sim_mat <- test$sims[[stat_name]]
     
     # Number of layers
     K <- dim(sim_mat)[2L]
     
     # Point estimates per layer: mean of simulations
     est_vec <- colMeans(sim_mat, na.rm = TRUE)
     
     # RMSE across layers: sqrt( mean( (estimate - observed)^2 ) )
     rmse <- sqrt(mean((est_vec - obs_vec[seq_len(K)])^2, na.rm = TRUE))
     
     rmse
}


compute_ppc_stats <- function(dataset, model) {
     # Source the C++ code for the model
     cpp_file <- paste0(model, "_mcmc_algorithm.cpp")
     Rcpp::sourceCpp(cpp_file)
     
     # Load data
     data_file <- paste0("data_", dataset, ".RData")
     load(data_file)
     
     # Load MCMC samples
     samples_file <- paste0("samples_", model, "_", dataset, ".RData")
     load(samples_file)
     
     # Default layer titles from y
     K <- dim(y)[3L]
     layer_titles <- seq_len(K)
     
     # Call the appropriate plot_test_stats_<model>() function
     plot_fun_name <- paste0("plot_test_stats_", model)
     plot_fun <- get(plot_fun_name)
     
     if (model == "model_1") {
          # No covariates
          test <- plot_fun(samples, y, layer_titles = layer_titles)
     } else if (model %in% paste0("model_", 2:5)) {
          # Models 2–5 use covariates X
          test <- plot_fun(samples, y, X, layer_titles = layer_titles)
     } else {
          stop("Unknown model: ", model)
     }
     
     # Compute PPC coverage for each statistic
     stats <- c(
          "density",
          "transitivity",
          "assortativity",
          "mean_degree",
          "sd_degree",
          "mean_distance",
          "diameter"
     )
     
     out <- sapply(stats, ppc_rmse, test = test)
     out
}

compute_ppc_all_models <- function(dataset) {
     # Models 1 to 5
     models <- paste0("model_", 1:5)
     
     # Compute PPC stats for each model
     out_list <- lapply(models, 
                        function(m) {
                              compute_ppc_stats(dataset = dataset, model = m)
     })
     
     # Bind rows
     out <- do.call(rbind, out_list)
     
     # Set row names
     rownames(out) <- c("SMN", "SMN-C", "SMN-C-BG", "SMN-C-LD", "SMN-C-SB")
     
     round(out, 3)
}

# Results ======================================================================

compute_ppc_all_models("metallica")
compute_ppc_all_models("slayer")
compute_ppc_all_models("megadeth")
compute_ppc_all_models("anthrax")

# End ==========================================================================

