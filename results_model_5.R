# Setup ========================================================================

setwd("~/Dropbox/PAPERS/projects/MFDA")

# Clear workspace
rm(list = ls())

# Inference on zeta ============================================================

# Bands and file pattern
bands <- c("metallica", "slayer", "megadeth", "anthrax")
files <- paste0("samples_model_5_", bands, ".RData")

# Extract zeta samples for each band
zeta_list <- lapply(files, function(f) {
     load(f)
     samples$zeta
})
names(zeta_list) <- bands

# Posterior means and 95% credible intervals
zeta_mean  <- round(sapply(zeta_list, mean), 3)
zeta_lower <- round(sapply(zeta_list, quantile, probs = 0.025), 3)
zeta_upper <- round(sapply(zeta_list, quantile, probs = 0.975), 3)

zeta_summary <- data.frame(
     Band  = toupper(bands),
     Mean  = zeta_mean,
     L95   = zeta_lower,
     U95   = zeta_upper,
     row.names = NULL
)

zeta_summary

# Densities with common smoothing
dens_list <- lapply(zeta_list, density, adjust = 1.5, n = 5000)

# Common x- and y-ranges
range_x <- range(unlist(zeta_list))
range_y <- range(unlist(lapply(dens_list, function(d) d$y)))

# Colors and labels
cols <- c("black", "red", "blue", "darkgreen")
labels <- toupper(bands)

# Plot
pdf(file = paste0("model_5_posterior_zeta.pdf"), width = 7, height = 7, pointsize = 20)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

plot(
     dens_list[[1]],
     col  = cols[1],
     lwd  = 2,
     xlim = range_x,
     ylim = range_y,
     xlab = expression(zeta),
     ylab = "Density",
     main = "")

for (i in 2:length(dens_list))
     lines(dens_list[[i]], col = cols[i], lwd = 2)

legend("topleft", legend = labels, fill = cols, border = cols,  bty = "n")

dev.off()

# Inference on mu ==============================================================

plot_mu_posterior_k <- function(k, ylim , plot_legend = F) {
     # Bands and file names
     bands <- c("metallica", "slayer", "megadeth", "anthrax")
     files <- paste0("samples_model_5_", bands, ".RData")
     nb    <- length(bands)
     
     zeta_list <- vector("list", nb)
     mu_list   <- vector("list", nb)
     
     # Load samples and extract zeta and mu for each band
     for (b in seq_along(files)) {
          load(files[b])
          zeta_list[[b]] <- samples$zeta
          mu_list[[b]]   <- samples$mu
     }
     
     # Number of layers
     K <- ncol(mu_list[[1]])
     
     # Add global intercept zeta to each column of mu (all layers)
     for (b in seq_along(mu_list))
          mu_list[[b]] <- sweep(mu_list[[b]], 1, zeta_list[[b]], "+")
     
     # Posterior summaries for layer k
     mu_mean  <- sapply(mu_list, function(m) mean(m[, k]))
     mu_lower <- sapply(mu_list, function(m) quantile(m[, k], probs = 0.025))
     mu_upper <- sapply(mu_list, function(m) quantile(m[, k], probs = 0.975))
     
     mu_summary <- data.frame(
          Band = toupper(bands),
          Mean = mu_mean,
          L95  = mu_lower,
          U95  = mu_upper,
          row.names = NULL
     )
     
     # Densities with common smoothing for layer k 
     dens_list <- lapply(mu_list, function(m) density(m[, k], adjust = 1.5, n = 5000))
     
     range_x <- range(unlist(lapply(mu_list, function(m) m[, k])))
     
     # Colors and labels
     cols   <- c("black", "red", "blue", "darkgreen")
     labels <- c("MET", "SLA", "MEG", "ANT")
     
     # Plot to PDF 
     pdf(file = sprintf("model_5_posterior_mu_%d.pdf", k), width = 5, height = 5, pointsize = 20)
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     
     plot(
          dens_list[[1]],
          col  = cols[1],
          lwd  = 2,
          xlim = range_x,
          ylim = ylim,
          xlab = bquote(zeta + mu[.(k)]),
          ylab = "Density",
          main = ""
     )
     
     for (b in 2:nb)
          lines(dens_list[[b]], col = cols[b], lwd = 2)
     
     if(plot_legend) 
          legend("topright", legend = labels, fill = cols, border = cols, ncol = 1, bty = "n")
     
     dev.off()
     
     # Return summary table
     mu_summary
}

plot_mu_posterior_k(k = 1, ylim = c(0, 1.2), TRUE)
plot_mu_posterior_k(k = 2, ylim = c(0, 1.2))
plot_mu_posterior_k(k = 3, ylim = c(0, 1.2))
plot_mu_posterior_k(k = 4, ylim = c(0, 1.2))

# Variance components ==========================================================

plot_param_posterior <- function(param, range_x, range_y, xlab_expr, plot_legend = F) {
     # Settings
     model <- "model_5"
     bands <- c("metallica", "slayer", "megadeth", "anthrax")
     out_dir <- "."
     
     # 1Build file names and extract parameter samples for each band
     files <- file.path(out_dir, paste0("samples_", model, "_", bands, ".RData"))
     
     param_list <- lapply(files, function(f) {
          env <- new.env()
          load(f, envir = env)
          if (!("samples" %in% ls(env))) {
               stop("Object 'samples' not found in file: ", f)
          }
          if (!(param %in% names(env$samples))) {
               stop("Parameter '", param, "' not found in samples in file: ", f)
          }
          env$samples[[param]]
     })
     names(param_list) <- bands
     
     # Posterior means and 95% credible intervals
     param_mean  <- round(sapply(param_list, mean), 3)
     param_lower <- round(sapply(param_list, quantile, probs = 0.025), 3)
     param_upper <- round(sapply(param_list, quantile, probs = 0.975), 3)
     
     summary_df <- data.frame(
          Band = toupper(bands),
          Mean = param_mean,
          L95  = param_lower,
          U95  = param_upper,
          row.names = NULL
     )
     
     # Densities with common smoothing
     dens_list <- lapply(param_list, density, adjust = 1.5, n = 5000)
     
     # Common x- and y-ranges
     if (missing(range_x)) range_x <- range(unlist(param_list))
     if (missing(range_y)) range_y <- range(unlist(lapply(dens_list, function(d) d$y)))
     
     # Colors, labels, and x-axis label as math expression
     cols <- c("black", "red", "blue", "darkgreen")
     labels <- c("MET", "SLA", "MEG", "ANT")
     
     # Try to interpret `param` as a math expression, fall back to plain text
     if (missing(xlab_expr)) {
          xlab_expr <- tryCatch(
               parse(text = param)[[1]],
               error = function(e) param
          )
     }
     
     # Plot to PDF
     pdf_file <- file.path(out_dir, paste0(model, "_posterior_", param, ".pdf"))
     pdf(file = pdf_file, width = 5, height = 5, pointsize = 22)
     on.exit(dev.off(), add = TRUE)
     
     op <- par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     on.exit(par(op), add = TRUE)
     
     plot(
          dens_list[[1]],
          col  = cols[1],
          lwd  = 2,
          xlim = range_x,
          ylim = range_y,
          xlab = xlab_expr,
          ylab = "Density",
          main = ""
     )
     
     if (length(dens_list) > 1L) {
          for (i in 2:length(dens_list)) {
               lines(dens_list[[i]], col = cols[i], lwd = 2)
          }
     }
     
     if (plot_legend) 
          legend("topright", legend = labels, fill = cols, border = cols, bty = "n")
     
     # Return summary table
     return(summary_df)
}

plot_param_posterior("omega2", range_x = c(0, 10),  xlab_expr = expression(omega^2), plot_legend = T)
plot_param_posterior("sigma2", range_x = c(0, 3),   xlab_expr = expression(sigma^2))
plot_param_posterior("tau2",   range_x = c(0, 0.2), xlab_expr = expression(tau^2))
plot_param_posterior("kappa2", range_x = c(0, 0.3), xlab_expr = expression(kappa^2))
plot_param_posterior("rho2",   range_x = c(0, 15),  xlab_expr = expression(rho^2))

# Regression coefficients ======================================================

plot_beta_posterior <- function(layer, band, range_y, coef_names) {
     # Settings
     model <- "model_5"
     out_dir <- "."
     
     # Build file name for model 5
     rdata_file <- paste0("samples_model_5_", band, ".RData")
     
     # Load 
     load(rdata_file)
     
     B <- dim(samples$beta)[1L]
     p <- dim(samples$beta)[2L]
     K <- dim(samples$beta)[3L]
     
     # Extract samples for the chosen layer: n_iter x p
     beta_layer <- samples$beta[ , , layer, drop = FALSE]
     beta_layer <- beta_layer[ , , 1L]
     
     # Coefficient names (if available)
     if (missing(coef_names)) coef_names <- paste0("beta", 1:p)
     
     # Posterior means and 95% credible intervals
     beta_mean  <- apply(beta_layer, 2L, mean)
     beta_l     <- apply(beta_layer, 2L, quantile, probs = 0.05)
     beta_u     <- apply(beta_layer, 2L, quantile, probs = 0.95)
     
     # Sign classification for coloring
     sign_class <- ifelse(beta_l > 0, "positive",
                          ifelse(beta_u < 0, "negative", "includes_zero"))
     cols <- ifelse(sign_class == "positive", "darkgreen",
                    ifelse(sign_class == "negative", "red", "black"))
     
     # x positions for regressors
     x_pos <- seq_len(p)
     
     # y-limits with a bit of padding
     if (missing(range_y)) {
          y_min <- min(beta_l)
          y_max <- max(beta_u)
          y_pad <- 0.05 * (y_max - y_min)
          range_y <- c(y_min - y_pad, y_max + y_pad)
     }
     
     # Base plot
     pdf_file <- file.path(out_dir, paste0(model, "_posterior_betas_", band, "_", layer, ".pdf"))
     pdf(file = pdf_file, width = 5, height = 5, pointsize = 23)
     on.exit(dev.off(), add = TRUE)
     
     op <- par(mar = c(2.55, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     on.exit(par(op), add = TRUE)
     
     plot(
          x_pos, beta_mean,
          ylim = range_y,
          xlab = "",
          ylab = expression(beta),
          xaxt = "n",
          pch  = 19,
          col  = cols,
          main = "",
          cex = 1,
          cex.axis = 0.7
     )
     
     # Horizontal line at 0
     abline(h = 0, lwd = 2, col = "gray")
     
     # Credible interval segments
     segments(x0 = x_pos, y0 = beta_l,
              x1 = x_pos, y1 = beta_u,
              col = cols, lwd = 1)
     
     # X-axis with regressor names
     axis(1, at = x_pos, labels = coef_names, las = 2)
     
     # Summary data frame to return
     beta_summary <- data.frame(
          Band   = toupper(band),
          Layer  = layer,
          Coef   = coef_names,
          Mean   = beta_mean,
          L95    = beta_l,
          U95    = beta_u,
          Sign   = sign_class,
          row.names = NULL
     )
     
     return(beta_summary)
}

coef_names <- c("year", "bmp", "dur", "alb", "emo", "vad")

plot_beta_posterior(layer = 1, band = "metallica", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 2, band = "metallica", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 3, band = "metallica", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 4, band = "metallica", range_y = c(-1,1)*0.55, coef_names = coef_names)

plot_beta_posterior(layer = 1, band = "slayer", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 2, band = "slayer", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 3, band = "slayer", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 4, band = "slayer", range_y = c(-1,1)*0.55, coef_names = coef_names)

plot_beta_posterior(layer = 1, band = "megadeth", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 2, band = "megadeth", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 3, band = "megadeth", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 4, band = "megadeth", range_y = c(-1,1)*0.55, coef_names = coef_names)

plot_beta_posterior(layer = 1, band = "anthrax", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 2, band = "anthrax", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 3, band = "anthrax", range_y = c(-1,1)*0.55, coef_names = coef_names)
plot_beta_posterior(layer = 4, band = "anthrax", range_y = c(-1,1)*0.55, coef_names = coef_names)

# Songs with significant effects ===============================================

summarize_vartheta <- function(catalog,
                               samples,
                               ci = 0.95,
                               digits = 3,
                               filter_positive = TRUE,
                               sort_decreasing = TRUE) {
     # vartheta samples
     vartheta_draws <- samples$vartheta
     
     # Posterior means
     vartheta_means <- colMeans(vartheta_draws, na.rm = TRUE)
     
     # Credible intervals
     alpha   <- (1 - ci) / 2
     probs   <- c(alpha, 1 - alpha)
     vartheta_ci <- t(apply(vartheta_draws, 2, stats::quantile, probs = probs, na.rm = TRUE))
     
     # Build data.frame
     df <- data.frame(
          Song  = as.character(catalog$song),
          Album = as.character(catalog$album),
          Mean  = round(vartheta_means, digits),
          Lower = round(vartheta_ci[, 1], digits),
          Upper = round(vartheta_ci[, 2], digits),
          stringsAsFactors = FALSE
     )
     
     # Sort by Mean
     o  <- order(df$Mean, decreasing = sort_decreasing, na.last = TRUE)
     df <- df[o, , drop = FALSE]
     row.names(df) <- NULL
     
     # Optional filter: keep only “significantly positive” baselines
     if (filter_positive) {
          df <- df[df$Lower > 0, , drop = FALSE]
          row.names(df) <- NULL
     }
     
     df
}

summarize_vartheta_band <- function(band) {
     # Build file names
     sample_file <- paste0("samples_model_5_", band, ".RData")
     data_file   <- paste0("data_", band, ".RData")
     
     # Load samples and data
     load(sample_file)  # should create `samples`
     load(data_file)    # should create `res_<band>`
     
     # Get the corresponding `res_*` object
     res_obj_name <- paste0("res_", band)
     res_obj      <- get(res_obj_name)
     
     # Extract catalog (logRMS layer, as in your example)
     catalog <- res_obj$logRMS$catalog
     
     # Call the existing helper
     out <- summarize_vartheta(catalog, samples)
     
     return(out)
}

summarize_vartheta_band("metallica")[1:10, ]
summarize_vartheta_band("slayer")
summarize_vartheta_band("megadeth") [1:10, ]
summarize_vartheta_band("anthrax")

# Community detection ==========================================================

get_dahl_partition <- function(xi_mat) {
     S <- nrow(xi_mat)
     n <- ncol(xi_mat)
     
     # 1. Posterior similarity matrix (PSM)
     psm <- matrix(0, n, n)
     for (s in 1:S) {
          z <- xi_mat[s, ]
          psm <- psm + outer(z, z, "==")
     }
     psm <- psm / S
     
     # 2. Dahl's least-squares partition
     best_s <- 1L
     best_d <- Inf
     
     for (s in 1:S) {
          z   <- xi_mat[s, ]
          A_s <- outer(z, z, "==")
          d_s <- sum((A_s - psm)^2)
          if (d_s < best_d) {
               best_d <- d_s
               best_s <- s
          }
     }
     
     list(
          cluster = xi_mat[best_s, ],  
          psm     = psm
     )
}

plot_band_posterior_networks <- function(band,
                                         album_col    = "album",
                                         out_prefix   = "model_5_network_",
                                         layer_names  = c("logRMS", "logSC", "logitSFM", "logFlux"),
                                         layer_labels = c("loudness", "brightness", "tonality", "rhythm"),
                                         thin_xi      = 5) {
     
     # Required packages
     suppressMessages(suppressWarnings(library(igraph)))
     suppressMessages(suppressWarnings(library(mclust)))
     suppressMessages(suppressWarnings(library(aricode)))
     
     ##  Load band-specific objects
     samples_file <- paste0("samples_model_5_", band, ".RData")
     data_file    <- paste0("data_", band, ".RData")
     
     if (!file.exists(samples_file)) {
          stop("File not found: ", samples_file)
     }
     if (!file.exists(data_file)) {
          stop("File not found: ", data_file)
     }
     
     load(samples_file)  # provides 'samples'
     load(data_file)     # provides 'res_<band>'
     
     res_name <- paste0("res_", band)
     if (!exists(res_name)) {
          stop("Object ", res_name, " not found in the workspace after loading ", data_file)
     }
     res_obj <- get(res_name)
     
     ##  Dimensions
     dim_delta <- dim(samples$delta)
     dim_xi    <- dim(samples$xi)
     
     S_tot <- dim_delta[1]
     n     <- dim_delta[2]
     K     <- dim_delta[3]
     
     ##  Posterior means of delta_{i,k} (n x K)
     delta_mean <- apply(samples$delta, c(2, 3), mean)
     
     ##  Posterior Dahl partitions for xi_{i,k} (n x K)
     xi_hat   <- matrix(NA_integer_, nrow = n, ncol = K)
     psm_list <- vector("list", K)
     
     for (k in 1:K) {
          # xi for layer k: S_tot x n
          xi_k <- samples$xi[ , , k]
          
          # Thinning if S_tot is large
          idx  <- seq(1, S_tot, by = thin_xi)
          xi_k <- xi_k[idx, , drop = FALSE]
          
          dahl_res      <- get_dahl_partition(xi_k)
          xi_hat[, k]   <- dahl_res$cluster
          psm_list[[k]] <- dahl_res$psm
     }
     
     ##  Album-based partition (same across layers)
     cat_obj   <- res_obj[[layer_names[1]]]$catalog
     album_vec <- as.factor(cat_obj[[album_col]])
     album_cl  <- as.integer(album_vec)
     
     ##  Loop over layers: plots + ARI
     ari_vec <- numeric(K)
     
     for (k in 1:K) {
          layer_name  <- layer_names[k]
          layer_label <- layer_labels[k]
          
          # Graph for this layer
          g <- res_obj[[layer_name]]$graph
          
          # Layout
          set.seed(42)
          layout_mat <- layout_with_kk(g)
          
          # Node sizes from delta_mean[, k] scaled to a reasonable range
          d_raw    <- delta_mean[, k]
          size_min <- 3
          size_max <- 9
          
          if (diff(range(d_raw)) == 0) {
               vsize <- rep((size_min + size_max) / 2, n)
          } else {
               vsize <- size_min + (d_raw - min(d_raw)) / diff(range(d_raw)) * (size_max - size_min)
          }
          
          # Node colors from communities xi_hat[, k]
          comm_k   <- xi_hat[, k]
          n_clust  <- max(comm_k, na.rm = TRUE)
          comm_cols <- rainbow(n_clust)
          vcol      <- comm_cols[comm_k]
          
          # Plot to PDF
          pdf(
               file      = paste0(out_prefix, band, "_", layer_label, ".pdf"),
               width     = 5,
               height    = 5,
               pointsize = 20
          )
          par(mar = 0 * c(2, 2, 2, 2))
          
          plot(
               g,
               layout             = layout_mat,
               vertex.size        = vsize,
               vertex.color       = grDevices::adjustcolor(vcol, 0.9),
               vertex.frame.color = grDevices::adjustcolor(vcol, 0.9),
               vertex.label       = NA,
               edge.color         = grDevices::adjustcolor("black", 0.3),
               main               = ""
          )
          
          dev.off()
          
          # ARI between community partition and album partition
          ari_vec[k] <- ARI(comm_k, album_cl)
     }
     
     ##  Return ARI as data frame
     ari_out <- data.frame(
          layer = layer_labels,
          ARI   = round(ari_vec, 3),
          row.names = NULL
     )
     
     ari_out
}

plot_band_posterior_networks("metallica")
plot_band_posterior_networks("slayer")
plot_band_posterior_networks("megadeth")
plot_band_posterior_networks("anthrax")

# End ==========================================================================