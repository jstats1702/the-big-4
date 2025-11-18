setwd("~/Dropbox/PAPERS/projects/MFDA")

rm(list = ls())

# variar k en 3,5, 10
# variar la distancia (cor, cos, y las alternativas de la funcion dist)

# Helper functions -------------------------------------------------------------

build_song_similarity_network <- function(years,
                                          albums,
                                          bands,
                                          k = 3,
                                          M = 1000,
                                          metrics = c("logRMS", "logSC", "logitSFM", "logFlux"),
                                          data_dir) {
     # Initialize a catalog of songs
     catalog <- data.frame(
          year  = integer(0),
          band  = character(0),
          album = character(0),
          song  = character(0),
          stringsAsFactors = FALSE
     )
     
     # Initialize a list to store curves for each metric
     curves_by_metric <- setNames(vector("list", length(metrics)), metrics)
     for (m in metrics) curves_by_metric[[m]] <- list()
     
     # Load each album file and extract audio features for all songs
     for (idx in seq_along(years)) {
          year  <- years[idx]
          band  <- bands[idx]
          album <- albums[idx]
          file_path <- file.path(data_dir, sprintf("%s_%s_%s.RData", year, band, album))
          load(file_path)
          
          for (i in seq_along(dat)) {
               af <- dat[[i]]$audio_features
               
               # Add song metadata to catalog
               catalog <- rbind(
                    catalog,
                    data.frame(year = year, band = band, album = album, song = dat[[i]]$song,
                               stringsAsFactors = FALSE)
               )
               
               # Collect the smoothed + z-scored curve for each metric (resampled to common grid)
               for (m in metrics) {
                    y <- stats::approx(
                         x    = af$time,
                         y    = af[[paste0(m, "_smooth_z")]],
                         xout = seq(0, 1, length.out = M),
                         rule = 2
                    )$y
                    curves_by_metric[[m]][[length(curves_by_metric[[m]]) + 1]] <- as.numeric(y)
               }
          }
     }
     
     # Number of total songs
     rownames(catalog) <- NULL
     n <- nrow(catalog)
     
     # Build similarity networks for each metric
     out <- list()
     for (m in metrics) {
          # Gather all curves for metric m into a matrix (songs × timepoints)
          # Euclidean distances between songs
          D <- as.matrix(stats::dist(do.call(rbind, curves_by_metric[[m]]), method = "canberra"))
          D[lower.tri(D)] <- t(D)[lower.tri(D)]
          diag(D) <- 0
          
          # Correlation distance
          # X <- do.call(rbind, curves_by_metric[[m]])
          # S <- cor(t(X))
          # D <- 1 - S
          # diag(D) <- 0
          
          # Cosine distance
          # X <- do.call(rbind, curves_by_metric[[m]])
          # nrms <- sqrt(rowSums(X^2)) + 1e-12
          # C <- tcrossprod(X / nrms)
          # D <- 1 - C
          # diag(D) <- 0
          
          # k-NN pruning + symmetrization
          # Keep top-k per row (ties broken by order)
          W <- 1 / (D + 1e-12)
          A <- matrix(0L, n, n)
          for (i in seq_len(n)) {
               nz <- which(W[i, ] > 0)
               if (length(nz)) {
                    kk  <- min(k, length(nz))
                    ord <- order(W[i, nz], decreasing = TRUE)
                    top <- nz[ord][seq_len(kk)]
                    A[i, top] <- 1L
               }
          }
          A <- (pmax(A, t(A)) > 0L) * 1L
          diag(A) <- 0L
          
          # Build igraph object for easy network analysis/plotting
          g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
          igraph::V(g)$name  <- catalog$song
          igraph::V(g)$song  <- catalog$song
          igraph::V(g)$band  <- catalog$band
          igraph::V(g)$album <- catalog$album
          igraph::V(g)$year  <- catalog$year
          
          # Store results for this metric
          out[[m]] <- list(
               metric      = m,
               adjacency   = A,   # binary adjacency
               distances   = D,       # distance matrix
               # sigma       = sig,     # kernel bandwidth
               graph       = g,       # igraph object
               catalog     = catalog  # song metadata
          )
     }
     
     out
}

plot_album_network <- function(g, main_title, show_legend) {
     # Degree and album labels
     deg    <- igraph::degree(g)
     albums <- igraph::V(g)$album
     ua     <- unique(albums)
     
     # Colors per album (Set3 up to 12; fallback to rainbow if more)
     pal <- grDevices::rainbow(length(ua))
     album_colors  <- setNames(pal, ua)
     vertex_color  <- album_colors[albums]
     
     # Edge styling (neutral)
     edge_color <- grDevices::adjustcolor("black", alpha.f = 0.3)
     edge_width <- 1
     
     # Plot
     set.seed(42)
     plot(g,
          layout             = igraph::layout_with_kk,
          vertex.label       = NA,
          vertex.color       = grDevices::adjustcolor(vertex_color, 0.9),
          vertex.frame.color = grDevices::adjustcolor(vertex_color, 0.9),
          vertex.size        = 2.0 * sqrt(deg),
          edge.color         = edge_color,
          edge.width         = edge_width,
          main               = main_title)
     
     # Auto legend (albums; colored words, no lines)
     if (show_legend) {
          legend("topleft",
                 legend   = ua,
                 text.col = album_colors[ua],
                 cex = 0.9, text.font = 2,
                 lty = 0, lwd = 0, pch = NA, bty = "n",
                 x.intersp = 0, seg.len = 0)
     }
}

read_album_txt <- function(year, album, band, data_dir) {
     # Build filename
     f   <- file.path(data_dir, sprintf("%d_%s_%s.txt", year, band, album))
     # Read metadata CSV for the album (expects columns listed in 'req')
     df  <- utils::read.csv(f, stringsAsFactors = FALSE)
     # Keep only the required columns (song-level metadata)
     # Expected semantics:
     #   - song: unique song id/name
     #   - album, band, year: identifiers
     #   - duration: seconds (numeric)
     #   - bpm: beats per minute (numeric)
     req <- c("song", "album", "band", "year", "duration", "bpm")
     df[, req]
}

read_album_rdata <- function(year, album, band, data_dir) {
     f <- file.path(data_dir, sprintf("%d_%s_%s.RData", year, band, album))
     
     # Load 'dat' from the RData file
     e <- new.env(parent = emptyenv())
     load(f, envir = e)  # should create 'dat'
     dat <- e$dat
     
     # Helper to extract numeric values from lyrics_features; return NA if missing
     get_lf_num <- function(lf, nm) {
          v <- lf[[nm]]
          if (is.null(v)) NA_real_ else as.numeric(v)
     }
     
     rows <- lapply(dat, function(x) {
          lf <- x$lyrics_features
          data.frame(
               song  = if (!is.null(x$song))  x$song  else NA_character_,
               album = if (!is.null(x$album)) x$album else album,
               band  = if (!is.null(x$band))  x$band  else band,
               year  = if (!is.null(x$year))  x$year  else year,
               prop_anger        = get_lf_num(lf, "prop_anger"),
               prop_anticipation = get_lf_num(lf, "prop_anticipation"),
               prop_disgust      = get_lf_num(lf, "prop_disgust"),
               prop_fear         = get_lf_num(lf, "prop_fear"),
               prop_joy          = get_lf_num(lf, "prop_joy"),
               prop_sadness      = get_lf_num(lf, "prop_sadness"),
               prop_surprise     = get_lf_num(lf, "prop_surprise"),
               prop_trust        = get_lf_num(lf, "prop_trust"),
               mean_valence      = get_lf_num(lf, "mean_valence"),
               mean_arousal      = get_lf_num(lf, "mean_arousal"),
               mean_dominance    = get_lf_num(lf, "mean_dominance"),
               stringsAsFactors = FALSE
          )
     })
     
     do.call(rbind, rows)
}

standardize_upper <- function(M) {
     n  <- nrow(M)
     up <- upper.tri(M, diag = FALSE)
     v  <- M[up]
     m  <- mean(v, na.rm = TRUE)
     s  <- sd  (v, na.rm = TRUE)
     Z <- (M - m) / s
     diag(Z) <- 0
     Z[lower.tri(Z)] <- t(Z)[lower.tri(Z)] 
     Z
}

absdiff <- function(x) { 
     M <- abs(outer(x, x, "-"))
     diag(M) <- 0 
     M 
}

build_dyad_mats <- function(df) {
     # df: rows in the same order as the graph vertices
     n <- nrow(df)
     
     # Absolute differences
     YR  <- absdiff(df$year)      # year
     BPM <- absdiff(df$bpm)       # bpm
     DUR <- absdiff(df$duration)  # duration (sec)
     
     # Same album
     ALB <- outer(df$album, df$album, `==`) * 1.0
     diag(ALB) <- 0
     
     # emo_cos (NRC-8)
     E <- as.matrix(df[, c("prop_anger",
                           "prop_anticipation",
                           "prop_disgust",
                           "prop_fear","prop_joy",
                           "prop_sadness",
                           "prop_surprise",
                           "prop_trust")])
     # Normalize to unit norm (avoid division by zero)
     nr <- sqrt(rowSums(E^2))
     nr[nr == 0 | is.na(nr)] <- 1
     EU <- E / nr
     EMO <- EU %*% t(EU)
     diag(EMO) <- 0
     
     # VAD distance (valence, arousal, dominance) with column-wise z-score
     VAD <- scale(as.matrix(df[, c("mean_valence",
                                   "mean_arousal",
                                   "mean_dominance")]))
     # Euclidean distance
     VAD <- as.matrix(dist(VAD, method = "euclidean"))
     diag(VAD) <- 0
     
     # Standardize each slice over i<j (do not standardize ALB)
     YR  <- standardize_upper(YR)
     BPM <- standardize_upper(BPM)
     DUR <- standardize_upper(DUR)
     EMO <- standardize_upper(EMO)
     VAD <- standardize_upper(VAD)
     
     # Build the X cube (p=6) in order
     X <- array(0, dim = c(n, n, 6))
     X[, , 1] <- YR
     X[, , 2] <- BPM
     X[, , 3] <- DUR
     X[, , 4] <- ALB
     X[, , 5] <- EMO
     X[, , 6] <- VAD
     dimnames(X)[[3]] <- c("year","bpm","dur","alb","emo","vad")
     X
}

build_X <- function(g,
                    band,
                    year_vec,
                    album_vec,
                    data_dir = "~/Dropbox/PAPERS/projects/MFDA/data") {
     require(igraph)
     
     # Graph vertex order (key)
     song_g  <- igraph::V(g)$song
     album_g <- igraph::V(g)$album
     year_g  <- igraph::V(g)$year
     band_g  <- igraph::V(g)$band
     
     # Collect per-album
     pieces <- lapply(seq_along(year_vec), function(h) {
          year  <- year_vec[h]
          album <- album_vec[h]
          df_1  <- read_album_txt  (year, album, band, data_dir = data_dir)
          df_2  <- read_album_rdata(year, album, band, data_dir = data_dir)
          df    <- merge(df_1, df_2, by = c("song","album", "year"), all = TRUE)
          df
     })
     df <- do.call(rbind, pieces)
     df[is.na(df)] <- 0  # instrumental songs
     
     # Align to the graph order (and filter to songs present)
     df <- df[match(song_g, df$song), , drop = FALSE]
     
     # Build dyadic matrices and standardize
     X <- build_dyad_mats(df)
     
     list(
          X = X, 
          df = df,
          slice_names = dimnames(X)[[3]],
          p = dim(X)[3],
          n = nrow(df)
     )
}

# Metallica, Slayer, Megadeth, Anthrax (in order of founding) ------------------

year_data_metallica  <- c(
     1983,
     1984,
     1986,
     1987,
     1988,
     1991,
     1996,
     1997,
     1998,
     2003,
     2008,
     2016,
     2023
)

album_data_metallica <- c(
     "kill_em_all",
     "ride_the_lightning",
     "master_of_puppets",
     "garage_days_rerevisited",
     "and_justice_for_all",
     "metallica",
     "load",
     "reload",
     "garage_inc",
     "st_anger",
     "death_magnetic",
     "hardwired_to_self_destruct",
     "72_seasons"
)

year_data_slayer <- c(
     1983,
     1985,
     1986,
     1988,
     1990,
     1994,
     1998,
     2001,
     2006,
     2009,
     2015
)

album_data_slayer <- c(
     "show_no_mercy",
     "hell_awaits",
     "reign_in_blood",
     "south_of_heaven",
     "seasons_in_the_abyss",
     "divine_intervention",
     "diabolus_in_musica",
     "god_hates_us_all",
     "christ_illusion",
     "world_painted_blood",
     "repentless"
)

year_data_megadeth  <- c(
     1985,
     1986,
     1988,
     1990,
     1992,
     1994,
     1997,
     1999,
     2001,
     2004,
     2007,
     2009,
     2011,
     2013,
     2016,
     2022
)

album_data_megadeth <- c(
     "killing_is_my_business_and_business_is_good",
     "peace_sells_but_whos_buying",
     "so_far_so_good_so_what",
     "rust_in_peace",
     "countdown_to_extinction",
     "youthanasia",
     "cryptic_writings",
     "risk",
     "the_world_needs_a_hero",
     "the_system_has_failed",
     "united_abominations",
     "endgame",
     "th1rt3en",
     "super_collider",
     "dystopia",
     "the_sick_the_dying_and_the_dead"
)

year_data_anthrax  <- c(
     1984,
     1985,
     1987,
     1988,
     1990,
     1993,
     1995,
     1998,
     2003,
     2011,
     2016
)

album_data_anthrax <- c(
     "fistful_of_metal",
     "spreading_the_disease",
     "among_the_living",
     "state_of_euphoria",
     "persistence_of_time",
     "sound_of_white_noise",
     "stomp_442",
     "volume_8_the_threat_is_real",
     "weve_come_for_you_all",
     "worship_music",
     "for_all_kings"
)

# Network data -----------------------------------------------------------------

# Metallica
res_metallica <- build_song_similarity_network(
     years   = year_data_metallica,
     albums  = album_data_metallica,
     bands   = rep("metallica", length(year_data_metallica)),
     metrics = c("logRMS", "logSC", "logitSFM", "logFlux"),
     data_dir = "~/Dropbox/PAPERS/projects/MFDA/data" 
     # data_dir = "C:/Users/User/Dropbox/PAPERS/projects/MFDA/data/" 
)

# Slayer
res_slayer <- build_song_similarity_network(
     years   = year_data_slayer,
     albums  = album_data_slayer,
     bands   = rep("slayer", length(year_data_slayer)),
     metrics = c("logRMS", "logSC", "logitSFM", "logFlux"),
     data_dir = "~/Dropbox/PAPERS/projects/MFDA/data" 
     # data_dir = "C:/Users/User/Dropbox/PAPERS/projects/MFDA/data/" 
)

# Megadeth
res_megadeth <- build_song_similarity_network(
     years   = year_data_megadeth,
     albums  = album_data_megadeth,
     bands   = rep("megadeth", length(year_data_megadeth)),
     metrics = c("logRMS", "logSC", "logitSFM", "logFlux"),
     data_dir = "~/Dropbox/PAPERS/projects/MFDA/data" 
     # data_dir = "C:/Users/User/Dropbox/PAPERS/projects/MFDA/data/" 
)

# Anthrax
res_anthrax <- build_song_similarity_network(
     years   = year_data_anthrax,
     albums  = album_data_anthrax,
     bands   = rep("anthrax", length(year_data_anthrax)),
     metrics = c("logRMS", "logSC", "logitSFM", "logFlux"),
     data_dir = "~/Dropbox/PAPERS/projects/MFDA/data" 
     # data_dir = "C:/Users/User/Dropbox/PAPERS/projects/MFDA/data/" 
)

# Covariates data --------------------------------------------------------------

cov_metallica <- build_X(g         = res_metallica$logRMS$graph,
                         band      = "metallica",
                         year_vec  = year_data_metallica,
                         album_vec = album_data_metallica,
                         data_dir  = "~/Dropbox/PAPERS/projects/MFDA/data")

cov_slayer    <- build_X(g         = res_slayer$logRMS$graph,
                         band      = "slayer",
                         year_vec  = year_data_slayer,
                         album_vec = album_data_slayer,
                         data_dir  = "~/Dropbox/PAPERS/projects/MFDA/data")

cov_megadeth  <- build_X(g         = res_megadeth$logRMS$graph,
                         band      = "megadeth",
                         year_vec  = year_data_megadeth,
                         album_vec = album_data_megadeth,
                         data_dir  = "~/Dropbox/PAPERS/projects/MFDA/data")

cov_anthrax   <- build_X(g         = res_anthrax$logRMS$graph,
                         band      = "anthrax",
                         year_vec  = year_data_anthrax,
                         album_vec = album_data_anthrax,
                         data_dir  = "~/Dropbox/PAPERS/projects/MFDA/data")

# Save -------------------------------------------------------------------------

save(res_metallica, cov_metallica, file = "data_metallica.RData")
save(res_slayer,    cov_slayer,    file = "data_slayer.RData")
save(res_megadeth,  cov_megadeth,  file = "data_megadeth.RData")
save(res_anthrax,   cov_anthrax,   file = "data_anthrax.RData")

# Plots: Metallica -------------------------------------------------------------

pdf(file = paste0("network_metallica_loudness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_metallica$logRMS$graph,
     main_title = "",
     show_legend = F
)
dev.off()


pdf(file = paste0("network_metallica_brightness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_metallica$logSC$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_metallica_tonality.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_metallica$logitSFM$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_metallica_rhythm.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_metallica$logFlux$graph,
     main_title = "",
     show_legend = F
)
dev.off()

# Plots: Slayer ----------------------------------------------------------------

pdf(file = paste0("network_slayer_loudness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_slayer$logRMS$graph,
     main_title = "",
     show_legend = F
)
dev.off()


pdf(file = paste0("network_slayer_brightness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_slayer$logSC$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_slayer_tonality.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_slayer$logitSFM$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_slayer_rhythm.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_slayer$logFlux$graph,
     main_title = "",
     show_legend = F
)
dev.off()

# Plots: Megadeth --------------------------------------------------------------

pdf(file = paste0("network_megadeth_loudness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_megadeth$logRMS$graph,
     main_title = "",
     show_legend = F
)
dev.off()


pdf(file = paste0("network_megadeth_brightness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_megadeth$logSC$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_megadeth_tonality.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_megadeth$logitSFM$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_megadeth_rhythm.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_megadeth$logFlux$graph,
     main_title = "",
     show_legend = F
)
dev.off()

# Plots: Anthrax --------------------------------------------------------------

pdf(file = paste0("network_anthrax_loudness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_anthrax$logRMS$graph,
     main_title = "",
     show_legend = F
)
dev.off()


pdf(file = paste0("network_anthrax_brightness.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_anthrax$logSC$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_anthrax_tonality.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_anthrax$logitSFM$graph,
     main_title = "",
     show_legend = F
)
dev.off()

pdf(file = paste0("network_anthrax_rhythm.pdf"), width = 5, height = 5, pointsize = 20)
par(mar = 0*c(2, 2, 2, 2))
plot_album_network(
     g = res_anthrax$logFlux$graph,
     main_title = "",
     show_legend = F
)
dev.off()

# End --------------------------------------------------------------------------
