# Read audio file and always return as a Wave object
# file   : path to audio file (wav, mp3, au, flac, etc.)
# Output : tuneR::Wave using decoder defaults (original sample rate, native channels)
read_audio_as_wave <- function(file) {
     suppressMessages(suppressWarnings(library(tuneR)))
     suppressMessages(suppressWarnings(library(av)))
     
     # Robust input checks
     if (!file.exists(file)) stop("File not found: ", file)
     
     # Always decode via FFmpeg so defaults are consistent across formats
     tmp_wav <- tempfile(fileext = ".wav")
     on.exit(unlink(tmp_wav), add = TRUE)
     
     av_audio_convert(file, output = tmp_wav, channels = 2, sample_rate = NULL, verbose = F)
     readWave(tmp_wav)
}

# Smooth a metric curve using cubic splines
# y    : numeric vector (log-RMS values)
# spar : smoothing parameter
smooth_metric <- function(y, spar = 0.4) {
     # Fit/evaluate spline on the same index grid
     x   <- seq_along(y)
     fit <- stats::smooth.spline(x = x, y = y, spar = spar)
     stats::predict(fit, x = x)$y
}

# Center or scale a variable in a data frame
# df     : data frame containing the target variable
# var    : name of the variable to transform (string)
# method : "none" (default), "median", or "zscore"
center_apply <- function(df, var, method = c("none", "median", "zscore")) {
     method <- match.arg(method)
     if (method == "none") return(df)
     
     # Check variable exists
     if (!var %in% names(df)) stop("Variable '", var, "' not found in data frame.")
     
     # Compute transformation
     if (method == "median") {
          base <- median(df[[var]], na.rm = TRUE)
          df[[var]] <- df[[var]] - base
     }
     
     if (method == "zscore") {
          mu  <- mean(df[[var]], na.rm = TRUE)
          sdv <- sd(df[[var]],  na.rm = TRUE)
          if (!is.finite(sdv) || sdv == 0) sdv <- 1
          df[[var]] <- (df[[var]] - mu) / sdv
     }
     
     df
}

# Compute multiple audio metrics for an audio file
#     - logRMS     : proxy for loudness 
#                    higher = more energy, lower = quieter
#     - logSC      : proxy for spectral brightness 
#                    higher = brighter, lower = darker/warmer
#     - logitSFM   : spectral flatness in logit scale 
#                    higher = noisier/flatter, lower = more harmonic/peaked
#     - logFlux    : proxy for onset strength and rhythmic activity 
#                    higher = beats/accents, lower = sustained or silent
#
# file   : path to audio file (wav, mp3, au, etc.)
# win_ms : window size in ms for frame analysis
# hop_ms : hop size in ms (frame step)
# eps    : small constant to avoid log(0) and stabilize ratios
# spar   : smoothing parameter for cubic spline
#
# Output (data frame) columns:
#     - time             : normalized frame time in [0, 1]
#     - logRMS           : log(RMS + eps) per frame
#     - logRMS_smooth    : smoothed logRMS via cubic spline
#     - logSC            : log-frequency spectral centroid per frame
#     - logSC_smooth     : smoothed logSC via cubic spline
#     - logitSFM         : logit(spectral flatness)
#     - logitSFM_smooth  : smoothed logitSFM via cubic spline
#     - logFlux          : log(spectral flux + eps) per frame
#     - logFlux_smooth   : smoothed logFlux via cubic spline
audio_metrics_song <- function(file, win_ms = 46, hop_ms = 23, eps = 1e-8, spar = 0.4) {
     # Decode audio
     wav <- read_audio_as_wave(file)
     fs  <- wav@samp.rate
     
     # Mono + amplitude normalization
     x <- wav@left
     if (!is.null(wav@right)) x <- (wav@left + wav@right) / 2
     if (isTRUE(wav@pcm)) x <- x / (2^(wav@bit - 1))
     x <- as.numeric(x)
     
     # Frame geometry (in samples)
     win <- max(1L, round(win_ms * 1e-3 * fs))
     hop <- max(1L, round(hop_ms * 1e-3 * fs))
     if (length(x) < win) stop("Audio shorter than one analysis window.")
     nframes <- floor((length(x) - win) / hop) + 1L
     
     # Per-frame RMS
     rms <- numeric(nframes)
     for (i in seq_len(nframes)) {
          start <- (i - 1L) * hop + 1L
          frame <- x[start:(start + win - 1L)]
          rms[i] <- sqrt(mean(frame^2))
     }
     
     # Per-frame spectral measures (log spectral centroid, spectral flatness, spectral flux)
     half_len <- floor(win / 2L)
     hann <- 0.5 - 0.5 * cos(2 * pi * (0:(win - 1L)) / (win - 1L))
     k <- 1:half_len                 # exclude DC (k = 0)
     f <- (k / win) * fs             # Hz > 0
     
     lsc_vals   <- numeric(nframes)  # E[log f] with magnitude weights
     logit_sfm  <- numeric(nframes)  # logit of spectral flatness
     flux_vals  <- numeric(nframes)  # spectral flux (positive differences)
     prevM <- NULL
     
     for (i in seq_len(nframes)) {
          start   <- (i - 1L) * hop + 1L
          frame_w <- x[start:(start + win - 1L)] * hann
          
          X <- stats::fft(frame_w)
          M <- Mod(X)[k + 1L]  # one-sided, DC excluded via k
          
          # Log-frequency spectral centroid
          lsc_vals[i] <- sum(log(f) * M) / (sum(M) + eps)
          
          # Spectral flatness (geometric mean / arithmetic mean), then logit
          gm <- exp(mean(log(M + eps)))
          am <- mean(M) + eps
          sfm <- gm / am
          p <- min(max(sfm, eps), 1 - eps)  # clip for logit stability
          logit_sfm[i] <- log(p) - log(1 - p)
          
          # Spectral flux (half-wave rectified positive spectral change)
          if (is.null(prevM)) {
               flux_vals[i] <- 0
          } else {
               d <- M - prevM
               flux_vals[i] <- sum(pmax(d, 0))
          }
          prevM <- M
     }
     
     # Smoothing
     log_rms          <- log(rms + eps)
     log_rms_smooth   <- smooth_metric(log_rms,   spar = spar)
     log_sc_smooth    <- smooth_metric(lsc_vals,  spar = spar)
     logit_sfm_smooth <- smooth_metric(logit_sfm, spar = spar)
     log_flux         <- log(flux_vals + eps)
     log_flux_smooth  <- smooth_metric(log_flux,  spar = spar)
     
     # Original (untransformed) metrics
     out <- data.frame(
          time            = if (nframes > 1L) seq(0, 1, length.out = nframes) else 0,
          logRMS          = log_rms,
          logRMS_smooth   = log_rms_smooth,
          logSC           = lsc_vals,
          logSC_smooth    = log_sc_smooth,
          logitSFM        = logit_sfm,
          logitSFM_smooth = logit_sfm_smooth,
          logFlux         = log_flux,
          logFlux_smooth  = log_flux_smooth,
          stringsAsFactors = FALSE
     )
     
     # Append z-score transformed counterparts
     vars <- c("logRMS","logRMS_smooth","logSC","logSC_smooth",
               "logitSFM","logitSFM_smooth","logFlux","logFlux_smooth")
     
     for (v in vars) {
          z <- scale(out[[v]], center = TRUE, scale = TRUE)[, 1]
          out[[paste0(v, "_z")]] <- as.numeric(z)
     }
     
     out
}


# Read and tokenize lyrics from a text file, returning both line-level text
# and token-level data.
# file_path : path to lyrics text file
# stop_tbl  : tibble of stopwords to remove
read_and_tokenize_lyrics <- function(file_path, stop_tbl = tidytext::stop_words) {
     # Dependencies (quiet)
     suppressMessages({
          suppressWarnings({
               library(readr)
               library(dplyr)
               library(stringr)
               library(tidytext)
               library(tibble)
          })
     })
     
     # Read file and build line-wise tibble
     raw_lines <- read_lines(file_path)
     lyrics    <- tibble(line = seq_along(raw_lines), text = raw_lines)
     
     # Tokenize, remove stopwords, normalize accents
     # Map accented vowels to ASCII: áéíóú -> aeiou
     replacement_map <- c("á" = "a", "é" = "e", "í" = "i", "ó" = "o", "ú" = "u")
     
     tokens <- lyrics %>%
          unnest_tokens(input = text, output = word) %>%
          filter(!is.na(word)) %>%
          anti_join(stop_tbl, by = "word") %>%
          mutate(
               word = chartr(
                    old = paste(names(replacement_map), collapse = ""),
                    new = paste(unname(replacement_map), collapse = ""),
                    x   = word
               )
          )
     
     # Return
     list(lyrics = lyrics, tokens = tokens)
}

# Compute lyric-level text metrics from preprocessed inputs.
# lyrics_tbl : a tibble with columns (line, text), where each row is a lyric line
# tokens_tbl : a tibble with at least column (word) containing tokenized words
#
# Output (one-row tibble) columns:
#   coverage_bing     = proportion of tokens found in Bing lexicon
#   prop_negative     = proportion of all tokens labeled as negative (unconditional)
#   prop_anger        = share of tokens expressing anger (NRC)
#   prop_anticipation = share of tokens expressing anticipation (NRC)
#   prop_disgust      = share of tokens expressing disgust (NRC)
#   prop_fear         = share of tokens expressing fear (NRC)
#   prop_joy          = share of tokens expressing joy (NRC)
#   prop_sadness      = share of tokens expressing sadness (NRC)
#   prop_surprise     = share of tokens expressing surprise (NRC)
#   prop_trust        = share of tokens expressing trust (NRC)
#   mean_sent         = average sentiment score across lines (sentimentr)
#   sd_sent           = standard deviation of sentiment across lines
#   slope_sent        = linear trend of sentiment over line index
#   arc_range         = range (max – min) of sentiment across lines
#   ttr               = type–token ratio (unique words / total words)
#   hapax_prop        = proportion of words that appear exactly once
#   mean_valence      = average valence score (pleasantness, NRC VAD)
#   mean_arousal      = average arousal score (activation, NRC VAD)
#   mean_dominance    = average dominance score (control, NRC VAD)
#   coverage_vad      = proportion of tokens matched in NRC VAD lexicon
compute_lyric_metrics <- function(lyrics_tbl, tokens_tbl) {
     # Dependencies (quiet)
     suppressMessages({
          suppressWarnings({
               library(dplyr)
               library(stringr)
               library(tidytext)
               library(tibble)
               library(sentimentr)
          })
     })
     
     # Defensive guards
     if (nrow(tokens_tbl) == 0L) {
          out <- tibble(
               coverage_bing = 0, prop_negative = NA_real_,
               prop_anger = NA_real_, prop_anticipation = NA_real_, prop_disgust = NA_real_,
               prop_fear = NA_real_, prop_joy = NA_real_, prop_sadness = NA_real_,
               prop_surprise = NA_real_, prop_trust = NA_real_,
               mean_sent = NA_real_, sd_sent = NA_real_, slope_sent = NA_real_, arc_range = NA_real_,
               ttr = NA_real_, hapax_prop = NA_real_,
               mean_valence = NA_real_, mean_arousal = NA_real_, mean_dominance = NA_real_, coverage_vad = 0
          )
          
          return(out)
     }
     
     toks <- tokens_tbl |> mutate(w = str_to_lower(word))
     total_tokens <- nrow(toks)
     
     # Bing polarity: coverage + unconditional negativity
     bing_lex <- tidytext::get_sentiments("bing") |> select(word, sentiment)
     neg_vocab <- bing_lex |> filter(sentiment == "negative") |> pull(word)
     
     bing_stats <- toks |>
          mutate(
               in_lex = w %in% bing_lex$word,
               is_neg = w %in% neg_vocab
          ) |>
          summarise(
               coverage_bing = mean(in_lex),
               prop_negative = mean(is_neg)
          )
     
     # NRC emotions (basic 8 proportions, conditional on hits)
     # sentiment includes basic emotions + positive/negative 
     # many-to-many is expected: repeated words across tokens
     nrc <- tidytext::get_sentiments("nrc") 
     basic_emos <- c("anger","anticipation","disgust","fear","joy","sadness","surprise","trust")
     
     nrc_counts <- toks %>%
          inner_join(nrc, by = c("w" = "word"), relationship = "many-to-many") %>%
          filter(sentiment %in% basic_emos) %>%
          count(sentiment, name = "n")
     
     emo_props <- nrc_counts %>%
          dplyr::mutate(prop = n / sum(n)) %>%
          dplyr::select(sentiment, prop) %>%
          tidyr::complete(sentiment = basic_emos, fill = list(prop = 0)) %>%
          tidyr::pivot_wider(names_from = sentiment, values_from = prop) %>%
          dplyr::rename_with(~ paste0("prop_", .x), .cols = tidyselect::everything())
     
     # Sentiment arc (context-aware, negators/intensifiers)
     sline <- sentiment(lyrics_tbl$text) |>
          mutate(idx = dplyr::row_number())
     
     arc <- tibble(
          mean_sent  = mean(sline$sentiment),
          sd_sent    = stats::sd(sline$sentiment),
          slope_sent = coef(stats::lm(sentiment ~ idx, data = sline))[2],
          arc_range  = diff(range(sline$sentiment))
     )
     
     # Lexical richness
     lex <- toks |>
          summarise(
               n_tokens   = n(),
               n_types    = n_distinct(w),
               ttr        = n_types / n_tokens,
               hapax_prop = mean(!duplicated(w) & !duplicated(w, fromLast = TRUE)) # freq==1
          ) |>
          select(ttr, hapax_prop)
     
     # VAD (NRC VAD) via textdata
     # Expected columns: Word, Valence, Arousal, Dominance
     # Preload and normalize VAD lexicon
     vad_tbl <- textdata::lexicon_nrc_vad() %>%
          dplyr::transmute(
               word      = stringr::str_to_lower(Word),
               valence   = Valence,
               arousal   = Arousal,
               dominance = Dominance
          )
     
     # Join tokens to VAD lexicon 
     # many-to-many is expected: repeated words across tokens
     vad_joined <- toks %>%
          dplyr::inner_join(vad_tbl, by = c("w" = "word"), relationship = "many-to-many")
     
     if (nrow(vad_joined) > 0L) {
          vad <- vad_joined %>%
               dplyr::summarise(
                    mean_valence   = mean(valence),
                    mean_arousal   = mean(arousal),
                    mean_dominance = mean(dominance),
                    coverage_vad   = dplyr::n() / total_tokens
               )
     } else {
          vad <- tibble::tibble(
               mean_valence   = NA_real_,
               mean_arousal   = NA_real_,
               mean_dominance = NA_real_,
               coverage_vad   = 0
          )
     }
     
     # Return
     dplyr::bind_cols(bing_stats, emo_props, arc, lex, vad)
}

# Process one song and return a list with metadata + features
process_song <- function(band, album, song) {
     # Paths
     stub <- sprintf("%s_%s_%s", band, album, song)
     base <- file.path(band, stub)
     audio_file  <- paste0(base, ".mp3")
     lyrics_file <- paste0(base, ".txt")
     
     # Validate inputs
     stopifnot(is.character(band), nzchar(band))
     stopifnot(is.character(album), nzchar(album))
     stopifnot(is.character(song), nzchar(song))
     
     if (!file.exists(audio_file)) {
          stop(sprintf("Audio file not found: %s", audio_file))
     }
     if (!file.exists(lyrics_file)) {
          stop(sprintf("Lyrics file not found: %s", lyrics_file))
     }
     
     # Audio features
     audio_features <- audio_metrics_song(
          file    = audio_file,
          win_ms  = 46,
          hop_ms  = 23,
          eps     = 1e-8,
          spar    = 0.5
     )
     
     # Lyrics features
     lyrics_data <- read_and_tokenize_lyrics(lyrics_file)
     lyrics_features <- compute_lyric_metrics(
          lyrics_tbl = lyrics_data$lyrics,
          tokens_tbl = lyrics_data$tokens
     )
     
     # Output (keep structure and names)
     list(
          band            = band,
          album           = album,
          song            = song,
          audio_features  = audio_features,
          lyrics_data     = lyrics_data,
          lyrics_features = lyrics_features
     )
}

pick_C<- function(y, 
                  algo = c("louvain", "fast_greedy", "walktrap"), 
                  seed = 42) {
     stopifnot(length(dim(y)) == 3)
     set.seed(seed)
     K <- dim(y)[3]
     
     algo <- match.arg(algo)
     n_comms <- integer(K)
     
     for (k in seq_len(K)) {
          A <- y[, , k]
          g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
          
          if (igraph::gsize(g) == 0) {
               n_comms[k] <- igraph::vcount(g) 
               next
          }
          
          com <- switch(algo,
                        louvain     = igraph::cluster_louvain(g),
                        fast_greedy = igraph::cluster_fast_greedy(g),
                        walktrap    = igraph::cluster_walktrap(g)
          )
          n_comms[k] <- length(igraph::sizes(com))
     }
     
     list(C = max(n_comms), per_layer = n_comms)
}

get_xi <- function(y, seed = 42) {
     set.seed(seed)
     
     n <- dim(y)[1]
     K <- dim(y)[3]
     
     xi <- matrix(0L, nrow = n, ncol = K)
     
     for (k in seq_len(K)) {
          A <- y[, , k]
          g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
          
          com <- igraph::cluster_louvain(g)
          mem <- igraph::membership(com)  # 1-based, possibly non-contiguous
          xi[, k] <- as.integer(factor(mem, levels = sort(unique(mem)))) - 1L
     }
     
     storage.mode(xi) <- "integer"
     
     xi
}

mce <- function(x) {
     (sd(x, na.rm = T) / coda::effectiveSize(x)) / abs(mean(x, na.rm = T))
}

summ_vec <- function(x) {
     c(
          mean(x, na.rm = TRUE),
          sd(x, na.rm = TRUE),
          as.numeric(quantile(x, 0.025, na.rm = TRUE)),
          as.numeric(quantile(x, 0.975, na.rm = TRUE))
     )
}

summ_mat <- function(M, alpha = 0.05) {
     data.frame(
          Mean  = apply(M, 2, function(x) mean(x, na.rm = TRUE)),
          Lower = apply(M, 2, function(x) as.numeric(quantile(x, alpha/2,      na.rm = TRUE))),
          Upper = apply(M, 2, function(x) as.numeric(quantile(x, 1  - alpha/2, na.rm = TRUE)))
     )
}

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

layers_correlation <- function(delta_samples,
                               layer_titles = c("Loudness","Brightness","Tonality","Rhythm")) {
     # delta_samples: 3D array [iter, n, K] of delta draws
     B <- dim(delta_samples)[1]
     K <- dim(delta_samples)[3]
     
     # Correlation matrix per MCMC iteration
     R_list <- vector("list", B)
     for (iter in seq_len(B)) {
          delta_iter <- delta_samples[iter, , ]  # n x K
          R_list[[iter]] <- cor(delta_iter)
     }     
     # Stack to array [K x K x iter]
     R_array <- simplify2array(R_list)
     
     # Posterior summaries per (k,ℓ)
     mean_R  <- apply(R_array, c(1, 2), mean, na.rm = TRUE)
     sd_R    <- apply(R_array, c(1, 2), sd,   na.rm = TRUE)
     ci_all  <- apply(R_array, c(1, 2),
                      function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
     lower_R <- ci_all[1, , ]
     upper_R <- ci_all[2, , ]
     
     # Build tidy table for upper-triangular pairs
     idx <- t(combn(K, 2))
     pairs <- paste(layer_titles[idx[, 1]], "-", layer_titles[idx[, 2]])
     
     data.frame(
          Layers = pairs,
          Mean   = round(mean_R[idx],  3),
          SD     = round(sd_R[idx],    3),
          Q2.5   = round(lower_R[idx], 3),
          Q97.5  = round(upper_R[idx], 3),
          stringsAsFactors = FALSE
     )
}

plot_log_likelihood <- function(log_lik, 
                                point_size = 0.6,
                                alpha = 0.4) {
     
     log_lik_df <- data.frame(
          Iteration     = seq_along(log_lik),
          LogLikelihood = as.numeric(log_lik)
     )
     
     p <- ggplot2::ggplot(log_lik_df, ggplot2::aes(x = Iteration, y = LogLikelihood)) +
          ggplot2::geom_point(size = point_size, color = "black", alpha = alpha) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "", x = "Iteration", y = "Log-Likelihood") +
          ggplot2::theme(
               axis.title = ggplot2::element_text(size = 12),
               axis.text  = ggplot2::element_text(size = 10)
          )
     
     p
}

plot_mu_by_layer <- function(mu_samples,
                             layer_titles = c("Loudness", "Brightness", "Tonality", "Rhythm"),
                             ci = 0.95) {
     K <- ncol(mu_samples)
     
     # Layer titles
     if (is.null(layer_titles)) layer_titles <- paste0("Layer ", seq_len(K))
     stopifnot(length(layer_titles) == K)
     
     # Symmetric global y-limits
     q_lim <- range(mu_samples, na.rm = TRUE)
     y_max <- max(abs(q_lim))
     y_lim <- c(-y_max, y_max)
     
     # CI probs
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     
     # Summaries per layer
     means    <- apply(mu_samples, 2, mean, na.rm = TRUE)
     ci_lower <- apply(mu_samples, 2, quantile, probs = probs[1], na.rm = TRUE)
     ci_upper <- apply(mu_samples, 2, quantile, probs = probs[2], na.rm = TRUE)
     
     df <- data.frame(
          Layer = factor(layer_titles, levels = layer_titles),
          Mean  = means,
          Lower = ci_lower,
          Upper = ci_upper,
          stringsAsFactors = FALSE
     )
     
     # Color: green if CI > 0, red if CI < 0, black otherwise
     df$Color <- ifelse(df$Upper < 0, "red",
                        ifelse(df$Lower > 0, "darkgreen", "black"))
     
     p <- ggplot2::ggplot(df, ggplot2::aes(x = Layer, y = Mean, color = Color)) +
          ggplot2::geom_segment(ggplot2::aes(xend = Layer, y = Lower, yend = Upper), size = 0.7) +
          ggplot2::geom_point(size = 1.8) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
          ggplot2::scale_color_identity() +
          ggplot2::coord_cartesian(ylim = y_lim) +
          ggplot2::labs(title = "Layer intercepts",
                        x = "Layer",
                        y = expression(mu)) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               axis.text.x  = ggplot2::element_text(angle = 0, hjust = 0.5),
               axis.ticks.x = ggplot2::element_line()
          )
     
     print(p)
     invisible(list(plot = p, summary = df))
}

plot_delta_by_layer <- function(delta_samples,
                                layer_titles = c("Loudness", "Brightness", "Tonality", "Rhythm"),
                                ci = 0.95) {
     # delta_samples: 3D array [iter, n, K] of posterior samples for delta
     n <- dim(delta_samples)[2]
     K <- dim(delta_samples)[3]
     
     # Global symmetric Y limits using overall range
     q_lim <- range(delta_samples, na.rm = TRUE)
     y_max <- max(abs(q_lim))
     y_lim <- c(-y_max, y_max)
     
     alpha <- 1 - ci
     probs <- c(alpha/2, 1 - alpha/2)
     
     # Plot list
     plot_list <- vector("list", K)
     
     for (k in seq_len(K)) {
          # Matrix [iter, n] for layer k
          delta_k <- delta_samples[, , k, drop = FALSE][, , 1]
          
          # Per-individual summaries (no matrixStats)
          means    <- colMeans(delta_k, na.rm = TRUE)
          ci_lower <- apply(delta_k, 2, quantile, probs = probs[1], na.rm = TRUE)
          ci_upper <- apply(delta_k, 2, quantile, probs = probs[2], na.rm = TRUE)
          
          delta_df <- data.frame(
               ID    = seq_len(n),
               Mean  = means,
               Lower = ci_lower,
               Upper = ci_upper
          )
          
          # Order by posterior mean and create plotting order
          delta_df <- delta_df[order(delta_df$Mean), , drop = FALSE]
          delta_df$Order <- factor(seq_len(n))
          
          # Colors by whether the 95% CI excludes 0
          delta_df$Color <- ifelse(
               delta_df$Upper < 0, "red",
               ifelse(delta_df$Lower > 0, "darkgreen", "black")
          )
          
          # Plot for this layer
          p <- ggplot2::ggplot(delta_df, ggplot2::aes(x = Order, y = Mean, color = Color)) +
               ggplot2::geom_segment(ggplot2::aes(xend = Order, y = Lower, yend = Upper), size = 0.3) +
               ggplot2::geom_point(size = 0.3) +
               ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
               ggplot2::scale_color_identity() +
               ggplot2::coord_cartesian(ylim = y_lim) +
               ggplot2::labs(title = layer_titles[k], x = "Actors", y = expression(delta)) +
               ggplot2::theme_minimal(base_size = 11) +
               ggplot2::theme(
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    axis.text.x  = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank()
               )
          
          plot_list[[k]] <- p
     }
     
     # 2 x 2 (or auto-fit if K != 4)
     gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
     invisible(plot_list)
}

plot_vartheta <- function(vartheta_samples, ci = 0.95) {
     # vartheta_samples: matrix/array of posterior samples [iter, n]
     
     # Posterior means per component
     vartheta_means <- apply(vartheta_samples, 2, mean)
     
     # y range
     q_lim <- range(vartheta_samples, na.rm = TRUE)
     y_max <- max(abs(q_lim))
     y_lim <- c(-y_max, y_max)
     
     # Credible intervals per component
     alpha <- 1 - ci
     probs <- c(alpha/2, 1 - alpha/2)
     ci_95 <- t(apply(vartheta_samples, 2, stats::quantile, probs = probs))
     
     # Build ordered data frame
     df <- data.frame(
          Index       = seq_along(vartheta_means),
          Mean        = vartheta_means,
          CI_95_Lower = ci_95[, 1],
          CI_95_Upper = ci_95[, 2]
     )
     df <- df[order(df$Mean), , drop = FALSE]  # order by mean
     df$Rank  <- seq_len(nrow(df))             # ordered index
     
     # Color by whether the 95% CI is entirely above/below 0 or crosses 0
     df$Color <- ifelse(df$CI_95_Upper < 0, "red",
                        ifelse(df$CI_95_Lower > 0, "darkgreen", "black"))
     
     # Plot
     p <- ggplot2::ggplot(df) +
          ggplot2::geom_segment(
               ggplot2::aes(x = Rank, xend = Rank, y = CI_95_Lower, yend = CI_95_Upper, color = Color),
               linewidth = 0.6
          ) +
          ggplot2::geom_point(
               ggplot2::aes(x = Rank, y = Mean, color = Color),
               size = 0.6
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
          ggplot2::scale_x_continuous(breaks = df$Rank, labels = df$Index) +
          ggplot2::labs(title = "", y = NULL, x = "Actors") +
          ggplot2::coord_cartesian(ylim = y_lim) +
          ggplot2::scale_color_identity() +
          ggplot2::theme_minimal() +
          ggplot2::theme(
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               axis.text.x  = ggplot2::element_blank(),
               axis.ticks.x = ggplot2::element_blank()
          )
     
     print(p)
     invisible(p)
}

plot_beta_by_layer <- function(beta_samples,
                               covariate_labels = c("year", "bpm", "dur", "alb", "emo", "vad"),
                               layer_titles = c("Loudness", "Brightness", "Tonality", "Rhythm"),
                               ci = 0.95) {
     p <- dim(beta_samples)[2]
     K <- dim(beta_samples)[3]
     
     # Symmetric global y-limits across layers
     q_lim <- range(beta_samples, na.rm = TRUE)
     y_max <- max(abs(q_lim))
     y_lim <- c(-y_max, y_max)
     
     # CI probs
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     
     plot_list <- vector("list", K)
     
     for (k in seq_len(K)) {
          # Matrix [iter, p] for layer k
          beta_k <- beta_samples[, , k, drop = FALSE][, , 1]
          
          means    <- colMeans(beta_k, na.rm = TRUE)
          ci_lower <- apply(beta_k, 2, quantile, probs = probs[1], na.rm = TRUE)
          ci_upper <- apply(beta_k, 2, quantile, probs = probs[2], na.rm = TRUE)
          
          df <- data.frame(
               Covariate = factor(covariate_labels, levels = covariate_labels), # keep original order
               Mean  = means,
               Lower = ci_lower,
               Upper = ci_upper,
               stringsAsFactors = FALSE
          )
          
          p_k <- ggplot2::ggplot(df, ggplot2::aes(x = Covariate, y = Mean, color = ifelse(Upper < 0, "red",
                                                                                          ifelse(Lower > 0, "darkgreen", "black")))) +
               ggplot2::geom_segment(ggplot2::aes(xend = Covariate, y = Lower, yend = Upper), size = 0.6) +
               ggplot2::geom_point(size = 1.2) +
               ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
               ggplot2::scale_color_identity() +
               ggplot2::coord_cartesian(ylim = y_lim) +
               ggplot2::labs(title = layer_titles[k], x = "Covariates", y = expression(beta)) +
               ggplot2::theme_minimal(base_size = 11) +
               ggplot2::theme(
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    axis.text.x  = ggplot2::element_text(angle = 0, hjust = 0.5),  # show labels
                    axis.ticks.x = ggplot2::element_line()
               )
          
          plot_list[[k]] <- p_k
     }
     
     # Arrange in a grid (2 columns by default like the example)
     gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
     invisible(plot_list)
}

plot_lambda_by_layer <- function(lambda_samples,
                                 layer_titles = c("Loudness", "Brightness", "Tonality", "Rhythm"),
                                 ci = 0.95) {
     K <- ncol(lambda_samples)
     
     # Layer titles
     if (is.null(layer_titles)) layer_titles <- paste0("Layer ", seq_len(K))
     stopifnot(length(layer_titles) == K)
     
     # Symmetric global y-limits
     q_lim <- range(lambda_samples, na.rm = TRUE)
     y_max <- max(abs(q_lim))
     y_lim <- c(-y_max, y_max)
     
     # CI probs
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     
     # Summaries per layer
     means    <- apply(lambda_samples, 2, mean, na.rm = TRUE)
     ci_lower <- apply(lambda_samples, 2, quantile, probs = probs[1], na.rm = TRUE)
     ci_upper <- apply(lambda_samples, 2, quantile, probs = probs[2], na.rm = TRUE)
     
     df <- data.frame(
          Layer = factor(layer_titles, levels = layer_titles),
          Mean  = means,
          Lower = ci_lower,
          Upper = ci_upper,
          stringsAsFactors = FALSE
     )
     
     # Color: green if CI > 0, red if CI < 0, black otherwise
     df$Color <- ifelse(df$Upper < 0, "red",
                        ifelse(df$Lower > 0, "darkgreen", "black"))
     
     p <- ggplot2::ggplot(df, ggplot2::aes(x = Layer, y = Mean, color = Color)) +
          ggplot2::geom_segment(ggplot2::aes(xend = Layer, y = Lower, yend = Upper), size = 0.7) +
          ggplot2::geom_point(size = 1.8) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
          ggplot2::scale_color_identity() +
          ggplot2::coord_cartesian(ylim = y_lim) +
          ggplot2::labs(title = "Layer-specific bilinear scales",
                        x = "Layer",
                        y = expression(lambda)) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               axis.text.x  = ggplot2::element_text(angle = 0, hjust = 0.5),
               axis.ticks.x = ggplot2::element_line()
          )
     
     print(p)
     invisible(list(plot = p, summary = df))
}

plot_u_norms <- function(U_samples,
                         ci = 0.95) {
     # Expect U_samples as [iterations x n x d]
     S <- dim(U_samples)[1]
     n <- dim(U_samples)[2]
     d <- dim(U_samples)[3]
     
     # R[s, i] = ||u_i|| at iteration s
     R <- sqrt(apply(U_samples^2, c(1, 2), sum))
     
     # Posterior summaries
     means <- colMeans(R, na.rm = TRUE)
     alpha <- 1 - ci
     ci_mat <- t(apply(R, 2, stats::quantile,
                       probs = c(alpha/2, 1 - alpha/2),
                       na.rm = TRUE))
     
     # Ordered data frame (by mean norm)
     df <- data.frame(
          Index = seq_len(n),
          Mean  = means,
          Lower = ci_mat[, 1],
          Upper = ci_mat[, 2],
          stringsAsFactors = FALSE
     )
     df <- df[order(df$Mean), , drop = FALSE]
     df$Rank <- seq_len(nrow(df))
     
     # Threshold: r_ci = sqrt(qchisq(ci, df = d))
     eps <- sqrt(stats::qchisq(p = 0.5, df = d))
     
     # Color: CI strictly above eps => darkgreen; else black
     df$Color <- ifelse(df$Lower > eps, "darkgreen", "black")
     
     # y-limits (nonnegative)
     y_lim <- c(0, max(df$Upper, na.rm = TRUE))
     
     # Plot (no actor labels on x-axis)
     p <- ggplot2::ggplot(df) +
          ggplot2::geom_segment(
               ggplot2::aes(x = Rank, xend = Rank, y = Lower, yend = Upper, color = Color),
               linewidth = 0.6
          ) +
          ggplot2::geom_point(
               ggplot2::aes(x = Rank, y = Mean, color = Color),
               size = 0.6
          ) +
          ggplot2::geom_hline(yintercept = eps, linetype = "dashed", color = "blue") +
          ggplot2::scale_x_continuous(breaks = NULL, labels = NULL) +  # no labels
          ggplot2::labs(
               title = "",
               x = "Actors (ordered by mean norm)",
               y = expression(plain("||") * u[i] * plain("||"))
          ) +
          ggplot2::coord_cartesian(ylim = y_lim) +
          ggplot2::scale_color_identity() +
          ggplot2::theme_minimal() +
          ggplot2::theme(
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank()
          )
     
     print(p)
     invisible(list(plot = p, summary = df, threshold = eps))
}

plot_geometry_strength_by_layer <- function(U_samples,
                                            lambda_samples,
                                            layer_titles = c("Loudness","Brightness","Tonality","Rhythm"),
                                            ci = 0.95,
                                            ncol = 2) {
     S <- dim(U_samples)[1]
     n <- dim(U_samples)[2]
     d <- dim(U_samples)[3]
     K <- ncol(lambda_samples)
     
     if (is.null(layer_titles)) layer_titles <- paste0("Layer ", seq_len(K))
     
     R <- sqrt(apply(U_samples^2, c(1, 2), sum))  # S x n
     
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     
     plot_list <- vector("list", K)
     summaries <- vector("list", K)
     
     # Global y-limits across layers
     all_vals <- c()
     for (k in seq_len(K)) {
          lam_k    <- abs(lambda_samples[, k])
          Gk       <- R * matrix(lam_k, nrow = S, ncol = n)
          all_vals <- c(all_vals, as.vector(Gk))
     }
     y_lim <- c(0, max(all_vals, na.rm = TRUE))
     
     # Per layer summaries & plots
     for (k in seq_len(K)) {
          lam_k <- abs(lambda_samples[, k])
          Gk <- R * matrix(lam_k, nrow = S, ncol = n)  # S x n
          
          means <- colMeans(Gk, na.rm = TRUE)
          ci_mat <- t(apply(Gk, 2, stats::quantile, probs = probs, na.rm = TRUE))
          
          df <- data.frame(
               Index = seq_len(n),
               Mean  = means,
               Lower = ci_mat[, 1],
               Upper = ci_mat[, 2],
               stringsAsFactors = FALSE
          )
          df <- df[order(df$Mean), , drop = FALSE]
          df$Rank <- seq_len(nrow(df))
          
          eps <- mean(abs(lambda_samples[, k])) * sqrt(qchisq(p = 0.50, df = d))
          
          # Colors (nonnegative)
          df$Color <- ifelse(df$Lower > eps, "darkgreen", "black")
          
          p_k <- ggplot2::ggplot(df) +
               ggplot2::geom_segment(
                    ggplot2::aes(x = Rank, xend = Rank, y = Lower, yend = Upper, color = Color),
                    linewidth = 0.3
               ) +
               ggplot2::geom_point(
                    ggplot2::aes(x = Rank, y = Mean, color = Color),
                    size = 0.3
               ) +
               ggplot2::geom_hline(yintercept = eps, linetype = "dashed", color = "blue") +
               ggplot2::scale_x_continuous(breaks = NULL, labels = NULL) +
               ggplot2::labs(
                    title = layer_titles[k],
                    x = "Actors (ordered by mean strength)",
                    y = expression(abs(lambda[k]) %*% plain("||") * u[i] * plain("||"))
               ) +
               ggplot2::coord_cartesian(ylim = y_lim) +
               ggplot2::scale_color_identity() +
               ggplot2::theme_minimal(base_size = 11) +
               ggplot2::theme(
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank()
               )
          
          plot_list[[k]] <- p_k
          summaries[[k]] <- df
     }
     
     # arrange
     gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
     
     invisible(list(plots = plot_list, summaries = summaries, y_lim = y_lim, ci = ci))
}

plot_latent_positions <- function(samples,
                                  mean_point_size = 3,
                                  label_points = FALSE,
                                  label_size = 3,
                                  alpha = 0.6,
                                  seed = 42) {
     set.seed(seed)
     
     suppressMessages(suppressWarnings(library(ggplot2)))
     suppressMessages(suppressWarnings(library(MCMCpack)))
     
     # Inputs & checks
     Uarr <- samples$U                    # dims: S x n x d
     stopifnot(length(dim(Uarr)) == 3)
     S <- dim(Uarr)[1]; n <- dim(Uarr)[2]; d <- dim(Uarr)[3]
     if (d < 2) stop("Need at least 2 latent dimensions to plot.")
     
     # Reference configuration: first draw, centered & scaled
     U0 <- scale(Uarr[1, , ], center = TRUE, scale = TRUE)  # n x d
     
     # Procrustes-align all draws to U0
     U_aligned <- array(NA_real_, c(S, n, d))
     for (b in seq_len(S)) {
          Ub <- scale(Uarr[b, , ], center = TRUE, scale = TRUE)
          U_aligned[b, , ] <- MCMCpack::procrustes(
               X = Ub, Xstar = U0, translation = TRUE, dilation = TRUE
          )$X.new
     }
     
     # Posterior mean positions (n x d)
     U_pm <- apply(U_aligned, c(2, 3), mean)
     
     # Choose two plotting dimensions (highest variance of U_pm)
     dims_idx <- if (d > 2) order(apply(U_pm, 2, stats::var), decreasing = TRUE)[1:2] else 1:2
     
     # Colors from reference U0 in the same two dims
     U0_plot <- U0[, dims_idx, drop = FALSE]
     rr <- atan2(U0_plot[, 2], U0_plot[, 1])
     rr <- (rr - min(rr)) / (max(rr) - min(rr) + 1e-12)
     gg <- 1 - rr
     bb <- rowSums(U0_plot^2)
     bb <- (bb - min(bb)) / (max(bb) + 1e-12)
     col_vec <- rgb(rr, gg, bb, alpha = alpha)
     
     # Plot data (mean positions)
     df <- data.frame(
          Dim1 = U_pm[, dims_idx[1]],
          Dim2 = U_pm[, dims_idx[2]],
          Color = col_vec,
          Vertex = seq_len(n)
     )
     
     p <- ggplot(df, aes(x = Dim1, y = Dim2)) +
          geom_point(aes(color = Color), size = mean_point_size) +
          { if (label_points) geom_text(aes(label = Vertex), color = "black", size = label_size, vjust = 1.4) else NULL } +
          scale_color_identity() +
          labs(
               x = "Dim 1", y = "Dim 2",
               title = "Posterior mean latent positions"
          ) +
          theme_minimal() +
          theme(panel.grid.minor = element_blank())
     
     print(p)
     invisible(list(plot = p, U_mean = U_pm, U0 = U0, dims_used = dims_idx, colors = col_vec))
}

ap_score <- function(label, score) {
     o <- order(score, decreasing = TRUE)
     y_ <- label[o]
     cum_tp <- cumsum(y_ == 1L)
     prec   <- cum_tp / seq_along(y_)
     mean(prec[y_ == 1L])
}

prob_metrics_per_layer_model_1 <- function(samples,
                                           y,
                                           thin = 50,
                                           ci = 0.95) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     
     # Helpers
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     qmat  <- function(M) t(apply(M, 2, stats::quantile, probs = probs, na.rm = TRUE))
     cmean <- function(M) colMeans(M, na.rm = TRUE)
     clamp <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
     
     # Pull observed upper-triangle labels once per layer
     y_up_list <- vector("list", K)
     for (k in seq_len(K)) {
          y_up_list[[k]] <- y[, , k][up]
     }
     
     # Storage
     auc     <- matrix(NA_real_, S_used, K)
     ap      <- matrix(NA_real_, S_used, K)
     brier   <- matrix(NA_real_, S_used, K)
     logloss <- matrix(NA_real_, S_used, K)
     
     # Loop over posterior draws
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          # Samples for this draw
          mu_s    <- samples$mu   [s, ] 
          delta_s <- samples$delta[s, , ]
          zeta_s  <- samples$zeta[s]
          
          # Predicted probabilities (now include zeta)
          P_s <- interaction_prob_cpp(mu_s, delta_s, zeta_s)
          
          # Per-layer metrics on upper triangle
          for (k in seq_len(K)) {
               pv <- clamp(P_s[, , k][up])
               yv <- y_up_list[[k]]
               
               roc_k <- pROC::roc(response = yv, predictor = pv,
                                  levels = c(0, 1), quiet = TRUE, direction = "<")
               
               auc    [ii, k] <- as.numeric(pROC::auc(roc_k))
               ap     [ii, k] <- ap_score(yv, pv)
               brier  [ii, k] <- mean((yv - pv)^2)
               logloss[ii, k] <- -mean(yv * log(pv) + (1 - yv) * log(1 - pv))
          }
     }
     
     list(
          auc     = list(samples = auc,     mean = cmean(auc),     ci = qmat(auc)),
          auprc   = list(samples = ap,      mean = cmean(ap),      ci = qmat(ap)),
          brier   = list(samples = brier,   mean = cmean(brier),   ci = qmat(brier)),
          logloss = list(samples = logloss, mean = cmean(logloss), ci = qmat(logloss))
     )
}

prob_metrics_per_layer_model_2 <- function(samples, 
                                           X, 
                                           y,
                                           thin = 50, 
                                           ci = 0.95) {
     
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Helpers
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     qmat  <- function(M) t(apply(M, 2, stats::quantile, probs = probs, na.rm = TRUE))
     cmean <- function(M) colMeans(M, na.rm = TRUE)
     clamp <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
     
     # Pull observed labels once per layer
     y_up_list <- vector("list", K)
     for (k in seq_len(K)) {
          y_up_list[[k]] <- y[ , , k][up]
     }
     
     # Storage
     auc     <- matrix(NA_real_, S_used, K)
     ap      <- matrix(NA_real_, S_used, K)
     brier   <- matrix(NA_real_, S_used, K)
     logloss <- matrix(NA_real_, S_used, K)
     
     # Loop over MCMC draws
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          # Samples
          mu_s    <- samples$mu   [s, ]
          delta_s <- samples$delta[s, , ]
          beta_s  <- samples$beta [s, , ]
          zeta_s  <- samples$zeta [s]
          
          # Predicted probabilities (now includes zeta)
          P_s <- interaction_prob_cpp(mu_s, delta_s, X, beta_s, zeta_s)
          
          # per-layer metrics
          for (k in seq_len(K)) {
               pv <- clamp(P_s[ , , k][up])
               yv <- y_up_list[[k]]
               
               roc_k <- pROC::roc(response = yv, predictor = pv, levels = c(0, 1), quiet = TRUE, direction = "<")
               
               auc    [ii, k] <- as.numeric(pROC::auc(roc_k))
               ap     [ii, k] <- ap_score(yv, pv)
               brier  [ii, k] <- mean((yv - pv)^2)
               logloss[ii, k] <- -mean(yv * log(pv) + (1 - yv) * log(1 - pv))
          }
     }
     
     list(
          auc     = list(samples = auc,     mean = cmean(auc),     ci = qmat(auc)),
          auprc   = list(samples = ap,      mean = cmean(ap),      ci = qmat(ap)),
          brier   = list(samples = brier,   mean = cmean(brier),   ci = qmat(brier)),
          logloss = list(samples = logloss, mean = cmean(logloss), ci = qmat(logloss))
     )
}

prob_metrics_per_layer_model_3 <- function(samples, 
                                           X, 
                                           y,
                                           thin = 50, 
                                           ci = 0.95) {
     
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$U)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Helpers
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     qmat  <- function(M) t(apply(M, 2, stats::quantile, probs = probs, na.rm = TRUE))
     cmean <- function(M) colMeans(M, na.rm = TRUE)
     clamp <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
     
     # Observed labels per layer
     y_up_list <- vector("list", K)
     for (k in seq_len(K)) {
          y_up_list[[k]] <- y[ , , k][up]
     }
     
     # Storage
     auc     <- matrix(NA_real_, S_used, K)
     ap      <- matrix(NA_real_, S_used, K)
     brier   <- matrix(NA_real_, S_used, K)
     logloss <- matrix(NA_real_, S_used, K)
     
     # Loop over MCMC draws
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          # Samples
          mu_s     <- samples$mu    [s, ]
          delta_s  <- samples$delta [s, , ]
          beta_s   <- samples$beta  [s, , ]
          U_s      <- samples$U     [s, , ]
          lambda_s <- samples$lambda[s, ]
          zeta_s   <- samples$zeta  [s]
          
          # Predicted probabilities
          P_s <- interaction_prob_cpp(mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s)
          
          # per-layer metrics
          for (k in seq_len(K)) {
               pv <- clamp(P_s[ , , k][up])
               yv <- y_up_list[[k]]
               
               roc_k <- pROC::roc(response = yv, predictor = pv, levels = c(0, 1),
                                  quiet = TRUE, direction = "<")
               
               auc    [ii, k] <- as.numeric(pROC::auc(roc_k))
               ap     [ii, k] <- ap_score(yv, pv)
               brier  [ii, k] <- mean((yv - pv)^2)
               logloss[ii, k] <- -mean(yv * log(pv) + (1 - yv) * log(1 - pv))
          }
     }
     
     list(
          auc     = list(samples = auc,     mean = cmean(auc),     ci = qmat(auc)),
          auprc   = list(samples = ap,      mean = cmean(ap),      ci = qmat(ap)),
          brier   = list(samples = brier,   mean = cmean(brier),   ci = qmat(brier)),
          logloss = list(samples = logloss, mean = cmean(logloss), ci = qmat(logloss))
     )
}

prob_metrics_per_layer_model_4 <- function(samples, 
                                           X, 
                                           y,
                                           thin = 50, 
                                           ci = 0.95) {
     
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$U)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Helpers
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     qmat  <- function(M) t(apply(M, 2, stats::quantile, probs = probs, na.rm = TRUE))
     cmean <- function(M) colMeans(M, na.rm = TRUE)
     clamp <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
     
     # Observed labels per layer
     y_up_list <- vector("list", K)
     for (k in seq_len(K)) {
          y_up_list[[k]] <- y[ , , k][up]
     }
     
     # Storage
     auc     <- matrix(NA_real_, S_used, K)
     ap      <- matrix(NA_real_, S_used, K)
     brier   <- matrix(NA_real_, S_used, K)
     logloss <- matrix(NA_real_, S_used, K)
     
     # Loop over MCMC draws
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          # Samples
          mu_s     <- samples$mu    [s, ]
          delta_s  <- samples$delta [s, , ]
          beta_s   <- samples$beta  [s, , ]
          U_s      <- samples$U     [s, , ]
          lambda_s <- samples$lambda[s, ]
          zeta_s   <- samples$zeta  [s]
          
          # Predicted probabilities
          P_s <- interaction_prob_cpp(mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s)
          
          # per-layer metrics
          for (k in seq_len(K)) {
               pv <- clamp(P_s[ , , k][up])
               yv <- y_up_list[[k]]
               
               roc_k <- pROC::roc(response = yv, predictor = pv, levels = c(0, 1),
                                  quiet = TRUE, direction = "<")
               
               auc    [ii, k] <- as.numeric(pROC::auc(roc_k))
               ap     [ii, k] <- ap_score(yv, pv)
               brier  [ii, k] <- mean((yv - pv)^2)
               logloss[ii, k] <- -mean(yv * log(pv) + (1 - yv) * log(1 - pv))
          }
     }
     
     list(
          auc     = list(samples = auc,     mean = cmean(auc),     ci = qmat(auc)),
          auprc   = list(samples = ap,      mean = cmean(ap),      ci = qmat(ap)),
          brier   = list(samples = brier,   mean = cmean(brier),   ci = qmat(brier)),
          logloss = list(samples = logloss, mean = cmean(logloss), ci = qmat(logloss))
     )
}

prob_metrics_per_layer_model_5 <- function(samples, 
                                           X, 
                                           y,
                                           thin = 50, 
                                           ci = 0.95) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     M  <- sum(up)
     
     # Helpers
     alpha <- (1 - ci) / 2
     probs <- c(alpha, 1 - alpha)
     qmat  <- function(M) t(apply(M, 2, stats::quantile, probs = probs, na.rm = TRUE))
     cmean <- function(M) colMeans(M, na.rm = TRUE)
     clamp <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
     
     # Observed labels per layer
     y_up_list <- vector("list", K)
     for (k in seq_len(K)) y_up_list[[k]] <- y[ , , k][up]
     
     # Storage
     auc     <- matrix(NA_real_, S_used, K)
     ap      <- matrix(NA_real_, S_used, K)
     brier   <- matrix(NA_real_, S_used, K)
     logloss <- matrix(NA_real_, S_used, K)
     
     # Loop over MCMC draws
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          # Samples
          mu_s    <- samples$mu   [s, ] 
          delta_s <- samples$delta[s, , ]
          beta_s  <- samples$beta [s, , ]
          zeta_s  <- samples$zeta [s] 
          Theta_s <- samples$Theta[ , , ((s - 1L) * K + 1L):(s * K)]
          xi_s    <- samples$xi   [s, , ]
          storage.mode(xi_s) <- "integer"
          
          # Predicted probabilities
          P_s <- interaction_prob_cpp(mu_s, delta_s, X, beta_s, Theta_s, xi_s, zeta_s)
          
          # per-layer metrics
          for (k in seq_len(K)) {
               pv <- clamp(P_s[ , , k][up])
               yv <- y_up_list[[k]]
               
               roc_k <- pROC::roc(response = yv, predictor = pv, levels = c(0, 1),
                                  quiet = TRUE, direction = "<")
               
               auc    [ii, k] <- as.numeric(pROC::auc(roc_k))
               ap     [ii, k] <- ap_score(yv, pv)
               brier  [ii, k] <- mean((yv - pv)^2)
               logloss[ii, k] <- -mean(yv * log(pv) + (1 - yv) * log(1 - pv))
          }
     }
     
     list(
          auc     = list(samples = auc,     mean = cmean(auc),     ci = qmat(auc)),
          auprc   = list(samples = ap,      mean = cmean(ap),      ci = qmat(ap)),
          brier   = list(samples = brier,   mean = cmean(brier),   ci = qmat(brier)),
          logloss = list(samples = logloss, mean = cmean(logloss), ci = qmat(logloss))
     )
}

log_mean_exp <- function(v) {
     m <- max(v)
     m + log(mean(exp(v - m)))
}

dic_model_1 <- function(log_lik, 
                        samples, 
                        y, 
                        thin = 50,
                        eps = 1e-12) {
     # Total draws
     S_tot <- nrow(samples$mu)
     
     # Thinning index
     s_index <- seq(from = thin, to = S_tot, by = thin)
     
     # Apply thinning to likelihood and samples
     log_lik_th <- log_lik      [s_index]
     mu_s_th    <- samples$mu   [s_index, , drop = FALSE]
     delta_s_th <- samples$delta[s_index, , , drop = FALSE]
     zeta_s_th  <- samples$zeta [s_index]
     
     # Dbar from per-iteration log-likelihoods
     D_s  <- -2 * log_lik_th
     Dbar <- mean(D_s)
     
     # Dhat: deviance at posterior means (theta-hat)
     mu_hat    <- colMeans(mu_s_th)
     delta_hat <- apply(delta_s_th, c(2, 3), mean)
     zeta_hat  <- mean(zeta_s_th)
     
     # Evaluate log-likelihood at (mu_hat, delta_hat, zeta_hat)
     log_lik_hat <- log_likelihood_iter_cpp(y, mu_hat, delta_hat, zeta_hat, eps)
     
     Dhat <- -2 * log_lik_hat
     pD   <- Dbar - Dhat
     DIC  <- Dbar + pD  # equivalently: 2*Dbar - Dhat
     
     list(
          DIC  = DIC,
          Dbar = Dbar,
          Dhat = Dhat,
          pD   = pD
     )
}

dic_model_2 <- function(log_lik, 
                        samples, 
                        y, 
                        X,
                        thin = 50,
                        eps = 1e-12) {
     # Total draws
     S_tot <- nrow(samples$mu)
     
     # Thinning index
     s_index <- seq(from = thin, to = S_tot, by = thin)
     
     # Apply thinning to likelihood and samples
     log_lik_th <- log_lik                  [s_index]
     mu_s_th    <- samples$mu               [s_index, , drop = FALSE]
     delta_s_th <- samples$delta            [s_index, , , drop = FALSE]
     beta_s_th  <- samples$beta             [s_index, , , drop = FALSE]
     zeta_s_th  <- as.numeric(samples$zeta) [s_index]
     
     # Dbar from per-iteration log-likelihoods
     D_s  <- -2 * log_lik_th
     Dbar <- mean(D_s)
     
     # Posterior means (theta-hat) using thinned draws
     mu_hat    <- colMeans(mu_s_th)                           # length K
     delta_hat <- apply(delta_s_th, c(2, 3), mean)            # n x K
     beta_hat  <- apply(beta_s_th,  c(2, 3), mean)            # p x K
     zeta_hat  <- mean(zeta_s_th)                             
     
     # Deviance at posterior means
     log_lik_hat <- log_likelihood_iter_cpp(y, mu_hat, delta_hat, X, beta_hat, zeta_hat, eps)
     Dhat <- -2 * log_lik_hat
     
     pD  <- Dbar - Dhat
     DIC <- Dbar + pD  # == 2*Dbar - Dhat
     
     list(
          DIC  = DIC,
          Dbar = Dbar,
          Dhat = Dhat,
          pD   = pD
     )
}

dic_model_3 <- function(log_lik,
                        samples,
                        y,
                        X,
                        thin = 1,
                        eps = 1e-12) {
     # Indices to use (thinning)
     S_tot   <- nrow(samples$mu)
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Dbar from per-iteration log-likelihoods (thinned)
     D_s  <- -2 * log_lik[s_index]
     Dbar <- mean(D_s)
     
     # Posterior means (theta-hat) using thinned draws
     mu_hat     <- colMeans(samples$mu[s_index, , drop = FALSE])                         # [K]
     delta_hat  <- apply(samples$delta[s_index, , , drop = FALSE], c(2, 3), mean)        # [n x K]
     beta_hat   <- apply(samples$beta [s_index, , , drop = FALSE], c(2, 3), mean)        # [p x K]
     lambda_hat <- colMeans(samples$lambda[s_index, , drop = FALSE])                     # [K]
     zeta_hat   <- mean(samples$zeta[s_index])                                           # scalar
     
     # S_hat = E[U U^T | data] (rotation-invariant), averaged over thinned draws
     n <- dim(samples$U)[2]
     d <- dim(samples$U)[3]
     S_hat <- matrix(0, n, n)
     for (s in s_index) {
          U_s   <- samples$U[s, , ]              # [n x d]
          S_hat <- S_hat + U_s %*% t(U_s)
     }
     S_hat <- S_hat / S_used
     
     # Deviance at posterior means
     log_lik_hat <- log_likelihood_iter_cpp(
          y, mu_hat, delta_hat, X, beta_hat, S_hat, lambda_hat, zeta_hat, eps
     )
     Dhat <- -2 * log_lik_hat
     
     # DIC components
     pD  <- Dbar - Dhat
     DIC <- Dbar + pD  # == 2*Dbar - Dhat
     
     list(
          DIC  = DIC,
          Dbar = Dbar,
          Dhat = Dhat,
          pD   = pD
     )
}

dic_model_4 <- function(log_lik,
                        samples,
                        y,
                        X,
                        thin = 1,
                        eps = 1e-12) {
     # Indices to use (thinning)
     S_tot   <- nrow(samples$mu)
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Dbar from per-iteration log-likelihoods (thinned)
     D_s  <- -2 * log_lik[s_index]
     Dbar <- mean(D_s)
     
     # Posterior means (theta-hat) using thinned draws
     mu_hat     <- colMeans(samples$mu[s_index, , drop = FALSE])                  # [K]
     delta_hat  <- apply(samples$delta[s_index, , , drop = FALSE], c(2, 3), mean) # [n x K]
     beta_hat   <- apply(samples$beta [s_index, , , drop = FALSE], c(2, 3), mean) # [p x K]
     lambda_hat <- colMeans(samples$lambda[s_index, , drop = FALSE])              # [K]
     zeta_hat   <- mean(samples$zeta[s_index])                                    # scalar
     
     # S_hat = E[U U^T | data] (rotation-invariant), averaged over thinned draws
     n <- dim(samples$U)[2]
     d <- dim(samples$U)[3]
     S_hat <- matrix(0, n, n)
     for (s in s_index) {
          U_s   <- samples$U[s, , ]             # [n x d]
          S_hat <- S_hat + U_s %*% t(U_s)
     }
     S_hat <- S_hat / S_used
     
     # Deviance at posterior means
     log_lik_hat <- log_likelihood_iter_cpp(
          y, mu_hat, delta_hat, X, beta_hat, S_hat, lambda_hat, zeta_hat, eps
     )
     Dhat <- -2 * log_lik_hat
     
     # DIC components
     pD  <- Dbar - Dhat
     DIC <- Dbar + pD  # == 2*Dbar - Dhat
     
     list(
          DIC  = DIC,
          Dbar = Dbar,
          Dhat = Dhat,
          pD   = pD
     )
}

dic_model_5 <- function(log_lik,
                        samples,
                        y,
                        X,
                        thin = 50,
                        eps = 1e-12) {
     # Thinning indices
     S_tot   <- nrow(samples$mu)
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Dbar from thinned draws
     D_s  <- -2 * log_lik[s_index]
     Dbar <- mean(D_s)
     
     # Build posterior-mean probabilities p_hat(i,j,k) = average_s p_ijk^(s)
     n <- dim(y)[1]
     K <- dim(y)[3]
     up <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
     
     P_sum <- array(0, dim = c(n, n, K))
     
     for (s in s_index) {
          mu_s    <- samples$mu   [s, ]
          delta_s <- samples$delta[s, , ]
          beta_s  <- samples$beta [s, , ]
          zeta_s  <- samples$zeta [s]
          # Theta is stored as C x C x (K * S); take this draw's K slices:
          Theta_s <- samples$Theta[ , , ((s - 1L) * K + 1L):(s * K), drop = FALSE]
          Theta_s <- array(Theta_s, dim = c(dim(samples$Theta)[1], dim(samples$Theta)[2], K))
          # xi for this draw: n x K (0-based ints are fine; arma::umat in C++)
          xi_s    <- samples$xi[s, , , drop = FALSE]
          xi_s    <- matrix(xi_s, nrow = n, ncol = K)
          
          # Probabilities for this draw
          P_s <- interaction_prob_cpp(mu_s, delta_s, X, beta_s, Theta_s, xi_s, zeta_s)
          P_sum <- P_sum + P_s
     }
     
     P_hat <- P_sum / S_used
     
     # Deviance at posterior-mean probabilities
     clamp <- function(p) pmin(pmax(p, eps), 1 - eps)
     Dhat <- 0
     for (k in seq_len(K)) {
          pv <- clamp(P_hat[ , , k][up])
          yv <- y     [ , , k][up]
          Dhat <- Dhat - 2 * sum(yv * log(pv) + (1 - yv) * log(1 - pv))
     }
     
     # DIC components (nonnegative pD by construction)
     pD  <- Dbar - Dhat
     DIC <- Dbar + pD
     
     list(DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD)
}

waic_model_1 <- function(samples, 
                         y, 
                         thin = 50, 
                         eps = 1e-12) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Build pointwise log-likelihood matrix: rows = draws, cols = data points
     M <- n * (n - 1) / 2
     N <- K * M
     ll_mat <- matrix(NA_real_, nrow = S_used, ncol = N)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]
          delta_s <- samples$delta[s, , ]
          zeta_s  <- as.numeric(samples$zeta[s]) 
          
          ll_mat[ii, ] <- log_likelihood_pointwise_cpp(y, mu_s, delta_s, zeta_s, eps)
     }
     
     # WAIC components
     lppd_i  <- apply(ll_mat, 2L, log_mean_exp) 
     pwaic_i <- apply(ll_mat, 2L, stats::var) 
     lppd    <- sum(lppd_i)
     p_waic  <- sum(pwaic_i)
     WAIC    <- -2 * (lppd - p_waic)
     
     # Rough SE for WAIC (Vehtari et al.)
     contrib <- -2 * (lppd_i - pwaic_i)
     waic_se <- sqrt(length(contrib) * var(contrib))
     
     list(
          WAIC   = WAIC,
          lppd   = lppd,
          p_waic = p_waic,
          se     = waic_se
     )
}

waic_model_2 <- function(samples, 
                         y, 
                         X,
                         thin = 50, 
                         eps = 1e-12) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Pointwise log-likelihood matrix: rows = draws, cols = data points
     M <- n * (n - 1) / 2
     N <- K * M
     ll_mat <- matrix(NA_real_, nrow = S_used, ncol = N)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]         # length K
          delta_s <- samples$delta[s, , ]       # n x K
          beta_s  <- samples$beta [s, , ]       # p x K
          zeta_s  <- samples$zeta [s]
          
          ll_mat[ii, ] <- log_likelihood_pointwise_cpp(y, mu_s, delta_s, X, beta_s, zeta_s, eps)
     }
     
     # WAIC components
     lppd_i  <- apply(ll_mat, 2L, log_mean_exp)  # length N
     pwaic_i <- apply(ll_mat, 2L, stats::var)    # length N
     lppd    <- sum(lppd_i)
     p_waic  <- sum(pwaic_i)
     WAIC    <- -2 * (lppd - p_waic)
     
     # Rough SE for WAIC (Vehtari et al.)
     contrib <- -2 * (lppd_i - pwaic_i)
     waic_se <- sqrt(length(contrib) * var(contrib))
     
     list(
          WAIC   = WAIC,
          lppd   = lppd,
          p_waic = p_waic,
          se     = waic_se
     )
}

waic_model_3 <- function(samples,
                         y,
                         X,
                         thin = 50,
                         eps = 1e-12) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Pointwise log-likelihood matrix: rows = draws, cols = data points
     M <- n * (n - 1) / 2
     N <- K * M
     ll_mat <- matrix(NA_real_, nrow = S_used, ncol = N)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s     <- samples$mu    [s, ]    # length K
          delta_s  <- samples$delta [s, , ]  # n x K
          beta_s   <- samples$beta  [s, , ]  # p x K
          U_s      <- samples$U     [s, , ]  # n x d
          lambda_s <- samples$lambda[s, ]    # length K
          zeta_s   <- samples$zeta  [s]
          
          ll_mat[ii, ] <- log_likelihood_pointwise_cpp(
               y, mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s, eps
          )
     }
     
     # WAIC components
     lppd_i  <- apply(ll_mat, 2L, log_mean_exp)  # length N
     pwaic_i <- apply(ll_mat, 2L, stats::var)    # length N
     lppd    <- sum(lppd_i)
     p_waic  <- sum(pwaic_i)
     WAIC    <- -2 * (lppd - p_waic)
     
     # Rough SE for WAIC (Vehtari et al.)
     contrib <- -2 * (lppd_i - pwaic_i)
     waic_se <- sqrt(length(contrib) * var(contrib))
     
     list(
          WAIC   = WAIC,
          lppd   = lppd,
          p_waic = p_waic,
          se     = waic_se
     )
}

waic_model_4 <- function(samples,
                         y,
                         X,
                         thin = 50,
                         eps = 1e-12) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Pointwise log-likelihood matrix: rows = draws, cols = data points
     M <- n * (n - 1) / 2
     N <- K * M
     ll_mat <- matrix(NA_real_, nrow = S_used, ncol = N)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s     <- samples$mu    [s, ]   # length K
          delta_s  <- samples$delta [s, , ] # n x K
          beta_s   <- samples$beta  [s, , ] # p x K
          U_s      <- samples$U     [s, , ] # n x d
          lambda_s <- samples$lambda[s, ]   # length K
          zeta_s   <- samples$zeta  [s]     # scalar
          
          # Model 4 pointwise log-likelihood (-exp(lambda)*||u_i-u_j|| + zeta)
          ll_mat[ii, ] <- log_likelihood_pointwise_cpp(
               y, mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s, eps
          )
     }
     
     # WAIC components
     lppd_i  <- apply(ll_mat, 2L, log_mean_exp)  # length N
     pwaic_i <- apply(ll_mat, 2L, stats::var)    # length N
     lppd    <- sum(lppd_i)
     p_waic  <- sum(pwaic_i)
     WAIC    <- -2 * (lppd - p_waic)
     
     # Rough SE for WAIC (Vehtari et al.)
     contrib <- -2 * (lppd_i - pwaic_i)
     waic_se <- sqrt(length(contrib) * var(contrib))
     
     list(
          WAIC   = WAIC,
          lppd   = lppd,
          p_waic = p_waic,
          se     = waic_se
     )
}

waic_model_5 <- function(samples,
                         y,
                         X,
                         thin = 50,
                         eps = 1e-12) {
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Pointwise log-likelihood matrix: rows = draws, cols = data points
     M <- n * (n - 1) / 2
     N <- K * M
     ll_mat <- matrix(NA_real_, nrow = S_used, ncol = N)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]       # length K
          delta_s <- samples$delta[s, , ]     # n x K
          beta_s  <- samples$beta [s, , ]     # p x K
          zeta_s  <- samples$zeta [s]         # scalar
          
          # Theta: C x C x (K * S) -> take this sample's K slices (C x C x K)
          Theta_s <- samples$Theta[ , , ((s - 1L) * K + 1L):(s * K)]
          
          # xi: expected S x n; take row s (length n, zero-based labels)
          xi_s <- samples$xi[s, , ]
          
          # Model 5 pointwise log-likelihood
          ll_mat[ii, ] <- log_likelihood_pointwise_cpp(
               y, mu_s, delta_s, X, beta_s, Theta_s, xi_s, zeta_s, eps
          )
     }
     
     # WAIC components
     lppd_i  <- apply(ll_mat, 2L, log_mean_exp)  # length N
     pwaic_i <- apply(ll_mat, 2L, stats::var)    # length N
     lppd    <- sum(lppd_i)
     p_waic  <- sum(pwaic_i)
     WAIC    <- -2 * (lppd - p_waic)
     
     # Rough SE for WAIC (Vehtari et al.)
     contrib <- -2 * (lppd_i - pwaic_i)
     waic_se <- sqrt(length(contrib) * var(contrib))
     
     list(
          WAIC   = WAIC,
          lppd   = lppd,
          p_waic = p_waic,
          se     = waic_se
     )
}

# Helpers for distances/diameter on disconnected graphs
safe_mean_distance <- function(g) {
     x <- suppressWarnings(try(igraph::mean_distance(g, directed = FALSE, unconnected = TRUE), silent = TRUE))
     if (inherits(x, "try-error") || !is.finite(x)) NA_real_ else as.numeric(x)
}

safe_diameter <- function(g) {
     x <- suppressWarnings(try(igraph::diameter(g, directed = FALSE, unconnected = TRUE), silent = TRUE))
     if (inherits(x, "try-error") || !is.finite(x)) NA_real_ else as.numeric(x)
}

plot_metrics <- function(alpha, layer_titles, summary, observed, y_label) {
     df <- data.frame(
          Layer    = factor(layer_titles, levels = layer_titles),  # keep given order
          Mean     = summary$Mean,
          Lower    = summary$Lower,
          Upper    = summary$Upper,
          Observed = observed
     )
     ggplot(df, aes(x = Layer)) +
          geom_linerange(aes(ymin = Lower, ymax = Upper), color = "black", linewidth = 0.7) +
          geom_point(aes(y = Mean),     color = "black", size = 1.5) +
          geom_point(aes(y = Observed), color = "red",  size = 2.5, shape = 18) +
          labs(title = "", y = y_label, x = "Layers") +
          theme_minimal(base_size = 12) +
          theme(
               panel.grid = element_blank(),
               axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
          )
}

plot_test_stats_model_1 <- function(samples, 
                                    y,
                                    layer_titles = c("Loudness","Brightness","Tonality","Rhythm"),
                                    ci = 0.95,
                                    thin = 50,
                                    seed = 42) {
     set.seed(seed)
     
     # Shapes
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Storage for simulated metrics (S_used x K)
     density_sim        <- matrix(NA_real_, nrow = S_used, ncol = K)
     transitivity_sim   <- matrix(NA_real_, nrow = S_used, ncol = K)
     assortativity_sim  <- matrix(NA_real_, nrow = S_used, ncol = K)
     loc_trans_mean_sim <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_mean_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_sd_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     mean_dist_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     diameter_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     n_comm_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     modularity_sim     <- matrix(NA_real_, nrow = S_used, ncol = K)
     gcc_frac_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     triangles_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     
     # Simulate graphs & compute metrics
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]
          delta_s <- samples$delta[s, , ]
          zeta_s  <- as.numeric(samples$zeta[s])
          
          y_sim <- simulate_multilayer_network_cpp(mu_s, delta_s, zeta_s)
          
          for (k in seq_len(K)) {
               g <- igraph::graph_from_adjacency_matrix(y_sim[ , , k], mode = "undirected", diag = FALSE)
               
               # Degree stats
               deg <- igraph::degree(g)
               deg_mean_sim[ii, k] <- mean(deg)
               deg_sd_sim  [ii, k] <- stats::sd(deg)
               
               # Density / clustering
               density_sim     [ii, k] <- igraph::edge_density(g, loops = FALSE)
               transitivity_sim[ii, k] <- igraph::transitivity(g, type = "global")
               loc_vec <- igraph::transitivity(g, type = "localundirected", isolates = "zero")
               loc_trans_mean_sim[ii, k] <- mean(loc_vec, na.rm = TRUE)
               
               # Assortativity
               assortativity_sim[ii, k] <- igraph::assortativity_degree(g, directed = FALSE)
               
               # Distances / diameter (safe)
               mean_dist_sim[ii, k] <- suppressWarnings(tryCatch(
                    igraph::mean_distance(g, directed = FALSE, unconnected = TRUE),
                    error = function(e) NA_real_))
               diameter_sim [ii, k] <- suppressWarnings(tryCatch(
                    igraph::diameter(g, directed = FALSE, unconnected = TRUE),
                    error = function(e) NA_real_))
               
               # Components
               comps <- igraph::components(g)
               gcc_frac_sim[ii, k] <- max(comps$csize) / igraph::vcount(g)
               
               # No. communities and modularity
               if (igraph::gsize(g) > 0) {
                    com <- igraph::cluster_louvain(g)
                    n_comm_sim    [ii, k] <- length(igraph::sizes(com))
                    modularity_sim[ii, k] <- igraph::modularity(com)
               } else {
                    n_comm_sim    [ii, k] <- igraph::vcount(g)
                    modularity_sim[ii, k] <- NA_real_
               }
               
               # Triangles (total)
               tri_per_v <- igraph::count_triangles(g)
               triangles_sim[ii, k] <- sum(tri_per_v) / 3
          }
     }
     
     # Observed metrics from y
     observed_density <- observed_transitivity <- observed_assortativity <- numeric(K)
     observed_loc_trans_mean <- observed_deg_mean <- observed_deg_sd <- numeric(K)
     observed_mean_dist <- observed_diameter <- observed_n_comm <- observed_modularity <- numeric(K)
     observed_gcc_frac <- observed_triangles <- numeric(K)
     
     for (k in seq_len(K)) {
          g_obs <- igraph::graph_from_adjacency_matrix(y[, , k], mode = "undirected", diag = FALSE)
          
          # Degree
          deg_o <- igraph::degree(g_obs)
          observed_deg_mean[k] <- mean(deg_o)
          observed_deg_sd  [k] <- stats::sd(deg_o)
          
          # Density / clustering
          observed_density     [k] <- igraph::edge_density(g_obs, loops = FALSE)
          observed_transitivity[k] <- igraph::transitivity(g_obs, type = "global")
          loc_o <- igraph::transitivity(g_obs, type = "localundirected", isolates = "zero")
          observed_loc_trans_mean[k] <- mean(loc_o, na.rm = TRUE)
          
          # Assortativity
          observed_assortativity[k] <- igraph::assortativity_degree(g_obs, directed = FALSE)
          
          # Distances / diameter
          observed_mean_dist[k] <- suppressWarnings(tryCatch(
               igraph::mean_distance(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          observed_diameter[k] <- suppressWarnings(tryCatch(
               igraph::diameter(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          
          # Components
          comps_o <- igraph::components(g_obs)
          observed_gcc_frac[k] <- max(comps_o$csize) / igraph::vcount(g_obs)
          
          # Communities
          if (igraph::gsize(g_obs) > 0) {
               com_o <- igraph::cluster_fast_greedy(g_obs)
               observed_n_comm   [k] <- length(igraph::sizes(com_o))
               observed_modularity[k] <- igraph::modularity(com_o)
          } else {
               observed_n_comm   [k] <- igraph::vcount(g_obs)
               observed_modularity[k] <- NA_real_
          }
          
          # Triangles
          tri_per_v_o <- igraph::count_triangles(g_obs)
          observed_triangles[k] <- sum(tri_per_v_o) / 3
     }
     
     # Summaries
     alpha <- (1 - ci) / 2
     sum_density        <- summ_mat(density_sim,        alpha)
     sum_transitivity   <- summ_mat(transitivity_sim,   alpha)
     sum_assortativity  <- summ_mat(assortativity_sim,  alpha)
     sum_loc_trans_mean <- summ_mat(loc_trans_mean_sim, alpha)
     sum_deg_mean       <- summ_mat(deg_mean_sim,       alpha)
     sum_deg_sd         <- summ_mat(deg_sd_sim,         alpha)
     sum_mean_dist      <- summ_mat(mean_dist_sim,      alpha)
     sum_diameter       <- summ_mat(diameter_sim,       alpha)
     sum_n_comm         <- summ_mat(n_comm_sim,         alpha)
     sum_modularity     <- summ_mat(modularity_sim,     alpha)
     sum_gcc_frac       <- summ_mat(gcc_frac_sim,       alpha)
     sum_triangles      <- summ_mat(triangles_sim,      alpha)
     
     # Plots
     p1  <- plot_metrics(alpha, layer_titles, sum_density,        observed_density,        "Density")
     p2  <- plot_metrics(alpha, layer_titles, sum_transitivity,   observed_transitivity,   "Transitivity")
     p3  <- plot_metrics(alpha, layer_titles, sum_loc_trans_mean, observed_loc_trans_mean, "Mean local clust.")
     p4  <- plot_metrics(alpha, layer_titles, sum_assortativity,  observed_assortativity,  "Assortativity")
     p5  <- plot_metrics(alpha, layer_titles, sum_deg_mean,       observed_deg_mean,       "Mean degree")
     p6  <- plot_metrics(alpha, layer_titles, sum_deg_sd,         observed_deg_sd,         "SD degree")
     p7  <- plot_metrics(alpha, layer_titles, sum_mean_dist,      observed_mean_dist,      "Mean geodesic dist.")
     p8  <- plot_metrics(alpha, layer_titles, sum_diameter,       observed_diameter,       "Diameter")
     p9  <- plot_metrics(alpha, layer_titles, sum_gcc_frac,       observed_gcc_frac,       "Giant fraction")
     p10 <- plot_metrics(alpha, layer_titles, sum_n_comm,         observed_n_comm,         "No. communities")
     p11 <- plot_metrics(alpha, layer_titles, sum_modularity,     observed_modularity,     "Modularity")
     p12 <- plot_metrics(alpha, layer_titles, sum_triangles,      observed_triangles,      "Triangle count")
     
     gridExtra::grid.arrange(
          p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
          ncol = 4
     )
     
     invisible(list(
          sims = list(
               density = density_sim, transitivity = transitivity_sim, assortativity = assortativity_sim,
               loc_clustering = loc_trans_mean_sim, mean_degree = deg_mean_sim, sd_degree = deg_sd_sim,
               mean_distance = mean_dist_sim, diameter = diameter_sim, gcc_fraction = gcc_frac_sim,
               n_communities = n_comm_sim, modularity = modularity_sim, triangles = triangles_sim
          ),
          observed = list(
               density = observed_density, transitivity = observed_transitivity, assortativity = observed_assortativity,
               loc_clustering = observed_loc_trans_mean, mean_degree = observed_deg_mean, sd_degree = observed_deg_sd,
               mean_distance = observed_mean_dist, diameter = observed_diameter, gcc_fraction = observed_gcc_frac,
               n_communities = observed_n_comm, modularity = observed_modularity, triangles = observed_triangles
          )
     ))
}

plot_test_stats_model_2 <- function(samples, 
                                    y,
                                    X,
                                    layer_titles = c("Loudness","Brightness","Tonality","Rhythm"),
                                    ci = 0.95,
                                    thin = 50,
                                    seed = 42) {
     set.seed(seed)
     
     # Shapes
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Storage for simulated metrics (S_used x K)
     density_sim        <- matrix(NA_real_, nrow = S_used, ncol = K)
     transitivity_sim   <- matrix(NA_real_, nrow = S_used, ncol = K)
     assortativity_sim  <- matrix(NA_real_, nrow = S_used, ncol = K)
     loc_trans_mean_sim <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_mean_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_sd_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     mean_dist_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     diameter_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     n_comm_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     modularity_sim     <- matrix(NA_real_, nrow = S_used, ncol = K)
     gcc_frac_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     triangles_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     
     # Simulate graphs & compute metrics
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]
          delta_s <- samples$delta[s, , ]
          beta_s  <- samples$beta [s, , ]   # <-- added for model 2
          zeta_s  <- as.numeric(samples$zeta[s])
          
          # model 2 simulator with covariates and zeta
          y_sim <- simulate_multilayer_network_cpp(mu_s, delta_s, X, beta_s, zeta_s)
          
          for (k in seq_len(K)) {
               g <- igraph::graph_from_adjacency_matrix(y_sim[ , , k], mode = "undirected", diag = FALSE)
               
               # Degree stats
               deg <- igraph::degree(g)
               deg_mean_sim[ii, k] <- mean(deg)
               deg_sd_sim  [ii, k] <- stats::sd(deg)
               
               # Density / clustering
               density_sim     [ii, k] <- igraph::edge_density(g, loops = FALSE)
               transitivity_sim[ii, k] <- igraph::transitivity(g, type = "global")
               loc_vec <- igraph::transitivity(g, type = "localundirected", isolates = "zero")
               loc_trans_mean_sim[ii, k] <- mean(loc_vec, na.rm = TRUE)
               
               # Assortativity
               assortativity_sim[ii, k] <- igraph::assortativity_degree(g, directed = FALSE)
               
               # Distances / diameter (safe)
               mean_dist_sim[ii, k] <- suppressWarnings(tryCatch(
                    igraph::mean_distance(g, directed = FALSE, unconnected = TRUE),
                    error = function(e) NA_real_))
               diameter_sim [ii, k] <- suppressWarnings(tryCatch(
                    igraph::diameter(g, directed = FALSE, unconnected = TRUE),
                    error = function(e) NA_real_))
               
               # Components
               comps <- igraph::components(g)
               gcc_frac_sim[ii, k] <- max(comps$csize) / igraph::vcount(g)
               
               # No. communities and modularity
               if (igraph::gsize(g) > 0) {
                    com <- igraph::cluster_louvain(g)
                    n_comm_sim    [ii, k] <- length(igraph::sizes(com))
                    modularity_sim[ii, k] <- igraph::modularity(com)
               } else {
                    n_comm_sim    [ii, k] <- igraph::vcount(g)
                    modularity_sim[ii, k] <- NA_real_
               }
               
               # Triangles (total)
               tri_per_v <- igraph::count_triangles(g)
               triangles_sim[ii, k] <- sum(tri_per_v) / 3
          }
     }
     
     # Observed metrics from y
     observed_density <- observed_transitivity <- observed_assortativity <- numeric(K)
     observed_loc_trans_mean <- observed_deg_mean <- observed_deg_sd <- numeric(K)
     observed_mean_dist <- observed_diameter <- observed_n_comm <- observed_modularity <- numeric(K)
     observed_gcc_frac <- observed_triangles <- numeric(K)
     
     for (k in seq_len(K)) {
          g_obs <- igraph::graph_from_adjacency_matrix(y[, , k], mode = "undirected", diag = FALSE)
          
          # Degree
          deg_o <- igraph::degree(g_obs)
          observed_deg_mean[k] <- mean(deg_o)
          observed_deg_sd  [k] <- stats::sd(deg_o)
          
          # Density / clustering
          observed_density     [k] <- igraph::edge_density(g_obs, loops = FALSE)
          observed_transitivity[k] <- igraph::transitivity(g_obs, type = "global")
          loc_o <- igraph::transitivity(g_obs, type = "localundirected", isolates = "zero")
          observed_loc_trans_mean[k] <- mean(loc_o, na.rm = TRUE)
          
          # Assortativity
          observed_assortativity[k] <- igraph::assortativity_degree(g_obs, directed = FALSE)
          
          # Distances / diameter
          observed_mean_dist[k] <- suppressWarnings(tryCatch(
               igraph::mean_distance(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          observed_diameter[k] <- suppressWarnings(tryCatch(
               igraph::diameter(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          
          # Components
          comps_o <- igraph::components(g_obs)
          observed_gcc_frac[k] <- max(comps_o$csize) / igraph::vcount(g_obs)
          
          # Communities
          if (igraph::gsize(g_obs) > 0) {
               com_o <- igraph::cluster_fast_greedy(g_obs)
               observed_n_comm   [k] <- length(igraph::sizes(com_o))
               observed_modularity[k] <- igraph::modularity(com_o)
          } else {
               observed_n_comm   [k] <- igraph::vcount(g_obs)
               observed_modularity[k] <- NA_real_
          }
          
          # Triangles
          tri_per_v_o <- igraph::count_triangles(g_obs)
          observed_triangles[k] <- sum(tri_per_v_o) / 3
     }
     
     # Summaries
     alpha <- (1 - ci) / 2
     sum_density        <- summ_mat(density_sim,        alpha)
     sum_transitivity   <- summ_mat(transitivity_sim,   alpha)
     sum_assortativity  <- summ_mat(assortativity_sim,  alpha)
     sum_loc_trans_mean <- summ_mat(loc_trans_mean_sim, alpha)
     sum_deg_mean       <- summ_mat(deg_mean_sim,       alpha)
     sum_deg_sd         <- summ_mat(deg_sd_sim,         alpha)
     sum_mean_dist      <- summ_mat(mean_dist_sim,      alpha)
     sum_diameter       <- summ_mat(diameter_sim,       alpha)
     sum_n_comm         <- summ_mat(n_comm_sim,         alpha)
     sum_modularity     <- summ_mat(modularity_sim,     alpha)
     sum_gcc_frac       <- summ_mat(gcc_frac_sim,       alpha)
     sum_triangles      <- summ_mat(triangles_sim,      alpha)
     
     # Plots
     p1  <- plot_metrics(alpha, layer_titles, sum_density,        observed_density,        "Density")
     p2  <- plot_metrics(alpha, layer_titles, sum_transitivity,   observed_transitivity,   "Transitivity")
     p3  <- plot_metrics(alpha, layer_titles, sum_loc_trans_mean, observed_loc_trans_mean, "Mean local clust.")
     p4  <- plot_metrics(alpha, layer_titles, sum_assortativity,  observed_assortativity,  "Assortativity")
     p5  <- plot_metrics(alpha, layer_titles, sum_deg_mean,       observed_deg_mean,       "Mean degree")
     p6  <- plot_metrics(alpha, layer_titles, sum_deg_sd,         observed_deg_sd,         "SD degree")
     p7  <- plot_metrics(alpha, layer_titles, sum_mean_dist,      observed_mean_dist,      "Mean geodesic dist.")
     p8  <- plot_metrics(alpha, layer_titles, sum_diameter,       observed_diameter,       "Diameter")
     p9  <- plot_metrics(alpha, layer_titles, sum_gcc_frac,       observed_gcc_frac,       "Giant fraction")
     p10 <- plot_metrics(alpha, layer_titles, sum_n_comm,         observed_n_comm,         "No. communities")
     p11 <- plot_metrics(alpha, layer_titles, sum_modularity,     observed_modularity,     "Modularity")
     p12 <- plot_metrics(alpha, layer_titles, sum_triangles,      observed_triangles,      "Triangle count")
     
     gridExtra::grid.arrange(
          p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
          ncol = 4
     )
     
     invisible(list(
          sims = list(
               density = density_sim, transitivity = transitivity_sim, assortativity = assortativity_sim,
               loc_clustering = loc_trans_mean_sim, mean_degree = deg_mean_sim, sd_degree = deg_sd_sim,
               mean_distance = mean_dist_sim, diameter = diameter_sim, gcc_fraction = gcc_frac_sim,
               n_communities = n_comm_sim, modularity = modularity_sim, triangles = triangles_sim
          ),
          observed = list(
               density = observed_density, transitivity = observed_transitivity, assortativity = observed_assortativity,
               loc_clustering = observed_loc_trans_mean, mean_degree = observed_deg_mean, sd_degree = observed_deg_sd,
               mean_distance = observed_mean_dist, diameter = observed_diameter, gcc_fraction = observed_gcc_frac,
               n_communities = observed_n_comm, modularity = observed_modularity, triangles = observed_triangles
          )
     ))
}

plot_test_stats_model_3 <- function(samples, 
                                    y,
                                    X,
                                    layer_titles = c("1","2","3","4"),
                                    ci = 0.95,
                                    thin = 50,
                                    seed = 42) {
     set.seed(seed)
     
     # Shapes
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Storage for simulated metrics (S_used x K)
     density_sim        <- matrix(NA_real_, nrow = S_used, ncol = K)
     transitivity_sim   <- matrix(NA_real_, nrow = S_used, ncol = K)
     assortativity_sim  <- matrix(NA_real_, nrow = S_used, ncol = K)
     loc_trans_mean_sim <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_mean_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_sd_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     mean_dist_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     diameter_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     n_comm_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     modularity_sim     <- matrix(NA_real_, nrow = S_used, ncol = K)
     gcc_frac_sim       <- matrix(NA_real_, nrow = S_used, ncol = K) 
     triangles_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     
     # Simulate graphs & compute metrics
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s     <- samples$mu    [s, ]
          delta_s  <- samples$delta [s, , ]
          beta_s   <- samples$beta  [s, , ]
          U_s      <- samples$U     [s, , ]
          lambda_s <- samples$lambda[s, ]
          zeta_s   <- as.numeric(samples$zeta[s])
          
          y_sim <- simulate_multilayer_network_cpp(mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s)
          
          for (k in seq_len(K)) {
               g <- igraph::graph_from_adjacency_matrix(y_sim[ , , k], mode = "undirected", diag = FALSE)
               
               # Degree stats
               deg <- igraph::degree(g)
               deg_mean_sim[ii, k] <- mean(deg)
               deg_sd_sim  [ii, k] <- stats::sd(deg)
               
               # Density / clustering
               density_sim     [ii, k] <- igraph::edge_density(g, loops = FALSE)
               transitivity_sim[ii, k] <- igraph::transitivity(g, type = "global")
               loc_vec <- igraph::transitivity(g, type = "localundirected", isolates = "zero")
               loc_trans_mean_sim[ii, k] <- mean(loc_vec, na.rm = TRUE)
               
               # Assortativity
               assortativity_sim[ii, k] <- igraph::assortativity_degree(g, directed = FALSE)
               
               # Distances
               mean_dist_sim[ii, k] <- safe_mean_distance(g)
               diameter_sim [ii, k] <- safe_diameter(g)
               
               # Components
               comps <- igraph::components(g)
               gcc_frac_sim[ii, k] <- max(comps$csize) / igraph::vcount(g)
               
               # No. communities and modularity
               # Works best on simple graphs; returns a dendrogram object
               if (igraph::gsize(g) > 0) {
                    com <- igraph::cluster_louvain(g)
                    n_comm_sim    [ii, k] <- length(igraph::sizes(com))
                    modularity_sim[ii, k] <- igraph::modularity(com)
               } else {
                    n_comm_sim    [ii, k] <- igraph::vcount(g)
                    modularity_sim[ii, k] <- NA_real_
               }
               
               # Triangles (total)
               tri_per_v <- igraph::count_triangles(g)
               triangles_sim[ii, k] <- sum(tri_per_v) / 3
          }
     }
     
     # Observed metrics from y
     observed_density <- observed_transitivity <- observed_assortativity <- numeric(K)
     observed_loc_trans_mean <- observed_deg_mean <- observed_deg_sd <- numeric(K)
     observed_mean_dist <- observed_diameter <- observed_n_comm <- observed_modularity <- numeric(K)
     observed_gcc_frac <- observed_triangles <- numeric(K)
     
     for (k in seq_len(K)) {
          g_obs <- igraph::graph_from_adjacency_matrix(y[, , k], mode = "undirected", diag = FALSE)
          
          # Degree
          deg_o <- igraph::degree(g_obs)
          observed_deg_mean[k] <- mean(deg_o)
          observed_deg_sd  [k] <- stats::sd(deg_o)
          
          # Density / clustering
          observed_density     [k] <- igraph::edge_density(g_obs, loops = FALSE)
          observed_transitivity[k] <- igraph::transitivity(g_obs, type = "global")
          loc_o <- igraph::transitivity(g_obs, type = "localundirected", isolates = "zero")
          observed_loc_trans_mean[k] <- mean(loc_o, na.rm = TRUE)
          
          # Assortativity
          observed_assortativity[k] <- igraph::assortativity_degree(g_obs, directed = FALSE)
          
          # Distances / diameter
          observed_mean_dist[k] <- suppressWarnings(tryCatch(
               igraph::mean_distance(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          observed_diameter[k] <- suppressWarnings(tryCatch(
               igraph::diameter(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          
          # Components
          comps_o <- igraph::components(g_obs)
          observed_gcc_frac[k] <- max(comps_o$csize) / igraph::vcount(g_obs)
          
          # Communities
          if (igraph::gsize(g_obs) > 0) {
               com_o <- igraph::cluster_fast_greedy(g_obs)
               observed_n_comm   [k] <- length(igraph::sizes(com_o))
               observed_modularity[k] <- igraph::modularity(com_o)
          } else {
               observed_n_comm   [k] <- igraph::vcount(g_obs)
               observed_modularity[k] <- NA_real_
          }
          
          # Triangles
          tri_per_v_o <- igraph::count_triangles(g_obs)
          observed_triangles[k] <- sum(tri_per_v_o) / 3
     }
     
     # Summaries
     alpha <- (1 - ci) / 2
     sum_density        <- summ_mat(density_sim,        alpha)
     sum_transitivity   <- summ_mat(transitivity_sim,   alpha)
     sum_assortativity  <- summ_mat(assortativity_sim,  alpha)
     sum_loc_trans_mean <- summ_mat(loc_trans_mean_sim, alpha)
     sum_deg_mean       <- summ_mat(deg_mean_sim,       alpha)
     sum_deg_sd         <- summ_mat(deg_sd_sim,         alpha)
     sum_mean_dist      <- summ_mat(mean_dist_sim,      alpha)
     sum_diameter       <- summ_mat(diameter_sim,       alpha)
     sum_n_comm         <- summ_mat(n_comm_sim,         alpha)
     sum_modularity     <- summ_mat(modularity_sim,     alpha)
     sum_gcc_frac       <- summ_mat(gcc_frac_sim,       alpha)
     sum_triangles      <- summ_mat(triangles_sim,      alpha)
     
     # Plots
     p1  <- plot_metrics(alpha, layer_titles, sum_density,        observed_density,        "Density")
     p2  <- plot_metrics(alpha, layer_titles, sum_transitivity,   observed_transitivity,   "Transitivity")
     p3  <- plot_metrics(alpha, layer_titles, sum_loc_trans_mean, observed_loc_trans_mean, "Mean local clust.")
     p4  <- plot_metrics(alpha, layer_titles, sum_assortativity,  observed_assortativity,  "Assortativity")
     p5  <- plot_metrics(alpha, layer_titles, sum_deg_mean,       observed_deg_mean,       "Mean degree")
     p6  <- plot_metrics(alpha, layer_titles, sum_deg_sd,         observed_deg_sd,         "SD degree")
     p7  <- plot_metrics(alpha, layer_titles, sum_mean_dist,      observed_mean_dist,      "Mean geodesic dist.")
     p8  <- plot_metrics(alpha, layer_titles, sum_diameter,       observed_diameter,       "Diameter")
     p9  <- plot_metrics(alpha, layer_titles, sum_gcc_frac,       observed_gcc_frac,       "Giant fraction")
     p10 <- plot_metrics(alpha, layer_titles, sum_n_comm,         observed_n_comm,         "No. communities")
     p11 <- plot_metrics(alpha, layer_titles, sum_modularity,     observed_modularity,     "Modularity")
     p12 <- plot_metrics(alpha, layer_titles, sum_triangles,      observed_triangles,      "Triangle count")
     
     gridExtra::grid.arrange(
          p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
          ncol = 4
     )
     
     invisible(list(
          sims = list(
               density = density_sim, transitivity = transitivity_sim, assortativity = assortativity_sim,
               loc_clustering = loc_trans_mean_sim, mean_degree = deg_mean_sim, sd_degree = deg_sd_sim,
               mean_distance = mean_dist_sim, diameter = diameter_sim, gcc_fraction = gcc_frac_sim,
               n_communities = n_comm_sim, modularity = modularity_sim, triangles = triangles_sim
          ),
          observed = list(
               density = observed_density, transitivity = observed_transitivity, assortativity = observed_assortativity,
               loc_clustering = observed_loc_trans_mean, mean_degree = observed_deg_mean, sd_degree = observed_deg_sd,
               mean_distance = observed_mean_dist, diameter = observed_diameter, gcc_fraction = observed_gcc_frac,
               n_communities = observed_n_comm, modularity = observed_modularity, triangles = observed_triangles
          )
     ))
}

plot_test_stats_model_4 <- function(samples, 
                                    y,
                                    X,
                                    layer_titles = c("1","2","3","4"),
                                    ci = 0.95,
                                    thin = 50,
                                    seed = 42) {
     set.seed(seed)
     
     # Shapes
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Storage for simulated metrics (S_used x K)
     density_sim        <- matrix(NA_real_, nrow = S_used, ncol = K)
     transitivity_sim   <- matrix(NA_real_, nrow = S_used, ncol = K)
     assortativity_sim  <- matrix(NA_real_, nrow = S_used, ncol = K)
     loc_trans_mean_sim <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_mean_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_sd_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     mean_dist_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     diameter_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     n_comm_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     modularity_sim     <- matrix(NA_real_, nrow = S_used, ncol = K)
     gcc_frac_sim       <- matrix(NA_real_, nrow = S_used, ncol = K) 
     triangles_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     
     # Simulate graphs & compute metrics (Model 4)
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s     <- samples$mu    [s, ]
          delta_s  <- samples$delta [s, , ]
          beta_s   <- samples$beta  [s, , ]
          U_s      <- samples$U     [s, , ]
          lambda_s <- samples$lambda[s, ]
          zeta_s   <- as.numeric(samples$zeta[s])
          
          y_sim <- simulate_multilayer_network_cpp(mu_s, delta_s, X, beta_s, U_s, lambda_s, zeta_s)
          
          for (k in seq_len(K)) {
               g <- igraph::graph_from_adjacency_matrix(y_sim[ , , k], mode = "undirected", diag = FALSE)
               
               # Degree stats
               deg <- igraph::degree(g)
               deg_mean_sim[ii, k] <- mean(deg)
               deg_sd_sim  [ii, k] <- stats::sd(deg)
               
               # Density / clustering
               density_sim     [ii, k] <- igraph::edge_density(g, loops = FALSE)
               transitivity_sim[ii, k] <- igraph::transitivity(g, type = "global")
               loc_vec <- igraph::transitivity(g, type = "localundirected", isolates = "zero")
               loc_trans_mean_sim[ii, k] <- mean(loc_vec, na.rm = TRUE)
               
               # Assortativity
               assortativity_sim[ii, k] <- igraph::assortativity_degree(g, directed = FALSE)
               
               # Distances
               mean_dist_sim[ii, k] <- safe_mean_distance(g)
               diameter_sim [ii, k] <- safe_diameter(g)
               
               # Components
               comps <- igraph::components(g)
               gcc_frac_sim[ii, k] <- max(comps$csize) / igraph::vcount(g)
               
               # No. communities and modularity
               if (igraph::gsize(g) > 0) {
                    com <- igraph::cluster_louvain(g)
                    n_comm_sim    [ii, k] <- length(igraph::sizes(com))
                    modularity_sim[ii, k] <- igraph::modularity(com)
               } else {
                    n_comm_sim    [ii, k] <- igraph::vcount(g)
                    modularity_sim[ii, k] <- NA_real_
               }
               
               # Triangles (total)
               tri_per_v <- igraph::count_triangles(g)
               triangles_sim[ii, k] <- sum(tri_per_v) / 3
          }
     }
     
     # Observed metrics from y
     observed_density <- observed_transitivity <- observed_assortativity <- numeric(K)
     observed_loc_trans_mean <- observed_deg_mean <- observed_deg_sd <- numeric(K)
     observed_mean_dist <- observed_diameter <- observed_n_comm <- observed_modularity <- numeric(K)
     observed_gcc_frac <- observed_triangles <- numeric(K)
     
     for (k in seq_len(K)) {
          g_obs <- igraph::graph_from_adjacency_matrix(y[, , k], mode = "undirected", diag = FALSE)
          
          # Degree
          deg_o <- igraph::degree(g_obs)
          observed_deg_mean[k] <- mean(deg_o)
          observed_deg_sd  [k] <- stats::sd(deg_o)
          
          # Density / clustering
          observed_density     [k] <- igraph::edge_density(g_obs, loops = FALSE)
          observed_transitivity[k] <- igraph::transitivity(g_obs, type = "global")
          loc_o <- igraph::transitivity(g_obs, type = "localundirected", isolates = "zero")
          observed_loc_trans_mean[k] <- mean(loc_o, na.rm = TRUE)
          
          # Assortativity
          observed_assortativity[k] <- igraph::assortativity_degree(g_obs, directed = FALSE)
          
          # Distances / diameter
          observed_mean_dist[k] <- suppressWarnings(tryCatch(
               igraph::mean_distance(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          observed_diameter[k] <- suppressWarnings(tryCatch(
               igraph::diameter(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          
          # Components
          comps_o <- igraph::components(g_obs)
          observed_gcc_frac[k] <- max(comps_o$csize) / igraph::vcount(g_obs)
          
          # Communities
          if (igraph::gsize(g_obs) > 0) {
               com_o <- igraph::cluster_fast_greedy(g_obs)
               observed_n_comm   [k] <- length(igraph::sizes(com_o))
               observed_modularity[k] <- igraph::modularity(com_o)
          } else {
               observed_n_comm   [k] <- igraph::vcount(g_obs)
               observed_modularity[k] <- NA_real_
          }
          
          # Triangles
          tri_per_v_o <- igraph::count_triangles(g_obs)
          observed_triangles[k] <- sum(tri_per_v_o) / 3
     }
     
     # Summaries
     alpha <- (1 - ci) / 2
     sum_density        <- summ_mat(density_sim,        alpha)
     sum_transitivity   <- summ_mat(transitivity_sim,   alpha)
     sum_assortativity  <- summ_mat(assortativity_sim,  alpha)
     sum_loc_trans_mean <- summ_mat(loc_trans_mean_sim, alpha)
     sum_deg_mean       <- summ_mat(deg_mean_sim,       alpha)
     sum_deg_sd         <- summ_mat(deg_sd_sim,         alpha)
     sum_mean_dist      <- summ_mat(mean_dist_sim,      alpha)
     sum_diameter       <- summ_mat(diameter_sim,       alpha)
     sum_n_comm         <- summ_mat(n_comm_sim,         alpha)
     sum_modularity     <- summ_mat(modularity_sim,     alpha)
     sum_gcc_frac       <- summ_mat(gcc_frac_sim,       alpha)
     sum_triangles      <- summ_mat(triangles_sim,      alpha)
     
     # Plots
     p1  <- plot_metrics(alpha, layer_titles, sum_density,        observed_density,        "Density")
     p2  <- plot_metrics(alpha, layer_titles, sum_transitivity,   observed_transitivity,   "Transitivity")
     p3  <- plot_metrics(alpha, layer_titles, sum_loc_trans_mean, observed_loc_trans_mean, "Mean local clust.")
     p4  <- plot_metrics(alpha, layer_titles, sum_assortativity,  observed_assortativity,  "Assortativity")
     p5  <- plot_metrics(alpha, layer_titles, sum_deg_mean,       observed_deg_mean,       "Mean degree")
     p6  <- plot_metrics(alpha, layer_titles, sum_deg_sd,         observed_deg_sd,         "SD degree")
     p7  <- plot_metrics(alpha, layer_titles, sum_mean_dist,      observed_mean_dist,      "Mean geodesic dist.")
     p8  <- plot_metrics(alpha, layer_titles, sum_diameter,       observed_diameter,       "Diameter")
     p9  <- plot_metrics(alpha, layer_titles, sum_gcc_frac,       observed_gcc_frac,       "Giant fraction")
     p10 <- plot_metrics(alpha, layer_titles, sum_n_comm,         observed_n_comm,         "No. communities")
     p11 <- plot_metrics(alpha, layer_titles, sum_modularity,     observed_modularity,     "Modularity")
     p12 <- plot_metrics(alpha, layer_titles, sum_triangles,      observed_triangles,      "Triangle count")
     
     gridExtra::grid.arrange(
          p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
          ncol = 4
     )
     
     invisible(list(
          sims = list(
               density = density_sim, transitivity = transitivity_sim, assortativity = assortativity_sim,
               loc_clustering = loc_trans_mean_sim, mean_degree = deg_mean_sim, sd_degree = deg_sd_sim,
               mean_distance = mean_dist_sim, diameter = diameter_sim, gcc_fraction = gcc_frac_sim,
               n_communities = n_comm_sim, modularity = modularity_sim, triangles = triangles_sim
          ),
          observed = list(
               density = observed_density, transitivity = observed_transitivity, assortativity = observed_assortativity,
               loc_clustering = observed_loc_trans_mean, mean_degree = observed_deg_mean, sd_degree = observed_deg_sd,
               mean_distance = observed_mean_dist, diameter = observed_diameter, gcc_fraction = observed_gcc_frac,
               n_communities = observed_n_comm, modularity = observed_modularity, triangles = observed_triangles
          )
     ))
}

plot_test_stats_model_5 <- function(samples, 
                                    y,
                                    X,
                                    layer_titles = c("1","2","3","4"),
                                    ci = 0.95,
                                    thin = 50,
                                    seed = 42) {
     set.seed(seed)
     
     # Shapes
     S_tot <- nrow(samples$mu)
     K     <- ncol(samples$mu)
     n     <- dim(samples$delta)[2]
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Storage for simulated metrics (S_used x K)
     density_sim        <- matrix(NA_real_, nrow = S_used, ncol = K)
     transitivity_sim   <- matrix(NA_real_, nrow = S_used, ncol = K)
     assortativity_sim  <- matrix(NA_real_, nrow = S_used, ncol = K)
     loc_trans_mean_sim <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_mean_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     deg_sd_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     mean_dist_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     diameter_sim       <- matrix(NA_real_, nrow = S_used, ncol = K)
     n_comm_sim         <- matrix(NA_real_, nrow = S_used, ncol = K)
     modularity_sim     <- matrix(NA_real_, nrow = S_used, ncol = K)
     gcc_frac_sim       <- matrix(NA_real_, nrow = S_used, ncol = K) 
     triangles_sim      <- matrix(NA_real_, nrow = S_used, ncol = K)
     
     # Simulate graphs & compute metrics (Model 5)
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          
          mu_s    <- samples$mu   [s, ]        # length K
          delta_s <- samples$delta[s, , ]      # n x K
          beta_s  <- samples$beta [s, , ]      # p x K
          
          # Theta stored as C x C x (K * S); take this sample’s K slices → C x C x K
          Theta_s <- samples$Theta[ , , ((s - 1L) * K + 1L):(s * K)]
          # Layer-specific cluster labels for this draw → n x K (0-based)
          xi_s    <- samples$xi[s, , ]
          
          zeta_s  <- as.numeric(samples$zeta[s])
          
          # Model 5 simulator
          y_sim <- simulate_multilayer_network_cpp(mu_s, delta_s, X, beta_s, Theta_s, xi_s, zeta_s)
          
          for (k in seq_len(K)) {
               g <- igraph::graph_from_adjacency_matrix(y_sim[ , , k], mode = "undirected", diag = FALSE)
               
               # Degree stats
               deg <- igraph::degree(g)
               deg_mean_sim[ii, k] <- mean(deg)
               deg_sd_sim  [ii, k] <- stats::sd(deg)
               
               # Density / clustering
               density_sim     [ii, k] <- igraph::edge_density(g, loops = FALSE)
               transitivity_sim[ii, k] <- igraph::transitivity(g, type = "global")
               loc_vec <- igraph::transitivity(g, type = "localundirected", isolates = "zero")
               loc_trans_mean_sim[ii, k] <- mean(loc_vec, na.rm = TRUE)
               
               # Assortativity
               assortativity_sim[ii, k] <- igraph::assortativity_degree(g, directed = FALSE)
               
               # Distances
               mean_dist_sim[ii, k] <- safe_mean_distance(g)
               diameter_sim [ii, k] <- safe_diameter(g)
               
               # Components
               comps <- igraph::components(g)
               gcc_frac_sim[ii, k] <- max(comps$csize) / igraph::vcount(g)
               
               # No. communities and modularity
               if (igraph::gsize(g) > 0) {
                    com <- igraph::cluster_louvain(g)
                    n_comm_sim    [ii, k] <- length(igraph::sizes(com))
                    modularity_sim[ii, k] <- igraph::modularity(com)
               } else {
                    n_comm_sim    [ii, k] <- igraph::vcount(g)
                    modularity_sim[ii, k] <- NA_real_
               }
               
               # Triangles (total)
               tri_per_v <- igraph::count_triangles(g)
               triangles_sim[ii, k] <- sum(tri_per_v) / 3
          }
     }
     
     # Observed metrics from y
     observed_density <- observed_transitivity <- observed_assortativity <- numeric(K)
     observed_loc_trans_mean <- observed_deg_mean <- observed_deg_sd <- numeric(K)
     observed_mean_dist <- observed_diameter <- observed_n_comm <- observed_modularity <- numeric(K)
     observed_gcc_frac <- observed_triangles <- numeric(K)
     
     for (k in seq_len(K)) {
          g_obs <- igraph::graph_from_adjacency_matrix(y[, , k], mode = "undirected", diag = FALSE)
          
          # Degree
          deg_o <- igraph::degree(g_obs)
          observed_deg_mean[k] <- mean(deg_o)
          observed_deg_sd  [k] <- stats::sd(deg_o)
          
          # Density / clustering
          observed_density     [k] <- igraph::edge_density(g_obs, loops = FALSE)
          observed_transitivity[k] <- igraph::transitivity(g_obs, type = "global")
          loc_o <- igraph::transitivity(g_obs, type = "localundirected", isolates = "zero")
          observed_loc_trans_mean[k] <- mean(loc_o, na.rm = TRUE)
          
          # Assortativity
          observed_assortativity[k] <- igraph::assortativity_degree(g_obs, directed = FALSE)
          
          # Distances / diameter
          observed_mean_dist[k] <- suppressWarnings(tryCatch(
               igraph::mean_distance(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          observed_diameter[k] <- suppressWarnings(tryCatch(
               igraph::diameter(g_obs, directed = FALSE, unconnected = TRUE),
               error = function(e) NA_real_))
          
          # Components
          comps_o <- igraph::components(g_obs)
          observed_gcc_frac[k] <- max(comps_o$csize) / igraph::vcount(g_obs)
          
          # Communities
          if (igraph::gsize(g_obs) > 0) {
               com_o <- igraph::cluster_fast_greedy(g_obs)
               observed_n_comm   [k] <- length(igraph::sizes(com_o))
               observed_modularity[k] <- igraph::modularity(com_o)
          } else {
               observed_n_comm   [k] <- igraph::vcount(g_obs)
               observed_modularity[k] <- NA_real_
          }
          
          # Triangles
          tri_per_v_o <- igraph::count_triangles(g_obs)
          observed_triangles[k] <- sum(tri_per_v_o) / 3
     }
     
     # Summaries
     alpha <- (1 - ci) / 2
     sum_density        <- summ_mat(density_sim,        alpha)
     sum_transitivity   <- summ_mat(transitivity_sim,   alpha)
     sum_assortativity  <- summ_mat(assortativity_sim,  alpha)
     sum_loc_trans_mean <- summ_mat(loc_trans_mean_sim, alpha)
     sum_deg_mean       <- summ_mat(deg_mean_sim,       alpha)
     sum_deg_sd         <- summ_mat(deg_sd_sim,         alpha)
     sum_mean_dist      <- summ_mat(mean_dist_sim,      alpha)
     sum_diameter       <- summ_mat(diameter_sim,       alpha)
     sum_n_comm         <- summ_mat(n_comm_sim,         alpha)
     sum_modularity     <- summ_mat(modularity_sim,     alpha)
     sum_gcc_frac       <- summ_mat(gcc_frac_sim,       alpha)
     sum_triangles      <- summ_mat(triangles_sim,      alpha)
     
     # Plots
     p1  <- plot_metrics(alpha, layer_titles, sum_density,        observed_density,        "Density")
     p2  <- plot_metrics(alpha, layer_titles, sum_transitivity,   observed_transitivity,   "Transitivity")
     p3  <- plot_metrics(alpha, layer_titles, sum_loc_trans_mean, observed_loc_trans_mean, "Mean local clust.")
     p4  <- plot_metrics(alpha, layer_titles, sum_assortativity,  observed_assortativity,  "Assortativity")
     p5  <- plot_metrics(alpha, layer_titles, sum_deg_mean,       observed_deg_mean,       "Mean degree")
     p6  <- plot_metrics(alpha, layer_titles, sum_deg_sd,         observed_deg_sd,         "SD degree")
     p7  <- plot_metrics(alpha, layer_titles, sum_mean_dist,      observed_mean_dist,      "Mean geodesic dist.")
     p8  <- plot_metrics(alpha, layer_titles, sum_diameter,       observed_diameter,       "Diameter")
     p9  <- plot_metrics(alpha, layer_titles, sum_gcc_frac,       observed_gcc_frac,       "Giant fraction")
     p10 <- plot_metrics(alpha, layer_titles, sum_n_comm,         observed_n_comm,         "No. communities")
     p11 <- plot_metrics(alpha, layer_titles, sum_modularity,     observed_modularity,     "Modularity")
     p12 <- plot_metrics(alpha, layer_titles, sum_triangles,      observed_triangles,      "Triangle count")
     
     gridExtra::grid.arrange(
          p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
          ncol = 4
     )
     
     invisible(list(
          sims = list(
               density = density_sim, transitivity = transitivity_sim, assortativity = assortativity_sim,
               loc_clustering = loc_trans_mean_sim, mean_degree = deg_mean_sim, sd_degree = deg_sd_sim,
               mean_distance = mean_dist_sim, diameter = diameter_sim, gcc_fraction = gcc_frac_sim,
               n_communities = n_comm_sim, modularity = modularity_sim, triangles = triangles_sim
          ),
          observed = list(
               density = observed_density, transitivity = observed_transitivity, assortativity = observed_assortativity,
               loc_clustering = observed_loc_trans_mean, mean_degree = observed_deg_mean, sd_degree = observed_deg_sd,
               mean_distance = observed_mean_dist, diameter = observed_diameter, gcc_fraction = observed_gcc_frac,
               n_communities = observed_n_comm, modularity = observed_modularity, triangles = observed_triangles
          )
     ))
}

originality_by_layer <- function(samples,
                                 layer_titles = c("Loudness","Brightness","Tonality","Rhythm"),
                                 thin = 50,
                                 gamma = 0.80) {
     S_tot <- dim(samples$delta)[1]
     n     <- dim(samples$delta)[2]
     K     <- dim(samples$delta)[3]
     
     if (is.null(layer_titles)) layer_titles <- paste0("Layer ", seq_len(K))
     
     s_index <- seq(from = thin, to = S_tot, by = thin)
     S_used  <- length(s_index)
     
     # Per-draw proportion negative per layer
     prop_neg <- matrix(NA_real_, S_used, K)
     # Per-draw indicator that a majority of nodes are negative
     maj_ind  <- matrix(NA_real_, S_used, K)
     
     # Accumulate node-wise negative counts to compute q_i,k = P(delta_i,k <= 0 | data)
     neg_counts <- matrix(0L, n, K)
     
     for (ii in seq_along(s_index)) {
          s <- s_index[ii]
          D <- samples$delta[s, , ]
          
          for (k in 1:K) {
               neg_vec <- (D[,k] <= 0)
               prop    <- mean(neg_vec)
               prop_neg[ii, k] <- prop
               maj_ind [ii, k] <- as.integer(prop >= 0.5)
               
               # accumulate counts
               neg_counts[, k] <- neg_counts[, k] + neg_vec
          }
     }
     
     # Posterior node-wise sign probabilities q_i,k
     q_node <- neg_counts / S_used   # n x K
     
     O_neg_mean_mean <- round(colMeans(prop_neg, na.rm = TRUE), 3)
     O_majority      <- round(colMeans(maj_ind,  na.rm = TRUE), 3)
     O_prev_gamma    <- round(colMeans(q_node >= gamma, na.rm = TRUE), 3)
     
     out <- data.frame(
          layer        = layer_titles,
          O_neg_mean   = O_neg_mean_mean,
          O_majority   = O_majority,
          O_prev_gamma = O_prev_gamma,
          row.names    = NULL,
          check.names  = FALSE
     )
     
     out
}