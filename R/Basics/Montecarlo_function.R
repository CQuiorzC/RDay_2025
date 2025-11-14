library(future.apply)
library(data.table)

source("Basics/Non_conformity_measures.R")
source("Basics/Betting_functions.R")
source("Basics/Helpers.R")

montecarlo_ICM <- function(n_sim        = 200,
                           h_vals       = seq(1, 6, 0.5),
                           theta_stream = 100,
                           mu1          = 1,
                           m            = 200,
                           ncm_fun      = Non_conformity_KNN,
                           bet_fun      = Constant_BF,
                           k            = NULL,
                           params_bf    = list(),
                           n_stream     = 1000) {
  
  if (is.null(k)) k <- 7
  pos_change <- theta_stream
  
  simulate_once <- function(sim_id) {
    training_set <- rnorm(m, 0, 1)
    stream <- c(rnorm(theta_stream - 1, 0, 1),
                rnorm(n_stream - theta_stream + 1, mu1, 1))
    
    alphas <- alphas_from_ncm(stream, training_set, ncm_fun, k = k)
    p <- pvals_from_alphas(alphas)
    
    S <- icm_C_path_from_p(p, bet_fun, params_bf = params_bf)
    
    tau <- first_cross_indices_linear(S, h_vals)
    
    fa <- as.integer(!is.na(tau) & tau < theta_stream)
    dl <- ifelse(is.na(tau) | tau < theta_stream, NA_real_, tau - theta_stream)
    
    list(fa = fa, delay = dl, tau = tau)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  TAUS <- do.call(rbind, lapply(sims, `[[`, "tau"))
  
  taus_wide <- as.data.frame(TAUS)
  colnames(taus_wide) <- paste0("h_", seq_along(h_vals))
  taus_wide$sim_id <- seq_len(nrow(taus_wide))
  
  taus_long <- tidyr::pivot_longer(
    taus_wide,
    cols = dplyr::starts_with("h_"),
    names_to = "h_idx",
    values_to = "tau"
  ) |>
    dplyr::mutate(
      threshold_idx = as.integer(sub("h_", "", h_idx)),
      threshold = h_vals[threshold_idx]
    ) |>
    dplyr::select(sim_id, threshold, tau)
  
  summary_df <- data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
  
  return(list(summary = summary_df, taus = taus_long))
}

montecarlo_ICM_MULTI <- function(
    n_sim        = 200,
    h_vals       = seq(1, 6, 0.5),
    m            = 200,
    n_stream     = 1000,
    theta_vec    = c(300, 700),
    mu_levels    = c(0, 1.5, 0),
    ncm_fun      = Non_conformity_KNN,
    bet_fun      = Constant_BF,
    params_bf    = list(),
    k            = NULL,
    window_mode  = c("abs","frac"),
    window_abs   = Inf,
    window_frac  = 1.0
){
  if (is.null(k)) k <- 7
  window_mode <- match.arg(window_mode)
  J <- length(theta_vec)
  
  one_sim <- function(sim_id){
    training_set <- rnorm(m, 0, 1)
    gen <- make_stream_mean_shifts(n_stream, theta_vec, mu_levels, sd = 1)
    
    alphas_full <- alphas_from_ncm(gen$stream, training_set,
                                   ncm_fun = ncm_fun, k = k)
    
    per_change_list <- vector("list", length(h_vals))
    alarms_list     <- vector("list", length(h_vals))
    stats_mat <- matrix(0, nrow = length(h_vals), ncol = 4)
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      
      out <- ICM_multi_fast(
        training_set           = training_set,
        stream_data            = gen$stream,
        non_conformity_measure = ncm_fun,
        betting_function       = bet_fun,
        th                     = h,
        params_bf              = params_bf,
        k                      = k,
        alphas_full            = alphas_full
      )
      
      eva <- match_alarms_to_changes_multi(
        alarms       = as.integer(out$change_points),
        true_changes = gen$true_changes,
        n_stream     = n_stream,
        window_mode  = window_mode,
        window_abs   = window_abs,
        window_frac  = window_frac
      )
      stats_mat[j,] <- c(eva$tp_count, eva$miss_count, eva$fa_count, eva$extras_after_tp)
      per_change_list[[j]] <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$per_change)
      alarms_list[[j]]     <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$alarms_df)
    }
    
    list(
      stats_mat  = stats_mat,
      per_change = dplyr::bind_rows(per_change_list),
      alarms     = dplyr::bind_rows(alarms_list)
    )
  }
  
  sims <- future.apply::future_lapply(seq_len(n_sim), one_sim, future.seed = TRUE)
  
  stats_arr <- simplify2array(lapply(sims, `[[`, "stats_mat"))
  stats_sum <- apply(stats_arr, c(1,2), sum)
  colnames(stats_sum) <- c("TP","MISS","FA","OSEG")
  TP   <- stats_sum[,"TP"]
  MISS <- stats_sum[,"MISS"]
  FA   <- stats_sum[,"FA"]
  OSEG <- stats_sum[,"OSEG"]
  
  per_change_all <- rbindlist(lapply(sims, `[[`, "per_change"))
  per_change_summ <- per_change_all[, .(
    seg_len          = mean(seg_len),
    recall_j         = mean(tp),
    miss_rate_j      = 1 - mean(tp),
    mean_delay_j     = mean(delay, na.rm = TRUE),
    median_delay_j   = median(delay, na.rm = TRUE),
    iqr_delay_j      = IQR(delay, na.rm = TRUE),
    norm_delay_j     = mean(delay / seg_len, na.rm = TRUE),
    extras_after_tp_j= mean(extras_after_tp),
    fa_in_segment_j  = mean(fa_in_segment)
  ),
  by = .(threshold, change_id)
  ]
  
  delay_by_h <- per_change_all[, .(mean_delay = mean(delay, na.rm = TRUE)), by = .(threshold)]
  
  dt_h <- data.table(threshold = h_vals)
  dt_h[, recall := TP / (J * n_sim)]
  dt_h[, miss_rate := 1 - recall]
  dt_h[, precision := ifelse(TP + FA > 0, TP / (TP + FA), NA_real_)]
  dt_h[, fa_per_1000 := FA / (n_stream * n_sim) * 1000]
  dt_h[, alarms_per_1000 := (TP + FA) / (n_stream * n_sim) * 1000]
  dt_h[, overseg_ratio := ifelse(TP > 0, OSEG / TP, NA_real_)]
  
  dt_h <- merge(dt_h, as.data.table(delay_by_h), by = "threshold", all.x = TRUE)
  dt_h[, log_delay := log10(1 + mean_delay)]
  
  list(
    summary    = as_tibble(dt_h),
    per_change = as_tibble(per_change_summ),
    alarms     = as_tibble(rbindlist(lapply(sims, `[[`, "alarms")))
  )
}

montecarlo_ICM_MULTI_ADAPTIVE <- function(
    n_sim        = 200,
    h_vals       = seq(1, 6, 0.5),
    m            = 200,
    n_stream     = 1000,
    theta_vec    = c(300, 700),
    mu_levels    = c(0, 1.5, 0),
    ncm_fun      = Non_conformity_KNN,
    bet_fun      = Constant_BF,
    params_bf    = list(),
    k            = NULL,
    m_retrain    = NULL,
    guard_band   = 0L,
    window_mode  = c("abs","frac"),
    window_abs   = Inf,
    window_frac  = 1.0
){
  if (is.null(k)) k <- 7
  if (is.null(m_retrain)) m_retrain <- m
  window_mode <- match.arg(window_mode)
  J <- length(theta_vec)
  
  one_sim <- function(sim_id){
    gen <- make_stream_mean_shifts(n_stream, theta_vec, mu_levels, sd = 1)
    training_set0 <- rnorm(m, 0, 1)
    
    per_change_list <- vector("list", length(h_vals))
    alarms_list     <- vector("list", length(h_vals))
    stats_mat <- matrix(0, nrow = length(h_vals), ncol = 4)
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      out <- ICM_multi_adaptive_fast(
        stream_data            = gen$stream,
        non_conformity_measure = ncm_fun,
        betting_function       = bet_fun,
        th                     = h,
        training_set           = training_set0,
        training_size          = NULL,
        m_retrain              = m_retrain,
        guard_band             = guard_band,
        shuffle_training       = TRUE,
        params_bf              = params_bf,
        k                      = k
      )
      eva <- match_alarms_to_changes_multi(
        alarms       = as.integer(out$change_points_stream),
        true_changes = gen$true_changes,
        n_stream     = n_stream,
        window_mode  = window_mode,
        window_abs   = window_abs,
        window_frac  = window_frac
      )
      stats_mat[j,] <- c(eva$tp_count, eva$miss_count, eva$fa_count, eva$extras_after_tp)
      per_change_list[[j]] <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$per_change)
      alarms_list[[j]]     <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$alarms_df)
    }
    
    list(
      stats_mat = stats_mat,
      per_change = dplyr::bind_rows(per_change_list),
      alarms     = dplyr::bind_rows(alarms_list)
    )
  }
  
  sims <- future.apply::future_lapply(seq_len(n_sim), one_sim, future.seed = TRUE)
  
  stats_arr <- simplify2array(lapply(sims, `[[`, "stats_mat"))
  stats_sum <- apply(stats_arr, c(1,2), sum)
  colnames(stats_sum) <- c("TP","MISS","FA","OSEG")
  TP   <- stats_sum[,"TP"]; MISS <- stats_sum[,"MISS"]; FA <- stats_sum[,"FA"]; OSEG <- stats_sum[,"OSEG"]
  
  per_change_all <- rbindlist(lapply(sims, `[[`, "per_change"))
  per_change_summ <- per_change_all[, .(
    seg_len          = mean(seg_len),
    recall_j         = mean(tp),
    miss_rate_j      = 1 - mean(tp),
    mean_delay_j     = mean(delay, na.rm = TRUE),
    median_delay_j   = median(delay, na.rm = TRUE),
    iqr_delay_j      = IQR(delay, na.rm = TRUE),
    norm_delay_j     = mean(delay / seg_len, na.rm = TRUE),
    extras_after_tp_j= mean(extras_after_tp),
    fa_in_segment_j  = mean(fa_in_segment)
  ),
  by = .(threshold, change_id)
  ]
  delay_by_h <- per_change_all[, .(mean_delay = mean(delay, na.rm = TRUE)), by = .(threshold)]
  
  dt_h <- data.table(threshold = h_vals)
  dt_h[, recall := TP / (J * n_sim)]
  dt_h[, miss_rate := 1 - recall]
  dt_h[, precision := ifelse(TP + FA > 0, TP / (TP + FA), NA_real_)]
  dt_h[, fa_per_1000 := FA / (n_stream * n_sim) * 1000]
  dt_h[, alarms_per_1000 := (TP + FA) / (n_stream * n_sim) * 1000]
  dt_h[, overseg_ratio := ifelse(TP > 0, OSEG / TP, NA_real_)]
  
  dt_h <- merge(dt_h, as.data.table(delay_by_h), by = "threshold", all.x = TRUE)
  dt_h[, log_delay := log10(1 + mean_delay)]
  
  list(
    summary    = as_tibble(dt_h),
    per_change = as_tibble(per_change_summ),
    alarms     = as_tibble(rbindlist(lapply(sims, `[[`, "alarms")))
  )
}