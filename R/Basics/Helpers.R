pvals_from_alphas <- function(alphas) {
  n <- length(alphas); p <- numeric(n)
  for (i in seq_len(n)) {
    ai <- alphas[i]; prev <- alphas[seq_len(i)]
    gt <- sum(prev > ai); eq <- sum(prev == ai)
    p[i] <- (gt + runif(1) * eq) / i
  }
  p
}

`%||%` <- function(a,b) if (is.null(a)) b else a

alphas_from_ncm <- function(stream, training_set, ncm_fun, k = NULL,
                            ..., shuffle_training = FALSE) {
  stopifnot(length(training_set) > 0)
  stream <- as.numeric(stream)
  training_set <- as.numeric(training_set)
  if (shuffle_training && length(training_set) > 1L) {
    training_set <- sample(training_set)
  }
  if (is.null(k)) k <- 7L
  
  if (identical(ncm_fun, Non_conformity_KNN)) {
    if (!requireNamespace("FNN", quietly = TRUE)) {
      stop("Para KNN vectorizado instala 'FNN' (install.packages('FNN')).")
    }
    Xtr <- matrix(training_set, ncol = 1L)
    Xte <- matrix(stream,       ncol = 1L)
    nn  <- FNN::get.knnx(Xtr, Xte, k = k)$nn.dist  # n x k
    
    return(rowMeans(nn))
  }
  
  if (identical(ncm_fun, Non_conformity_MAD)) {
    med  <- stats::median(training_set)
    s    <- stats::mad(training_set, constant = 1)
    s    <- if (s > 0) s else .Machine$double.eps
    return(abs(stream - med) / s)
  }
  
  if (identical(ncm_fun, Non_conformity_LNR)) {
    dots  <- list(...)
    mu_r  <- dots$mu_r %||% 1
    mu0   <- mean(training_set)
    
    num <- dnorm(stream, mean = mu_r, sd = sqrt(2))
    den <- dnorm(stream, mean = mu0, sd = 1)
    return(num / den)
  }
  
  if (identical(ncm_fun, Non_conformity_IQR)) {
    dots <- list(...)
    probs <- dots$probs %||% c(0.25, 0.50, 0.75)
    qtype <- dots$type  %||% 8L
    c_iqr <- dots$c_iqr %||% 1.349        
    qs    <- as.numeric(stats::quantile(training_set, probs = probs, type = qtype, names = FALSE))
    med   <- qs[2]
    width <- max(qs[3] - qs[1], 1e-12)
    return(abs(stream - med) / (width * c_iqr))
  }
  
  vapply(stream, function(xi) {
    ncm_fun(xi = xi, training_set = training_set, k = k, ...)
  }, numeric(1))
}

icm_C_path_from_p <- function(p, betting_function, params_bf = list()) {
  n <- length(p); C <- numeric(n); C_prev <- 0
  for (i in seq_len(n)) {
    p_hist <- if (i == 1) numeric(0) else p[1:(i - 1)]
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p[i], i = i), params_bf))
    if (!is.finite(g_i) || g_i <= 0) g_i <- 1e-12
    C_prev <- max(0, C_prev + log(g_i))
    C[i] <- C_prev
  }
  C
}

icm_cbf_S_path_from_p <- function(p, bet_fun, W = 100, epsilon = 1.5, params_bf = list()) {
  n <- length(p); S <- numeric(n); S_prev <- 1.0
  step_bf <- function(prev_p, new_p, i) {
    fi <- do.call(bet_fun, c(list(p_values = prev_p, new_p = new_p, i = i), params_bf))
    if (!is.finite(fi) || fi <= 0) fi <- 1e-12
    fi
  }
  for (i in seq_len(n)) {
    if (i == 1) {
      ratio <- 0
    } else {
      w <- min(W, i - 1)
      min_window <- min(S[(i - w):(i - 1)])
      if (min_window <= 0) min_window <- 1e-12
      ratio <- (S_prev - min_window) / min_window
    }
    fi <- if (ratio <= epsilon) 1.0 else {
      prev_p <- if (i == 1) numeric(0) else p[seq_len(i - 1)]
      step_bf(prev_p, p[i], i)
    }
    S_prev <- S_prev * fi
    S[i] <- S_prev
  }
  S
}

first_cross_indices_linear <- function(path, h_vals) {
  ord <- order(h_vals)
  h_sorted <- h_vals[ord]
  idx_sorted <- rep(NA_integer_, length(h_sorted))
  j <- 1L
  for (i in seq_along(path)) {
    while (j <= length(h_sorted) && path[i] >= h_sorted[j]) {
      if (is.na(idx_sorted[j])) idx_sorted[j] <- i
      j <- j + 1L
    }
    if (j > length(h_sorted)) break
  }
  idx <- idx_sorted[order(ord)]
  idx
}

make_stream_mean_shifts <- function(n_stream, theta_vec, mu_levels, sd = 1) {
  stopifnot(length(mu_levels) == length(theta_vec) + 1)
  if (length(theta_vec)) theta_vec <- sort(unique(as.integer(theta_vec)))
  
  seg_starts <- c(1L, theta_vec)
  seg_ends   <- c(if (length(theta_vec)) theta_vec - 1L else integer(0), n_stream)
  
  z <- numeric(n_stream)
  off <- 1L
  for (j in seq_along(mu_levels)) {
    a <- seg_starts[j]; b <- seg_ends[j]
    if (b >= a) z[a:b] <- rnorm(b - a + 1L, mean = mu_levels[j], sd = sd)
  }
  list(stream = z, true_changes = as.integer(theta_vec))
}

match_alarms_to_changes_multi <- function(
    alarms, true_changes, n_stream,
    window_mode = c("abs","frac"),
    window_abs = Inf, window_frac = 1.0
){
  window_mode <- match.arg(window_mode)
  alarms <- sort(unique(as.integer(alarms)))
  cj <- sort(as.integer(true_changes))
  J <- length(cj)
  
  seg_end <- if (J > 1) c(cj[-1] - 1L, n_stream) else n_stream
  seg_len <- if (J > 0) pmax(0L, seg_end - cj + 1L) else integer(0)
  
  if (window_mode == "abs") {
    win_raw <- rep(window_abs, J)
  } else {
    win_raw <- pmax(1, floor(as.numeric(seg_len) * window_frac))
  }
  win_eff <- pmin(win_raw, as.numeric(seg_len))
  
  if (length(alarms) == 0L || J == 0L) {
    return(list(
      alarms_df = tibble(
        alarm = integer(0), is_tp = logical(0),
        change_id = integer(0), delay = numeric(0)
      ),
      per_change = tibble(
        change_id = seq_len(J), seg_len = seg_len,
        tp = 0L, delay = NA_real_, extras_after_tp = 0L, fa_in_segment = 0L
      ),
      tp_count = 0L, fa_count = length(alarms), miss_count = J,
      delays = numeric(0), extras_after_tp = 0L
    ))
  }
  
  used_alarm <- rep(FALSE, length(alarms))
  change_id_of_alarm <- rep(NA_integer_, length(alarms))
  delay_of_alarm <- rep(NA_real_, length(alarms))
  tp_per_change <- integer(J)
  delay_per_change <- rep(NA_real_, J)
  extras_after_tp <- integer(J)
  fa_in_segment <- integer(J)
  
  for (j in seq_len(J)) {
    left_j  <- cj[j]
    right_j <- min(cj[j] + win_eff[j], seg_end[j])
    idx <- which(!used_alarm & alarms >= left_j & alarms <= right_j)
    if (length(idx)) {
      a_idx <- idx[1]
      used_alarm[a_idx] <- TRUE
      change_id_of_alarm[a_idx] <- j
      delay_j <- alarms[a_idx] - cj[j]
      delay_of_alarm[a_idx] <- delay_j
      tp_per_change[j] <- 1L
      delay_per_change[j] <- delay_j
    }
  }
  
  for (j in seq_len(J)) {
    seg_l <- cj[j]
    seg_r <- seg_end[j]
    in_seg <- which(alarms >= seg_l & alarms <= seg_r)
    
    if (tp_per_change[j] == 1L) {
      tp_pos <- alarms[which(used_alarm & change_id_of_alarm == j)][1]
      extras_after_tp[j] <- sum(alarms[in_seg] > tp_pos)
      fa_in_segment[j]   <- sum(!used_alarm[in_seg] & alarms[in_seg] < tp_pos)
    } else {
      extras_after_tp[j] <- 0L
      fa_in_segment[j]   <- length(in_seg)
    }
  }
  
  is_tp_vec <- !is.na(change_id_of_alarm)
  
  bad_delay <- which(!is.na(delay_of_alarm) & delay_of_alarm < 0)
  if (length(bad_delay)) {
    stop(sprintf("Delay negativo en alarmas: %s", paste(alarms[bad_delay], collapse = ",")))
  }
  if (any(tp_per_change == 1L)) {
    ok <- TRUE
    for (j in which(tp_per_change == 1L)) {
      if (delay_per_change[j] > win_eff[j] + 1e-9 ||
          delay_per_change[j] > seg_len[j] + 1e-9) {
        ok <- FALSE
      }
    }
    if (!ok) stop("Se detectó un delay > win_eff o > seg_len (esto no debería pasar).")
  }
  
  alarms_df <- tibble(
    alarm     = alarms,
    is_tp     = is_tp_vec,
    change_id = change_id_of_alarm,
    delay     = delay_of_alarm
  )
  
  per_change <- tibble(
    change_id = seq_len(J),
    seg_len   = seg_len,
    tp        = tp_per_change,
    delay     = delay_per_change,
    extras_after_tp = extras_after_tp,
    fa_in_segment   = fa_in_segment
  )
  
  list(
    alarms_df       = alarms_df,
    per_change      = per_change,
    tp_count        = sum(tp_per_change),
    fa_count        = sum(!is_tp_vec),                # sin doble conteo
    miss_count      = J - sum(tp_per_change),
    delays          = alarms_df$delay[alarms_df$is_tp],
    extras_after_tp = sum(extras_after_tp)
  )
}

# --- DROP-IN: detector sin order() -----------------------------------
.detect_first_alarm_from_alphas <- function(alphas_seg, bet_fun, th, params_bf) {
  # 0) Chequeos de tipo (si falla, te dirá exactamente qué llegó mal)
  if (!is.numeric(alphas_seg)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'alphas_seg' no es numeric (clase: %s)",
                 paste(class(alphas_seg), collapse = "/")))
  }
  if (!is.function(bet_fun)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'bet_fun' no es function (clase: %s)",
                 paste(class(bet_fun), collapse = "/")))
  }
  if (!is.list(params_bf)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'params_bf' no es list (clase: %s)",
                 paste(class(params_bf), collapse = "/")))
  }
  if (!is.numeric(th) || length(th) != 1L || !is.finite(th)) {
    stop(".detect_first_alarm_from_alphas: 'th' debe ser numérico escalar finito")
  }
  
  n <- length(alphas_seg)
  if (n == 0L) return(NA_integer_)
  
  bf_eval <- function(p_hist, new_p, i) {
    do.call(bet_fun, c(list(p_values = p_hist, new_p = new_p, i = i), params_bf))
  }
  
  C_prev <- 0.0
  p_seg  <- numeric(n)
  sorted_vals <- numeric(0L)
  
  for (i in seq_len  (n)) {
    ai  <- alphas_seg[i]
    pos <- findInterval(ai, sorted_vals)
    gt_prev <- length(sorted_vals) - pos
    
    eq_prev <- 0L
    if (length(sorted_vals)) {
      L <- pos;       while (L >= 1L && sorted_vals[L] == ai) { eq_prev <- eq_prev + 1L; L <- L - 1L }
      R <- pos + 1L;  while (R <= length(sorted_vals) && sorted_vals[R] == ai) { eq_prev <- eq_prev + 1L; R <- R + 1L }
    }
    eq <- eq_prev + 1L
    
    u <- runif(1)
    if (u <= 0) u <- .Machine$double.eps
    if (u >= 1) u <- 1 - .Machine$double.eps
    p_i <- (gt_prev + u * eq) / i
    p_seg[i] <- p_i
    
    if (pos == length(sorted_vals)) {
      sorted_vals <- c(sorted_vals, ai)
    } else if (pos == 0L) {
      sorted_vals <- c(ai, sorted_vals)
    } else {
      sorted_vals <- append(sorted_vals, ai, after = pos)
    }
    
    g_i <- bf_eval(if (i == 1L) numeric(0) else p_seg[1:(i-1L)], p_i, i)
    if (!is.finite(g_i) || g_i <= 0) {
      bf_name  <- if (!is.null(params_bf$bf_name)) params_bf$bf_name else "bet_fun"
      ncm_name <- if (!is.null(params_bf$ncm_name)) params_bf$ncm_name else "ncm"
      rng <- if (i > 1L) paste0("[", sprintf("%.6g", range(p_seg[1:(i-1L)], na.rm = TRUE)), "]") else "[]"
      stop(sprintf("g_i<=0 o no finito en i=%d (g_i=%s). BF='%s', NCM='%s', new_p=%.8g, len_hist=%d, range_hist=%s",
                   i, as.character(g_i), bf_name, ncm_name, p_i, i-1L, rng), call. = FALSE)
    }
    C_curr <- max(0, C_prev + log(g_i))
    if (C_curr >= th) return(i)
    C_prev <- C_curr
  }
  
  NA_integer_
}

ICM_multi_fast <- function(training_set,
                           stream_data,
                           non_conformity_measure,
                           betting_function,
                           th = 1,
                           params_bf = list(),
                           k = NULL,
                           alphas_full = NULL,
                           ...) {
  n <- length(stream_data)
  if (is.null(alphas_full)) {
    alphas_full <- alphas_from_ncm(stream_data, training_set,
                                   ncm_fun = non_conformity_measure, k = k, ...)
  }
  change_points <- integer(0); start_idx <- 1L
  while (start_idx <= n) {
    alphas_seg <- alphas_full[start_idx:n]
    tau_rel <- .detect_first_alarm_from_alphas(
      alphas_seg,
      bet_fun   = betting_function,
      th        = th,
      params_bf = params_bf
    )
    if (is.na(tau_rel)) break
    alarm_idx <- start_idx + tau_rel - 1L
    change_points <- c(change_points, alarm_idx)
    start_idx <- alarm_idx + 1L
  }
  list(change_points = change_points)
}

ICM_multi_adaptive_fast <- function(stream_data,
                                    non_conformity_measure,
                                    betting_function,
                                    th = 1,
                                    training_set = NULL,
                                    training_size = NULL,
                                    m_retrain = NULL,
                                    guard_band = 0L,
                                    shuffle_training = TRUE,
                                    params_bf = list(),
                                    k = NULL,
                                    ...) {
  n <- length(stream_data)
  if (!is.numeric(th) || length(th) != 1L || !is.finite(th)) {
    stop("'th' should be finite")
  }
  if (is.null(k)) k <- 7L
  
  if (!is.null(training_set)) {
    if (shuffle_training && length(training_set) > 1L) training_set <- sample(training_set)
    m0  <- length(training_set)
    pos <- 1L
  } else {
    if (is.null(training_size) || training_size <= 0L)
      stop("Set training_size > 0 or set training_set.")
    m0 <- as.integer(training_size)
    if (m0 >= n) stop("training_size cannot be >= length(stream_data).")
    training_set <- stream_data[1:m0]
    if (shuffle_training && m0 > 1L) training_set <- sample(training_set)
    pos <- m0 + 1L
  }
  if (is.null(m_retrain)) m_retrain <- m0
  m_retrain <- as.integer(m_retrain)
  guard_band <- as.integer(guard_band)
  if (m_retrain <= 0L) stop("'m_retrain' has to be > 0")
  
  cap   <- max(8L, min(n, ceiling(n / max(1L, m_retrain + guard_band + 1L))))
  cps   <- integer(cap)
  n_cps <- 0L
  
  while (pos <= n) {
    alphas_seg <- alphas_from_ncm(stream_data[pos:n], training_set,
                                  ncm_fun = non_conformity_measure, k = k, ...)
    
    tau_rel <- .detect_first_alarm_from_alphas(
      alphas_seg,
      bet_fun   = betting_function,
      th        = th,
      params_bf = params_bf
    )
    if (is.na(tau_rel)) break
    
    alarm_pos <- pos + tau_rel - 1L
    
    n_cps <- n_cps + 1L
    if (n_cps > length(cps)) cps <- c(cps, integer(length(cps)))
    cps[n_cps] <- alarm_pos
    
    train_start <- alarm_pos + 1L + guard_band
    if (train_start > n) break
    train_end <- min(n, train_start + m_retrain - 1L)
    if (train_end < train_start) break
    
    new_train <- stream_data[train_start:train_end]
    if (shuffle_training && length(new_train) > 1L) new_train <- sample(new_train)
    training_set <- new_train
    
    pos <- train_end + 1L
  }
  
  if (n_cps == 0L) {
    return(list(change_points_stream = integer(0)))
  } else {
    return(list(change_points_stream = cps[seq_len(n_cps)]))
  }
}


