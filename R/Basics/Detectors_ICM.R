#Basic method of ICM

#icm method
#refer to this article "Inductive conforma martingales for change point detection".
ICM <- function(training_set,
                stream_data,
                non_conformity_measure,
                betting_function,
                th = 1,
                params_bf = list(), #args for betting function
                ...) { #args for NCM
  
  training_set <- sample(training_set)
  m <- length(training_set)
  n <- length(stream_data)
  
  alphas <- numeric(n)
  p_vals <- numeric(n)
  Cn     <- numeric(n)
  C_prev <- 0
  
  alarm<-NA
  for (i in seq_len(n)) {
    xi <- stream_data[i]
    
    alphas[i] <- do.call(non_conformity_measure,
                         c(list(xi = xi, training_set = training_set), list(...)))
    
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    p_hist <- if (i == 1) numeric(0) else p_vals[1:(i - 1)]
    
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p_vals[i], i = i), params_bf))
    
    val <- if (g_i > 0) (C_prev + log(g_i)) else -Inf
    Cn[i] <- max(0, val)
    C_prev <- Cn[i]
    
    if (Cn[i] > th) {
      alarm <- i
      break
    }
  }
  
  if (!is.na(alarm)) {
    alphas <- alphas[1:alarm]
    p_vals <- p_vals[1:alarm]
    Cn     <- Cn    [1:alarm]
  }
  
  list(
    Cn           = Cn,
    p_vals       = p_vals,
    change_point = alarm
  )
}

#ICM multiple
ICM_multi <- function(training_set,
                      stream_data,
                      non_conformity_measure,
                      betting_function,
                      th = 1,
                      params_bf = list(),
                      ...) {
  
  training_set <- sample(training_set)
  n <- length(stream_data)
  
  all_Cn     <- numeric()
  all_p_vals <- numeric()
  change_points <- c()
  
  start_idx <- 1
  C_prev <- 0
  
  while (start_idx <= n) {
    alphas <- numeric(n - start_idx + 1)
    p_vals <- numeric(n - start_idx + 1)
    Cn     <- numeric(n - start_idx + 1)
    
    for (i in seq_along(alphas)) {
      xi <- stream_data[start_idx + i - 1]
      
      alphas[i] <- do.call(non_conformity_measure,
                           c(list(xi = xi, training_set = training_set), list(...)))
      
      greater <- sum(alphas[1:i] > alphas[i])
      equal   <- sum(alphas[1:i] == alphas[i])
      p_vals[i] <- (greater + runif(1) * equal) / i
      
      p_hist <- if (i == 1) numeric(0) else p_vals[1:(i - 1)]
      
      g_i <- do.call(betting_function,
                     c(list(p_values = p_hist, new_p = p_vals[i], i = i), params_bf))
      
      val <- if (g_i > 0) (C_prev + log(g_i)) else -Inf
      Cn[i] <- max(0, val)
      C_prev <- Cn[i]
      
      if (Cn[i] > th) {
        alarm_idx <- start_idx + i - 1
        change_points <- c(change_points, alarm_idx)
        
        all_Cn     <- c(all_Cn, Cn[1:i])
        all_p_vals <- c(all_p_vals, p_vals[1:i])
        
        start_idx <- alarm_idx + 1
        C_prev <- 0
        break
      }
    }
    
    if (start_idx <= n && length(Cn) + start_idx - 1 == n) {
      all_Cn     <- c(all_Cn, Cn)
      all_p_vals <- c(all_p_vals, p_vals)
      break
    }
  }
  
  list(
    Cn            = all_Cn,
    p_vals        = all_p_vals,
    change_points = change_points
  )
}


ICM_multi_adaptive <- function(stream_data,
                               non_conformity_measure,
                               betting_function,
                               th = 1,
                               training_set = NULL,
                               training_size = NULL,
                               m_retrain = NULL,       # re-training size
                               guard_band = 0,         # skip after an alarm
                               shuffle_training = TRUE,
                               params_bf = list(),
                               ...) {
  n <- length(stream_data)
  
  if (!is.null(training_set)) {
    if (shuffle_training) training_set <- sample(training_set)
    m0  <- length(training_set)
    pos <- 1L
  } else {
    if (is.null(training_size) || training_size <= 0)
      stop("Set training size > 0 or define a training_set")
    m0 <- training_size
    if (m0 >= n) stop("training_size cannot be >= length(stream_data).")
    training_set <- stream_data[1:m0]
    if (shuffle_training) training_set <- sample(training_set)
    pos <- m0 + 1L
  }
  mR <- if (is.null(m_retrain)) m0 else m_retrain
  
  all_Cn     <- rep(NA_real_, n)
  all_p_vals <- rep(NA_real_, n)
  
  change_points <- integer(0)
  events <- list()
  
  C_prev <- 0
  i_seg <- 0L
  alphas_seg <- numeric(0)
  p_seg <- numeric(0)
  Cn_seg <- numeric(0)
  seg_start_pos <- pos
  
  while (pos <= n) {
    i_seg <- i_seg + 1L
    xi <- stream_data[pos]
    
    alpha_i <- do.call(non_conformity_measure,
                       c(list(xi = xi, training_set = training_set), list(...)))
    alphas_seg[i_seg] <- alpha_i
    
    greater <- sum(alphas_seg[1:i_seg] >  alpha_i)
    equal   <- sum(alphas_seg[1:i_seg] == alpha_i)
    p_i <- (greater + runif(1) * equal) / i_seg
    p_seg[i_seg] <- p_i
    
    p_hist <- if (i_seg == 1L) numeric(0) else p_seg[1:(i_seg - 1L)]
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p_i, i = i_seg), params_bf))
    
    val <- if (is.finite(g_i) && g_i > 0) (C_prev + log(g_i)) else -Inf
    C_curr <- max(0, val)
    Cn_seg[i_seg] <- C_curr
    C_prev <- C_curr
    
    if (C_curr >= th) {
      alarm_pos <- pos
      change_points <- c(change_points, alarm_pos)
      
      all_Cn[seg_start_pos:alarm_pos]     <- Cn_seg[1:i_seg]
      all_p_vals[seg_start_pos:alarm_pos] <- p_seg[1:i_seg]
      
      train_start <- min(n, alarm_pos + 1L + as.integer(guard_band))
      if (train_start > n) break
      train_end <- min(n, train_start + mR - 1L)
      new_train <- stream_data[train_start:train_end]
      if (length(new_train) < 1L) break
      if (shuffle_training) new_train <- sample(new_train)
      training_set <- new_train
      
      all_Cn[train_start:train_end]     <- NA_real_
      all_p_vals[train_start:train_end] <- NA_real_
      
      events[[length(events) + 1L]] <- data.frame(
        alarm        = alarm_pos,
        seg_start    = seg_start_pos,
        train_start  = train_start,
        train_end    = train_end,
        next_start   = train_end + 1L
      )
      
      pos <- train_end + 1L
      seg_start_pos <- pos
      C_prev <- 0
      i_seg <- 0L
      alphas_seg <- numeric(0)
      p_seg <- numeric(0)
      Cn_seg <- numeric(0)
      next
    }
    
    pos <- pos + 1L
  }
  
  if (i_seg > 0L) {
    end_pos <- min(n, seg_start_pos + i_seg - 1L)
    all_Cn[seg_start_pos:end_pos]     <- Cn_seg[1:(end_pos - seg_start_pos + 1L)]
    all_p_vals[seg_start_pos:end_pos] <- p_seg[1:(end_pos - seg_start_pos + 1L)]
  }
  
  ev <- if (length(events)) {
    do.call(rbind, events)
  } else {
    data.frame(alarm=integer(), seg_start=integer(),
               train_start=integer(), train_end=integer(),
               next_start=integer())
  }
  
  list(
    Cn_aligned     = all_Cn,
    p_vals_aligned = all_p_vals,
    change_points  = change_points,
    events         = ev
  )
}
