#All the functions here correspond to Betting Functions

#Constant BF
#refer to this araticle "Inductive conforma martingales for change point detection".
Constant_BF <- function(p_values, new_p, i, ...) {
  if (new_p >= 0 && new_p < 0.5) {
    return(1.5)
  } else if (new_p >= 0.5 && new_p <= 1) {
    return(0.5)
  } 
}

#Mixture BF
#refer to this article "Inductive conforma martingales for change point detection".
Mixture_BF <- function(p_values, new_p, i, eps = 1e-12, ...) {
  p  <- as.numeric(new_p)
  p  <- pmin(1 - eps, pmax(eps, p))
  lp <- log1p(p - 1)
  
  if (abs(lp) < 1e-12) {
    return(0.5 + lp/12 + (lp*lp)/48)
  }
  
  num <- p * lp - expm1(lp)
  den <- p * (lp * lp)
  
  val <- num / den
  if (!is.finite(val) || val <= 0) val <- eps
  val
}

#Kernel BF
#refer to this article "Inductive conforma martingales for change point detection".
Kernel_BF <- function(p_values, new_p, i, ...) {
  dots <- list(...)
  L           <- if (!is.null(dots$L))           as.integer(dots$L)           else 100L
  n_grid      <- if (!is.null(dots$n_grid))      as.integer(dots$n_grid)      else 512L
  bw_floor    <- if (!is.null(dots$bw_floor))    as.numeric(dots$bw_floor)    else 1e-3
  min_history <- if (!is.null(dots$min_history)) as.integer(dots$min_history) else 2L
  
  p_hist <- as.numeric(p_values)
  m <- length(p_hist)
  
  if (m < min_history) return(1.0)
  
  if (m > L) p_hist <- p_hist[(m - L + 1L):m]
  
  h0 <- stats::bw.nrd0(p_hist)
  if (!is.finite(h0) || h0 <= 0) h0 <- bw_floor
  bw <- max(h0, bw_floor)
  
  ext <- c(-p_hist, p_hist, 2 - p_hist)
  
  dens <- try(stats::density(x = ext, kernel = "gaussian", bw = bw,
                             n = n_grid, from = -0.25, to = 1.25),
              silent = TRUE)
  if (inherits(dens, "try-error")) return(1.0)
  
  in01 <- dens$x >= 0 & dens$x <= 1
  x01  <- dens$x[in01]; y01 <- dens$y[in01]
  
  delta <- x01[2] - x01[1]
  area  <- sum(y01) * delta
  if (is.finite(area) && area > 0) {
    y01 <- y01 / area
  } else {
    y01 <- rep(1, length(x01))
    y01 <- y01 / (sum(y01) * delta)
  }
  
  p <- min(1, max(0, as.numeric(new_p)))
  g <- approx(x01, y01, xout = p, rule = 2)$y
  if (!is.finite(g) || g <= 0) g <- 1e-12
  g
}

#KDE BF
KDE <- function(p_values, n_grid = 512) {
  p_values <- as.numeric(p_values)
  n        <- length(p_values)
  
  # Bandwidth: regla nrd0 + fallback
  h0 <- stats::bw.nrd0(p_values)
  if (!is.finite(h0) || h0 <= 0) {
    h0 <- max(diff(range(p_values)) * n^(-1/5), 1e-3)
  }
  
  # Reflect samples to mitigate bordes
  extended <- c(-p_values, p_values, 2 - p_values)
  
  # Ajustar KDE en [-1.5, 2.5] y luego recortar a [0,1]
  kde <- stats::density(
    x      = extended,
    kernel = "gaussian",
    bw     = h0,
    n      = n_grid,
    from   = -1.5,
    to     = 2.5
  )
  in01  <- (kde$x >= 0 & kde$x <= 1)
  x01   <- kde$x[in01]
  y01   <- kde$y[in01]
  
  # Normalizar para que ∫₀¹g(p)dp = 1
  delta <- x01[2] - x01[1]
  area  <- sum(y01) * delta
  if (area > 0) {
    y01 <- y01 / area
  } else {
    y01 <- rep(1, length(y01))
    y01 <- y01 / (sum(y01) * delta)
  }
  
  # Devuelve una función estática g(p)
  function(p) {
    p <- as.numeric(p)
    out <- rep(0, length(p))
    inside <- which(p >= 0 & p <= 1)
    if (length(inside)>0) {
      out[inside] <- approx(x01, y01, xout = p[inside], rule = 2)$y
    }
    out
  }
}

#PRECOMPUTED KDE BF
#refer to this article "Inductive conforma martingales for change point detection".

Precomputed_KDE_BF <- function(training_set,
                               calibration_data,
                               non_conformity_measure,
                               k = 1,
                               n_grid = 512) {
  m <- length(calibration_data)
  alphas  <- numeric(m)
  pvals   <- numeric(m)
  
  # 1) Calcular α₁…αₘ
  for (i in seq_len(m)) {
    alphas[i] <- non_conformity_measure(calibration_data[i], training_set, k)
  }
  
  # 2) Calcular p₁…pₘ de forma secuencial
  for (i in seq_len(m)) {
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    pvals[i] <- (greater + runif(1) * equal) / i
  }
  
  # 3) Ajustar KDE
  g_p <- KDE(pvals, n_grid = n_grid)
  
  # 4) Retornar función que se ajusta al framework esperado
  function(p_values, new_p, i, ...) {
    g_p(new_p)
  }
}

#HISTOGRAM BF
#refer to this article "A histogram based betting function for conformal martingales".
histogram_betting_function <- function(p_values, new_p, num_bins = 2L, ...) {
  # Asegurar tipos y dominio
  k <- as.integer(num_bins)
  if (k < 1L) stop("k debe ser >= 1")
  clamp01 <- function(x) pmin(1, pmax(0, as.numeric(x)))
  p_prev <- clamp01(p_values)
  p_new  <- clamp01(new_p)
  
  # Si no hay historial aún, no apostar (regla del paper)
  n_prev <- length(p_prev)
  if (n_prev == 0L) return(1)
  
  # Asignación de bins: B1=[0,1/k), ..., Bk=[(k-1)/k, 1]
  bin_of <- function(p) pmax(1L, pmin(k, ceiling(p * k)))
  
  counts <- tabulate(bin_of(p_prev), nbins = k)
  j <- bin_of(p_new)
  
  # Si el bin de p_new está vacío, usar 1 para no anular la martingala (regla del paper)
  if (counts[j] == 0L) return(1)
  
  # Ecuación (4): densidad por histograma
  (counts[j] * k) / n_prev
}


