source("Basics/Non_conformity_measures.R")
source("Basics/Betting_functions.R")

set.seed(2025)
m0          <- 200
n_calib     <- 1000
theta_cal <- 500
mu1_cal     <- 1
k0          <- 7

train_for_kde <- rnorm(m0, 0, 1)
calib_stream <- numeric(n_calib)
for (i in seq_len(n_calib)) {
  if (i < theta_cal) calib_stream[i] <- rnorm(1, 0, 1) else calib_stream[i] <- rnorm(1, mu1_cal, 1)
}

kde_bf_fixed <- Precomputed_KDE_BF(
  training_set           = train_for_kde,
  calibration_data       = calib_stream,
  non_conformity_measure = Non_conformity_KNN,
  k                      = k0,
  n_grid                 = 512
)

saveRDS(kde_bf_fixed, file = "Data/kde_bf_fixed.rds")