library(tictoc)

proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"

# Source functions for simulating data and fit LASSO models
source(file.path(proj_dir, "R", "simulate_data.R"))
source(file.path(proj_dir, "R", "fit_models.R"))

# Global parameters for simulations
n_iter <- 500
n_samples <- 1000
n_row <- 16
n_col <- 16
seed <- 42
sparse_level <- 0.1
beta_effect <- 0.1
b_effect <- 0.2
p_train <- 0.8

sim1a_generated_data <- run_pixel_to_freq_simulation(
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  effect_size = beta_effect,
  c_adj = generate_exp_corr_matrix(n_row, n_col), # exp decay adj matrix
  seed = seed
)

sim1a_fitted_model <- run_lasso_pipeline(
  x = sim1a_generated_data$data$x,
  y = sim1a_generated_data$data$y,
  p = p_train,
  lambda = "lambda.min",
  seed = seed
)

sim1b_data <- run_simulation_1b(
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  sparse_level = sparse_level,
  effect_size = b_effect,
  seed = seed
)

sim2a_data <- run_pixel_to_freq_simulation(
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  effect_size = beta_effect,
  c_adj = generate_n_neighbor_matrix(nrow, n_col, d_max = 2), # 2-neighbor adj
  seed = seed
)
filename <- paste0("model_metrics_", format(Sys.Date(), "%y%m%d"), ".RData")
