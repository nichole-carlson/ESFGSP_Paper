library(tictoc)

proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"

# Source functions for fitting LASSO models
source(file.path(proj_dir, "R", "fit_models.R"))
source(file.path(proj_dir, "R", "simulate_data.R"))

p_train <- 0.8
seed <- 42

# Load simulation 1a list of data and fit LASSO models
load(file.path(proj_dir, "results", "sim1a_pixel_data_241229.RData"))

sim1a_fitted_results <- list()
for (dat in sim1a_exp_pixel_data) {
  # Get AUC, ACC, estimated coefs
  model_metrics <- run_lasso_pipeline(
    x = dat$x,
    y = dat$y,
    p = p_train,
    lambda = "lambda.min",
    seed = seed
  )

  # Transform estimated coefs to the freq space
  eigen_vectors <- eigen_decomp(dat$hparams$c_adj)$vectors
  model_metrics$coefs_transformed <- t(eigen_vectors) %*% model_metrics$coefs

  sim1a_fitted_results[[i]] <- model_metrics
}

save(
  sim1a_fitted_results,
  file = file.path(
    proj_dir, "results",
    paste0("sim1a_fitted_results_", format(Sys.Date(), "%y%m%d"), ".RData")
  )
)
