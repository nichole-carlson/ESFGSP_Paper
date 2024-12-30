library(tictoc)

proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"

# Source functions for fitting LASSO models
source(file.path(proj_dir, "R", "fit_models.R"))
source(file.path(proj_dir, "R", "simulate_data.R"))

p_train <- 0.8
seed <- 42

# Load simulation 1a list of data and fit LASSO models
cat("Loading Simulation 1a .RData file ...\n")
tic()
load(file = file.path(proj_dir, "results", "sim1a_pixel_data_241229.RData"))
toc()

cat("Fitting models on Simulation 1a data ... \n")
tic()
sim1a_lasso_results <- list()
for (dat in sim1a_exp_pixel_data) {
  # Fit model on data in the pixel space, get AUC, ACC and estimated coefs
  pixel_model_metrics <- run_lasso_pipeline(
    x = dat$x,
    y = dat$y,
    p = p_train,
    lambda = "lambda.min",
    seed = seed
  )

  # Transform estimated coefs to the freq space
  eigen_vectors <- dat$hparams$eigen_vectors
  pixel_model_metrics$coefs_freq <- t(eigen_vectors) %*% pixel_model_metrics$coefs

  sim1a_lasso_results[["pixel"]][[i]] <- pixel_model_metrics

  # Fit model on data in the freq space, get AUC, ACC and estimated coefs
  freq_model_metrics <- run_lasso_pipeline(
    x = dat$x_freq,
    y = dat$y,
    p = p_train,
    lambda = "lambda.min",
    seed = seed
  )

  # Transform estimated coefs to the pixel space
  freq_model_metrics$coefs_pixel <- eigen_vectors %*% freq_model_metrics$coefs

  sim1a_fitted_results[["freq"]][[i]] <- freq_model_metrics
}
toc()

cat("Saving fitted model as a .RData file ...\n")
tic()
save(
  sim1a_lasso_results,
  file = file.path(
    proj_dir, "results",
    paste0("sim1a_lasso_results_", format(Sys.Date(), "%y%m%d"), ".RData")
  )
)
toc()
