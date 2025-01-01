# Determine project and data dirs based on the operating system
sys_info <- Sys.info()
if (!is.null(sys_info) && sys_info["sysname"] == "Linux") {
  proj_dir <- "/projects/sren@xsede.org/ESFGSP_Paper"
  data_dir <- "/scratch/alpine/sren@xsede.org/esfgsp"
} else {
  proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"
  data_dir <- file.path(proj_dir, "results")
}

# Source functions for fitting LASSO models
source(file.path(proj_dir, "R", "fit_models.R"))
source(file.path(proj_dir, "R", "simulate_data.R"))
source(file.path(proj_dir, "R", "sim_wrapper.R"))

p_train <- 0.8
seed <- 42
packages <- c("glmnet", "pROC", "tictoc", "hdi", "MASS")


# Define the task function for parallel processing
fit_lasso_task <- function(dat, p_train, lambda, seed) {
  # Fit model on data in the pixel space
  pixel_model_metrics <- run_lasso_pipeline(
    x = dat$x,
    y = dat$y,
    p = p_train,
    lambda = lambda,
    seed = seed
  )

  # Transform estimated coefs to the freq space
  eigen_vectors <- dat$hparams$eigen_vectors
  pixel_model_metrics$coefs_freq <- t(eigen_vectors) %*% pixel_model_metrics$coefs

  # Fit model on data in the freq space
  freq_model_metrics <- run_lasso_pipeline(
    x = dat$x_freq,
    y = dat$y,
    p = p_train,
    lambda = lambda,
    seed = seed
  )

  # Transform estimated coefs to the pixel space
  freq_model_metrics$coefs_pixel <- eigen_vectors %*% freq_model_metrics$coefs

  # Return results
  results <- list(
    pixel = pixel_model_metrics,
    freq = freq_model_metrics
  )

  return(results)
}


# # ---------- Simulation 1 ----------
# # Load simulation 1a list of data and fit LASSO models
# cat("Loading Simulation 1a .RData file ...\n")
# tic()
# load(file = file.path(data_dir, "sim1a_pixel_data_241229.RData"))
# toc()
# 
# cat("Fitting models on Simulation 1a data ... \n")
# tic()
# # Create a list of arguments for parallel execution
# arg_list <- lapply(sim1a_exp_pixel_data, function(dat) {
#   list(dat = dat, p_train = p_train, lambda = "lambda.min", seed = seed)
# })
# # Use parallel_wrapper to fit models in parallel
# sim1a_lasso_results <- parallel_wrapper(
#   task_fn = fit_lasso_task,
#   args_list = arg_list,
#   cores = 12,
#   pkgs = packages
# )
# toc()
# 
# cat("Saving fitted model as a .RData file ...\n")
# tic()
# save(
#   sim1a_lasso_results,
#   file = file.path(
#     data_dir,
#     paste0("sim1a_lasso_results_", format(Sys.Date(), "%y%m%d"), ".RData")
#   )
# )
# toc()


# ---------- Simulation 2 ----------
# Load simulation 2a list of data and fit LASSO models
cat("Loading Simulation 2a .RData file ...\n")
tic()
load(file = file.path(data_dir, "sim2a_pixel_data_241231.RData"))
toc()

cat("Fitting models on Simulation 2a data ... \n")
tic()
# Create a list of arguments for parallel execution
arg_list <- lapply(sim2a_exp_pixel_data, function(dat) {
  list(dat = dat, p_train = p_train, lambda = "lambda.min", seed = seed)
})
# Use parallel_wrapper to fit models in parallel
sim2a_lasso_results <- parallel_wrapper(
  task_fn = fit_lasso_task,
  args_list = arg_list,
  cores = 64,
  pkgs = packages
)
toc()

cat("Saving fitted model as a .RData file ...\n")
tic()
save(
  sim2a_lasso_results,
  file = file.path(
    data_dir,
    paste0("sim2a_lasso_results_", format(Sys.Date(), "%y%m%d"), ".RData")
  )
)
toc()
