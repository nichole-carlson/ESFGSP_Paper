# ----- Model Fitting -----
# Fit LASSO models using pixel and frequency space covariates.
# Split data: 80% training, 20% test.
# Tune lambda via 10-fold cross-validation.
# Optimal lambda: lowest average binomial deviance.

# Evaluation
# Metrics: accuracy and AUC.
# Permutation test (100 iterations) for p-values.
# Calculate mean, std deviation, and significant p-values percentage.

# Load required libraries
library(glmnet)
library(pROC)
library(tictoc)
library(hdi)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
results_data_dir <- file.path(simulation_dir, "results", "data")

cat("Loading simulated data ... \n")
tic()
load(file = file.path(results_data_dir, "simulated_data_240826.RData"))
toc()
cat("Simulated data loaded  \n")

source(file.path(code_dir, "simWrapper.r"))


# ----- Define Functions for Model Fitting -----

# Perform Lasso Regression with Cross-Validation
#
# This function performs Lasso regression on the provided dataset using
# cross-validation to identify the optimal lambda values. It evaluates the
# model's performance using both the minimum lambda value (lambda.min) and the
# lambda value within one standard error of the minimum (lambda.1se).
#
# Arguments:
# - x: A matrix of predictor variables.
# - y: A vector of response variables.
# - p_train: A numeric value between 0 and 1 representing the proportion of
#   data to be used for training.
# - seed: An optional seed for reproducibility.
#
# Return:
# A matrix containing the accuracy and AUC for lambda.min and lambda.1se.

perform_lasso <- function(x, y, p_train, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Split the dataset
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]

  # Perform cross-validation to find the optimal lambda values
  cv_model <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
  lambda_min <- cv_model$lambda.min
  lambda_1se <- cv_model$lambda.1se

  # Function to evaluate model performance
  eval_perf <- function(l) {
    mod <- glmnet(x_train, y_train, alpha = 1, lambda = l, family = "binomial")
    preds <- predict(mod, newx = x_test, type = "response")[, 1]
    preds_bin <- ifelse(preds > 0.5, 1, 0)

    acc <- mean(preds_bin == y_test)
    auc <- pROC::auc(pROC::roc(y_test, preds))

    coefs <- as.vector(coef(mod))[-1]

    res <- list(
      metrics = c(acc = acc, auc = auc),
      coefs = coefs
    )
    return(res)
  }

  # Get performance metrics for both lambda.min and lambda.1se
  perf_min <- eval_perf(lambda_min)
  perf_1se <- eval_perf(lambda_1se)

  results <- list(lambda_min = perf_min, lambda_1se = perf_1se)
  return(results)
}

# Perform Permutation Test for Lasso Regression
#
# This function performs a permutation test for Lasso regression on the provided
# dataset. It estimates the significance of the model coefficients by comparing
# the original coefficients to those obtained from permuted response variables.
#
# Arguments:
# - x: A matrix of predictor variables.
# - y: A vector of response variables.
# - n_perm: An integer representing the number of permutations to perform.
# - seed: An optional seed for reproducibility.
#
# Return:
# A 1 * n_pixel matrix containing the p-values for each coefficient based on
# the permutation test.

perm_lasso <- function(x, y, n_perm, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  cv_lasso <- function(x, y) {
    cv_fit <- cv.glmnet(x, y, alpha = 1)
    coefs <- coef(cv_fit, s = "lambda.min")[-1, 1]
    return(coefs)
  }

  # get original coef estimates
  orig_coefs <- cv_lasso(x, y)

  # perform permutations
  perm_coefs <- matrix(NA, n_perm, length(orig_coefs))

  for (i in seq_len(n_perm)) {
    y_perm <- sample(y)
    perm_coefs[i, ] <- cv_lasso(x, y_perm)
  }

  # calculate p-values
  p_vals <- matrix(colMeans(abs(perm_coefs) > abs(orig_coefs)), nrow = 1)
  return(p_vals)
}

# evaluate p-values using the hdi package
hdi_pval <- function(x, y) {
  p_vals <- matrix(hdi::lasso.proj(x, y)$pval, nrow = 1)
  return(p_vals)
}


# ----- Fit Models -----

p_train <- 0.8
n_perm <- 100

# AUC and accuracy
eval_auc_acc <- function(i, sim_data, p_train, seed) {
  n_samples <- sim_data$meta_data$n_samples
  index_st <- (i - 1) * n_samples + 1
  index_end <- i * n_samples
  x <- sim_data$data$x[index_st:index_end, ]
  x_freq <- sim_data$data$x_freq[index_st:index_end, ]
  y <- sim_data$data$y[index_st:index_end]

  perf_metrics <- perform_lasso(x, y, p_train, seed + i)
  perf_metrics_freq <- perform_lasso(x_freq, y, p_train, seed + i)

  return(list(pixel = perf_metrics, freq = perf_metrics_freq))
}

cat("Calculating AUCs and ACCs on simulated data 1 ...\n")
tic()
sim1_auc_acc_coefs <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_auc_acc(i, sim1_data, p_train, seed),
  list_export = ls(),
  list_package = c("glmnet", "pROC")
)
toc()

cat("Calculating AUCs and ACCs on simulated data 2 ... \n")
tic()
sim2_auc_acc_coefs <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_auc_acc(i, sim2_data, p_train, seed),
  list_export = ls(),
  list_package = c("glmnet", "pROC")
)
toc()

extract_auc_acc <- function(.mat) {
  apply(.mat, 2, function(.ls) {
    lambda_min_metrics <- do.call("rbind", lapply(.ls, function(.lls) {
      .lls$lambda_min$metrics
    }))
    lambda_1se_metrics <- do.call("rbind", lapply(.ls, function(.lls) {
      .lls$lambda_1se$metrics
    }))
    metrics <- cbind(lambda_min_metrics, lambda_1se_metrics)
    rownames(metrics) <- NULL
    colnames(metrics) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")

    return(metrics)
  }, simplify = FALSE)
}

extract_coefs <- function(.mat) {
  apply(.mat, 2, function(.ls) {
    lambda_min_coefs <- do.call("rbind", lapply(.ls, function(.lls) {
      .lls$lambda_min$coefs
    }))
    lambda_1se_coefs <- do.call("rbind", lapply(.ls, function(.lls) {
      .lls$lambda_1se$coefs
    }))

    coefs <- list(
      lambda_min = lambda_min_coefs,
      lambda_1se = lambda_1se_coefs
    )
    return(coefs)
  }, simplify = FALSE)
}

sim1_auc_acc <- extract_auc_acc(sim1_auc_acc_coefs)
sim2_auc_acc <- extract_auc_acc(sim2_auc_acc_coefs)
sim1_coefs <- extract_coefs(sim1_auc_acc_coefs)
sim2_coefs <- extract_coefs(sim2_auc_acc_coefs)


# Permutation test for p-values
eval_pvals <- function(i, sim_data, n_perm, seed) {
  n_samples <- sim_data$meta_data$n_samples
  index_st <- (i - 1) * n_samples + 1
  index_end <- i * n_samples
  x <- sim_data$data$x[index_st:index_end, ]
  x_freq <- sim_data$data$x_freq[index_st:index_end, ]
  y <- sim_data$data$y[index_st:index_end]

  p_vals <- perm_lasso(x, y, n_perm, seed + i)
  p_vals_freq <- perm_lasso(x_freq, y, n_perm, seed + i)

  return(list(pixel = p_vals, freq = p_vals_freq))
}

# hdi package for p-values
eval_pvals_hdi <- function(i, sim_data) {
  n_samples <- sim_data$meta_data$n_samples
  index_st <- (i - 1) * n_samples + 1
  index_end <- i * n_samples
  x <- sim_data$data$x[index_st:index_end, ]
  x_freq <- sim_data$data$x_freq[index_st:index_end, ]
  y <- sim_data$data$y[index_st:index_end]

  p_vals <- hdi_pval(x, y)
  p_vals_freq <- hdi_pval(x_freq, y)

  return(list(pixel = p_vals, freq = p_vals_freq))
}

cat("Calculating p-values on simulated data 1 ...\n")
tic()
sim1_pvals <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_pvals(i, sim1_data, n_perm, seed),
  list_export = ls(),
  list_package = c("glmnet")
)
toc()

cat("Calculating p-values on simulated data 2 ...\n")
tic()
sim2_pvals <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_pvals(i, sim2_data, n_perm, seed),
  list_export = ls(),
  list_package = c("glmnet")
)
toc()

sim1_pvals <- apply(sim1_pvals, 2, function(col) {
  as.data.frame(do.call("rbind", col))
}, simplify = FALSE)

sim2_pvals <- apply(sim2_pvals, 2, function(col) {
  as.data.frame(do.call("rbind", col))
}, simplify = FALSE)

cat("Calculating p-values by hdi on simulated data 1 ...\n")
tic()
sim1_pvals_hdi <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_pvals_hdi(i, sim1_data),
  list_export = ls(),
  list_package = c("hdi")
)
toc()

cat("Calculating p-values by hdi on simulated data 2 ...\n")
tic()
sim2_pvals_hdi <- simWrapper(
  n_sim = n_sim,
  f_sim = function(i) eval_pvals_hdi(i, sim2_data),
  list_export = ls(),
  list_package = c("hdi")
)
toc()

sim1_pvals_hdi <- apply(sim1_pvals_hdi, 2, function(col) {
  as.data.frame(do.call("rbind", col))
}, simplify = FALSE)

sim2_pvals_hdi <- apply(sim2_pvals_hdi, 2, function(col) {
  as.data.frame(do.call("rbind", col))
}, simplify = FALSE)

filename <- paste0("model_metrics_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  sim1_auc_acc, sim2_auc_acc, sim1_coefs, sim2_coefs, sim1_pvals, sim2_pvals,
  sim1_pvals_hdi,  sim2_pvals_hdi,
  file = file.path(results_data_dir, filename)
)
