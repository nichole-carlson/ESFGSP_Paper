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

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
results_dir <- file.path(simulation_dir, "results", "figures")

source(file.path(code_dir, "01_simulate_data.r"))
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
# A matrix containing the accuracy and AUC for models with lambda.min and lambda.1se.

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

        c(acc = acc, AUC = auc)
    }

    # Get performance metrics for both lambda.min and lambda.1se
    perf_min <- eval_perf(lambda_min)
    perf_1se <- eval_perf(lambda_1se)

    # Combine results into a coherent matrix
    results <- c(perf_min, perf_1se)
    results <- matrix(results, nrow = 1)
    colnames(results) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")

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
        best_lambda <- cv_fit$lambda.min
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

# ----- Simulate Data -----

n_sim <- 100
n_samples <- 1000
img_size <- 16
beta_effect <- 0.1
b_effect <- 0.4
b_sparsity <- 0.1
seed <- 42

tic()
cat("Simulating Data for", n_sim, "iterations ... \n")
simulated_data <- simulate_data(
    n_sim = n_sim,
    n_samples = n_samples,
    img_size = img_size,
    beta_effect = beta_effect,
    b_effect = b_effect,
    b_sparsity = b_sparsity,
    seed = seed
)
toc()

# ----- Fit Models -----
