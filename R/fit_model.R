# -------------------------------------------------------
# This script implements a LASSO-based modeling pipeline.
#
# It splits input data (x, y) into training and test sets.
# Cross-validation is performed on the training set to select
# the optimal lambda (either lambda.min or lambda.1se).
#
# The selected model is then evaluated on the test set,
# with performance metrics (accuracy, AUC) and HDI-adjusted
# p-values returned.
#
# A wrapper function is also provided to support simulations.
# It assumes x and y are from pixel space and allows fitting
# either in pixel space or frequency space (if an eigenvector matrix is given).
# -------------------------------------------------------

# List of required packages
packages <- c("glmnet", "pROC", "hdi")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load packages
sapply(packages, require, character.only = TRUE)


# Splits the dataset into training and testing sets.
#
# Args:
#  x: A data frame or matrix with observations as rows and features as columns.
#  y: A vector of target values with the same number of observations as 'x'.
#  p_train: A float value (0-1) specifying the proportion of data for training.
#  seed: An optional integer for reproducibility.
#
# Returns:
#  train: A list with 'x' and 'y'.
#  test: A list with 'x' and 'y'.

train_test_split <- function(x, y, p_train, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check whether x and y have the same # of obs
  if (nrow(x) != length(y)) {
    stop("Check that 'x' and 'y' have the same number of observations.")
  }

  # Calculate the train data size
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  # Split x and y into train and test
  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]

  list(
    train = list(x = x_train, y = y_train),
    test = list(x = x_test, y = y_test)
  )
}


# Fits a Lasso model using cross-validation and refits with the chosen lambda.
#
# Args:
#   x: Predictor matrix (rows: observations, columns: features).
#   y: Response vector, same length as rows in 'x'.
#   lambda: Lambda selection, either "lambda.min" or "lambda.1se".
#
# Returns:
#   A refitted glmnet model object.

fit_lasso_cv <- function(x, y, lambda = c("lambda.min", "lambda.1se")) {
  # Match the 'lambda' argement to ensure it's one of the allowed options
  lambda <- match.arg(lambda)

  if (!all(y %in% c(0, 1))) {
    stop("fit_lasso_cv expects a binary response vector with values 0 and 1.")
  }

  # Perform cross-validation to find the optimal lambda value
  cv_model <- glmnet::cv.glmnet(x, y, alpha = 1, family = "binomial")
  lambda_to_use <- cv_model[[lambda]]

  # Return the refitted model with the lambda chosen
  glmnet::glmnet(
    x, y,
    alpha = 1, lambda = lambda_to_use, family = "binomial"
  )
}


# Computes the threshold that maximizes Youden's J statistic.
#
# Args:
#   labels: Binary vector of true class labels (0/1).
#   probs: Numeric vector of predicted probabilities (same length as labels).
#
# Returns:
#   The cutoff value maximizing Youden’s J = sensitivity + specificity - 1

youden_cutoff <- function(labels, probs) {
  stopifnot(length(labels) == length(probs))

  thresholds <- sort(unique(probs), decreasing = TRUE)
  youdens <- numeric(length(thresholds))

  for (i in seq_along(thresholds)) {
    pred <- as.numeric(probs >= thresholds[i])
    tp <- sum(pred == 1 & labels == 1)
    fn <- sum(pred == 0 & labels == 1)
    tn <- sum(pred == 0 & labels == 0)
    fp <- sum(pred == 1 & labels == 0)

    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0

    youdens[i] <- sens + spec - 1
  }

  best_idx <- which.max(youdens)
  thresholds[best_idx]
}


# Evaluates a Lasso model's performance and extracts estimated coefficients.
#
# Args:
#   model: A fitted glmnet model object.
#   x: Predictor matrix (rows: observations, columns: features).
#   y: Response vector, same length as rows in 'x'.
#   threshold: Threshold for classification. If NULL, use Youden's index.
#
# Returns:
#   A list with:
#     - auc: Area under the ROC curve.
#     - acc: Classification accuracy.
#     - coefs: Estimated coefficients, excluding the intercept.

evaluate_lasso <- function(model, x, y, threshold = NULL) {
  # Check whether x and y have the same # of obs
  if (nrow(x) != length(y)) {
    stop("The number of rows in 'x' does not match the length of 'y'.")
  }

  # Calculate predicted probabilities
  pred_probas <- predict(model, newx = x, type = "response")[, 1]

  # Decide cutoff
  if (is.null(threshold)) {
    threshold <- youden_cutoff(y, pred_probas)
  }

  # Calculate predicted labels
  pred_labels <- ifelse(pred_probas > threshold, 1, 0)

  # Evaluate by accuracy and AUC
  acc <- mean(pred_labels == y)
  auc <- pROC::auc(pROC::roc(y, pred_probas))

  # Extract estimated coefficients
  coefs <- as.vector(coef(model))[-1]

  list(
    auc = auc,
    acc = acc,
    coefs = coefs
  )
}


# Executes a complete lasso analysis pipeline.
#
# Splits data into training and testing sets, performs cross-validation to
# select lambda, refits the model, and evaluates performance on the test set.
#
# Args:
#   x: Predictor matrix.
#   y: Response vector.
#   p_train: Proportion of data for training. Default 0.8.
#   lambda: Lambda selection ("lambda.min" or "lambda.1se")
#   threshold (optional): Classification threshold for binary outcome.
#   seed (optional): Seed for reproducibility.
#
# Returns:
#   A list with:
#     - auc: Area under the ROC curve.
#     - acc: Classification accuracy.
#     - coefs: Model coefficients.
#     - pvals: hdi adjusted p-values.

fit_evaluate_lasso <- function(x,
                               y,
                               p_train = 0.8,
                               lambda = c("lambda.min", "lambda.1se"),
                               threshold = NULL,
                               seed = NULL) {
  lambda <- match.arg(lambda)

  split <- train_test_split(x, y, p_train, seed)
  x_train <- split$train$x
  y_train <- split$train$y
  x_test <- split$test$x
  y_test <- split$test$y

  # Fit model with cross-validation
  model <- fit_lasso_cv(x_train, y_train, lambda)

  # Evaluate model on the test split
  evaluation <- evaluate_lasso(model, x_test, y_test, threshold)

  # Add p-values calculated by hdi
  evaluation$pvals <- as.vector(hdi::lasso.proj(x, y)$pval)

  return(evaluation)
}


# Fit a binary classification Lasso model in either pixel or frequency space.
#
# Args:
#   x: n x p matrix. Design matrix in PIXEL space
#   y: Binary outcome vector of length n.
#   in_pixel_space: Logical. If TRUE, fit model directly on x.
#                   If FALSE, x is transformed using matrix e before fitting.
#   e: p x p transformation matrix (required if in_pixel_space = FALSE).
#      Used to map pixel space to frequency space.
#   ...: Additional arguments passed to fit_evaluate_lasso().
#
# Returns:
#   A list with:
#     - coefs: Estimated coefficients (excluding intercept).
#     - pvals: p-values from HDI.
#     - auc: AUC on the test set.
#     - acc: Accuracy on the test set.
lasso_pixel_or_freq <- function(x, y, in_pixel_space = TRUE, e = NULL, ...) {
  if (!in_pixel_space) {
    if (is.null(e)) {
      stop(
        "Transformation matrix e is required when not in pixel space."
      )
    }
    x <- x %*% e # map pixel tp freq
  }

  fit_evaluate_lasso(x, y, ...)
}
