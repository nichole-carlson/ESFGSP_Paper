library(glmnet)
library(pROC)
library(tictoc)
library(hdi)


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
    stop("The number of rows in 'x' does not match the length of 'y'.")
  }

  # Calculate the train data size
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  # Split x and y into train and test
  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]

  results <- list(
    train = list(x = x_train, y = y_train),
    test = list(x = x_test, y = y_test)
  )

  return(results)
}


# Fits a Lasso model using cross-validation and refits with the chosen lambda.
#
# Args:
#   x: Predictor matrix (rows: observations, columns: features).
#   y: Response vector, same length as rows in 'x'.
#   lambda: Lambda selection, either "lambda.min" or "lambda.1se".
#   family: Model type ("binomial", "gaussian", "poisson"). Default: "binomial".
#
# Returns:
#   A refitted glmnet model object.

cv_lasso <- function(
    x, y, lambda = c("lambda.min", "lambda.1se"), family = "binomial") {
  # Match the 'lambda' argement to ensure it's one of the allowed options
  lambda <- match.arg(lambda)

  # Perform cross-validation to find the optimal lambda value
  cv_model <- glmnet::cv.glmnet(x, y, alpha = 1, family = family)
  lambda_to_use <- cv_model[[lambda]]

  # Refit the model with the lambda chosen
  model <- glmnet::glmnet(
    x, y,
    alpha = 1, lambda = lambda_to_use, family = family
  )

  return(model)
}


# Evaluates a Lasso model's performance and extracts estimated coefficients.
#
# Args:
#   model: A fitted glmnet model object.
#   x: Predictor matrix (rows: observations, columns: features).
#   y: Response vector, same length as rows in 'x'.
#   threshold: Threshold for classification. Default is 0.5.
#
# Returns:
#   A list with:
#     - auc: Area under the ROC curve.
#     - acc: Classification accuracy.
#     - coefs: Estimated coefficients, excluding the intercept.

evaluate_lasso <- function(model, x, y, threshold = 0.5) {
  # Check whether x and y have the same # of obs
  if (nrow(x) != length(y)) {
    stop("The number of rows in 'x' does not match the length of 'y'.")
  }

  # Calculate predicted probabilities and predicted labels
  pred_probas <- predict(model, newx = x, type = "response")[, 1]
  pred_labels <- ifelse(pred_probas > threshold, 1, 0)

  # Evaluate by accuracy and AUC
  acc <- mean(pred_labels == y)
  auc <- pROC::auc(pROC::roc(y, pred_probas))

  # Extract estimated coefficients
  coefs <- as.vector(coef(model))[-1]

  results <- list(
    auc = auc,
    acc = acc,
    coefs = coefs
  )

  return(results)
}


# Executes a complete lasso analysis pipeline.
#
# Splits data into training and testing sets, performs cross-validation to
# select lambda, refits the model, and evaluates performance on the test set.
#
# Args:
#   x: Predictor matrix.
#   y: Response vector.
#   p: Proportion of data for training (0-1).
#   lambda: Lambda selection ("lambda.min" or "lambda.1se")
#   family: Default "binomial".
#   threshold: Classification threshold for binary outcome. Default 0.5.
#   seed: Seed for reproducibility.
#
# Returns:
#   A list with:
#     - auc: Area under the ROC curve.
#     - acc: Classification accuracy.
#     - coefs: Model coefficients.

run_lasso_pipeline <- function(
    x, y, p, lambda = c("lambda.min", "lambda.1se"), family = "binomial",
    threshold = 0.5, seed = NULL) {
  lambda <- match.arg(lambda)

  split <- train_test_split(x, y, p, seed)
  x_train <- split$train$x
  y_train <- split$train$y
  x_test <- split$test$x
  y_test <- split$test$y

  # Fit model with cross-validation
  model <- cv_lasso(x_train, y_train, lambda, family)

  # Evaluate model on the test split
  evaluation <- evaluate_lasso(model, x_test, y_test, threshold)

  # Add p-values calculated by hdi
  evaluation$pvals <- as.vector(hdi::lasso.proj(x, y)$pval)

  return(evaluation)
}
