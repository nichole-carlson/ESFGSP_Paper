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

cv_lasso <- function(x, y, lambda = c("lambda.min", "lambda.1se"), family = "binomial") {
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
