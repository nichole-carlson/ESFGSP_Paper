# -------------------------------------------------------
# Binary Classification LASSO Pipeline
# -------------------------------------------------------

# Splits the dataset into training and testing sets.
.train_test_split <- function(x, y, p_train) {
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


# Computes the threshold that maximizes Youden's J statistic. Youden’s J =
# sensitivity + specificity - 1
.youden_cutoff <- function(labels, probs) {
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

  thresholds[which.max(youdens)]
}


# Calculate accuracy and AUC
.calculate_metrics <- function(true_labels, pred_probas, threshold) {
  pred_labels <- as.numeric(pred_probas > threshold)

  list(
    accuracy = mean(pred_labels == true_labels),
    auc = as.numeric(
      pROC::auc(pROC::roc(true_labels, pred_probas, quiet = TRUE))
    )
  )
}


# Each row of the returned matrix is the estimated coefs of a permutation.
.run_lasso_permutations <- function(x, y, lambda, n_perm) {
  p <- ncol(x)

  perm_coefs <- matrix(NA, nrow = n_perm, ncol = p)

  for (i in seq_len(n_perm)) {
    y_perm <- sample(y)
    model <- glmnet::glmnet(
      x, y_perm,
      alpha = 1, lambda = lambda, family = "binomial"
    )
    perm_coefs[i, ] <- as.vector(coef(model))[-1] # remove intercept
  }

  perm_coefs
}


# Calculate pvalue using projection based method
.get_hdi_pvals <- function(x, y, pval_adj = FALSE, ...) {
  if (!requireNamespace("hdi", quietly = TRUE)) {
    stop("Package 'hdi' is required")
  }

  result <- hdi::lasso.proj(x, y, family = "binomial", ...)

  if (pval_adj) {
    return(result$pval.corr)
  } else {
    return(result$pval)
  }
}


# --------------------
# Calculate two-sided p-values for LASSO coefficients using permutation test
compute_permutation_pvals <- function(perm_coefs, orig_coefs) {
  if (ncol(perm_coefs) != length(orig_coefs)) {
    stop("The dimensions of inputs do not match")
  }

  n_perm <- nrow(perm_coefs)

  perm_coefs <- abs(perm_coefs)
  orig_coefs <- abs(orig_coefs)

  (colSums(sweep(perm_coefs, 2, orig_coefs, FUN = ">=")) + 1) / (n_perm + 1)
}


# For binary outcome, use cross-validation to select optimal lambda,
# then evaluate on the test split
fit_evaluate_lasso <- function(
    x, y,
    lambda_choice = c("lambda.min", "lambda.1se"),
    p_train = 0.8,
    n_perm = 1000) {
  lambda_choice <- match.arg(lambda_choice)

  # Split data
  split <- .train_test_split(x, y, p_train)
  x_train <- split$train$x
  y_train <- split$train$y
  x_test <- split$test$x
  y_test <- split$test$y

  # Find optimal lambda with CV, then refit on the train split
  cv_model <- glmnet::cv.glmnet(
    x_train, y_train,
    alpha = 1, family = "binomial"
  )

  lambda <- cv_model[[lambda_choice]]

  model <- glmnet::glmnet(
    x_train, y_train,
    alpha = 1, lambda = lambda, family = "binomial"
  )
  coefs <- as.vector(coef(model))[-1]

  # Threshold from training
  train_probs <- predict(model, newx = x_train, type = "response")[, 1]
  threshold <- .youden_cutoff(y_train, train_probs)

  # Test performance
  test_probs <- predict(model, newx = x_test, type = "response")[, 1]
  metrics <- .calculate_metrics(y_test, test_probs, threshold)

  # Permutation coefs from training
  perm_coefs <- .run_lasso_permutations(x_train, y_train, lambda, n_perm)

  list(
    threshold = threshold,
    lambda = lambda,
    accuracy = metrics$accuracy,
    auc = metrics$auc,
    coefs = coefs,
    perm_coefs = perm_coefs
  )
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
# lasso_pixel_or_freq <- function(x, y, in_pixel_space = TRUE, e = NULL, ...) {
#   if (!in_pixel_space) {
#     if (is.null(e)) {
#       stop(
#         "Transformation matrix e is required when not in pixel space."
#       )
#     }
#     x <- x %*% e # map pixel to freq
#   }
#
#   fit_evaluate_lasso(x, y, ...)
# }
