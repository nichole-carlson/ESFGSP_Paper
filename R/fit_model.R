# -------------------------------------------------------
# Binary Classification LASSO Pipeline
# -------------------------------------------------------


# --------------------
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


# Transform coefs between spaces
.transform_coef <- function(coefs, e, to_freq = TRUE) {
  # Ensure coefs is a column matrix
  if (!is.matrix(coefs)) {
    coefs <- matrix(coefs, ncol = 1)
  }

  if (to_freq) {
    # pixel to freq: b = t(E) * beta
    t(e) %*% coefs
  } else {
    # freq to pixel: beta = E * b
    e %*% coefs
  }
}


# Calculate two-sided p-values for LASSO coefficients using permutation test
.compute_permutation_pvals <- function(perm_coefs, orig_coefs) {
  n_perm <- nrow(perm_coefs)

  perm_coefs <- abs(perm_coefs)
  orig_coefs <- abs(orig_coefs)

  (colSums(sweep(perm_coefs, 2, orig_coefs, FUN = ">=")) + 1) / (n_perm + 1)
}


# --------------------
# For binary outcome, use cross-validation to select optimal lambda,
# then evaluate on the test split.
#
# Return estiamted coefs in pixel and freq space, and corresponding pvals
fit_evaluate_lasso <- function(
    x, y,
    on_freq,
    transform_mat,
    lambda_choice = c("lambda.min", "lambda.1se"),
    p_train = 0.8,
    n_perm = 1000,
    ) {
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

  # Est coefs
  coef_orig <- matrix(as.vector(coef(model))[-1], ncol = 1)
  coef_trans <- .transform_coef(coef_orig, transform_mat, !on_freq)
  coef_arr <- cbind(coef_orig, coef_trans)
  colnames(coef_arr) <- c("orig", "trans")

  # Threshold from training
  train_probs <- predict(model, newx = x_train, type = "response")[, 1]
  threshold <- .youden_cutoff(y_train, train_probs)

  # Test performance
  test_probs <- predict(model, newx = x_test, type = "response")[, 1]
  metrics <- .calculate_metrics(y_test, test_probs, threshold)

  # Permutation coefs from training
  perm_coefs <- .run_lasso_permutations(x_train, y_train, lambda, n_perm)
  pval_orig <- .compute_permutation_pvals(perm_coefs, coef_orig)

  perm_trans <- t(.transform_coef(t(perm_coefs), transform_mat, !on_freq))
  pval_trans <- .compute_permutation_pvals(perm_trans, coef_trans)

  pval_arr <- cbind(pval_orig, pval_trans)
  colnames(pval_arr) <- c("orig", "trans")

  list(
    threshold = threshold,
    lambda = lambda,
    accuracy = metrics$accuracy,
    auc = metrics$auc,
    coef_arr = coef_arr,
    pval_arr = pval_arr,
    on_freq = on_freq
  )
}
