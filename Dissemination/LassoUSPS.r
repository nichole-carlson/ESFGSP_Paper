library(tidyverse)
library(gridExtra)
library(fields)
library(MASS)
library(Matrix)
library(glmnet)

# Use the USPS handwritten digit data for 3, 5 and 8. It contains 1756 images,
# each image contains 256 pixels, with -1 corresponds to white and 1 to black.
data(USPSdigits, package = "IMIFA")
data58 <- USPSdigits$train[USPSdigits$train[, 1] %in% c(5, 8), ]
y <- data58[, 1]
x <- data58[, -1] + 1 # William's data used 0-2, and my data used -1 to 1


# Data preprocessing steps --------------------------------------------
x_adj <- x / 2 # scale from 0-2 to 0-1
# 2. nudge 0 and 1 1/10th the distance towards the nearest non-boundary
# values (globally)
min_pos <- min(x_adj[x_adj > 0])
max_below_one <- max(x_adj[x_adj < 1])
x_adj[x_adj == 0] <- min_pos / 10
x_adj[x_adj == 1] <- 1 - (1 - max_below_one) / 10
# 3. logit transformation
x_logit <- gtools::logit(x_adj)
# 4. standardize each picture
x_logit <- t(apply(x_logit, 1, function(x) (x - mean(x)) / sd(x)))


# Build the C matrix ------------------------------------------------
# Euclidean distance between two points in the matrix
dist <- function(v1, v2, ncol) {
  r1 <- (v1 - 1) %/% ncol
  c1 <- (v1 - 1) %% ncol
  r2 <- (v2 - 1) %/% ncol
  c2 <- (v2 - 1) %% ncol

  sqrt((r1 - r2)^2 + (c1 - c2)^2)
}

# For adjacenecy matrix
# distance <= dmax are considered adjacent
adj_lattice <- function(nrow, ncol, dmax = 1) {
  p <- nrow * ncol
  dist_matrix <- outer(1:p, 1:p, Vectorize(function(i, j) {
    dist(i, j, ncol)
  }))
  adj_mat <- dist_matrix <= dmax
  diag(adj_mat) <- 0

  adj_mat
}

c_adj <- adj_lattice(16, 16, dmax = 2)

# For Matern correlation matrix
# for each distance d, summarize the average semivariance
tg_vg_stats <- map_dfr(seq_len(nrow(x_logit)), function(i) {
  vg <- vgram.matrix(matrix(x_logit[i, ], 16, 16), R = 8)
  data.frame(d = vg$d.full, g = vg$vgram.full) %>%
    filter(!is.na(g)) %>%
    mutate(img = i)
}) %>%
  group_by(d) %>%
  summarise(
    mean_g = mean(g),
    lower_g = quantile(g, 0.025),
    upper_g = quantile(g, 0.975)
  )

# predict semivariance given phi, nu, d (distance), mean_g (actual semivariance)
calc_prop_vg <- function(phi, nu, d, g) {
  # Scale the distance by the range parameter
  x <- phi * d

  # Compute theoretical semivariance using the Matérn formula
  v <- 1 - besselK(x, nu) * x^nu / ((2^(nu - 1)) * gamma(nu))

  # Calculate means necessary for variance components estimation
  v_bar <- mean(v) # Mean of theoretical semivariance
  g_bar <- mean(g) # Mean of observed semivariance
  vg_bar <- mean(v * g)
  vsquare_bar <- mean(v^2) # Mean of the square of theoretical semivariance

  # Estimate partial sill (sigma2) and nugget effect (tau2)
  sigma2 <- (vg_bar - v_bar * g_bar) / (vsquare_bar - v_bar^2)
  tau2 <- g_bar - v_bar * sigma2

  # Adjustment for negative nugget effect
  tau2 <- if (tau2 < 0) {
    tau2 <- 0
    sigma2 <- vg_bar / vsquare_bar
  } else {
    tau2
  }

  # Return predicted semivariance
  tau2 + sigma2 * v
}

# search the optimal sum sqaure diff between obs and predicted semivariances
params_grid <- expand.grid(
  phi = seq(0.05, 5, 0.01),
  nu = seq(0.05, 10, 0.01)
)
ss_vec <- apply(params_grid, 1, function(p) {
  prop_vg <- calc_prop_vg(p[1], p[2], tg_vg_stats$d, tg_vg_stats$mean_g)
  sum((tg_vg_stats$mean_g - prop_vg)^2)
})
phi <- params_grid[which.min(ss_vec), 1]
nu <- params_grid[which.min(ss_vec), 2]

matern_corr <- function(nrow, ncol, phi, nu = 0.5, zerodiag = FALSE) {
  p <- nrow * ncol
  vec_corr <- Vectorize(function(i, j) {
    x <- dist(i, j, ncol) * phi # scaled distance
    besselK(x, nu) * (x^nu) / (2^(nu - 1) * gamma(nu))
  })
  corr_mat <- outer(1:p, 1:p, vec_corr)

  if (!zerodiag) {
    diag(corr_mat) <- 1
  } else {
    diag(corr_mat) <- 0
  }

  return(corr_mat)
}

c_mat <- matern_corr(16, 16, phi, nu, zerodiag = TRUE)

# For empirical matrix
c_emp <- cor(x_logit)
diag(c_emp) <- 0


# Eigen-decompose the C matrix ----------------------------------------
decompose_transform <- function(C, X) {
  p <- ncol(C)
  M <- diag(p) - matrix(1, p, p) / p
  eig_data <- eigen(M %*% C %*% M, symmetric = TRUE)

  order_idx <- order(eig_data$values, decreasing = TRUE)
  eig_vecs <- eig_data$vectors[, order_idx]
  eig_vals <- eig_data$values[order_idx]

  x_trans <- X %*% eig_vecs

  list(
    eigenvectors = eig_vecs,
    eigenvalues = eig_vals,
    x_trans = x_trans
  )
}

eig_results <- lapply(
  list(c_adj, c_mat, c_emp),
  function(c) decompose_transform(c, x_logit)
)
names(eig_results) <- c("adj", "mat", "emp")


# LASSO ---------------------------------------------------
# train/test split. Within test, perform cross-validation to find the optimal
# lambda, then evaluate the model performance on the test set using MSE.

perform_lasso <- function(x, y, test_frac = 0.2) {
  # Split data into training and test sets
  set.seed(123)
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * (1 - test_frac)))

  x_train <- as.matrix(x[train_idx, ])
  y_train <- y[train_idx]
  x_test <- as.matrix(x[-train_idx, ])
  y_test <- y[-train_idx]

  # Cross-validation to find the optimal lambda
  cv_model <- cv.glmnet(x_train, y_train, alpha = 1)
  lambda_best <- cv_model$lambda.min
  lambda_1se <- cv_model$lambda.1se

  # Fit Lasso models on training set with both lambda.min and lambda.1se
  model_min <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_best)
  model_1se <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_1se)

  # Evaluate model performance on the test set
  preds_min <- predict(model_min, s = lambda_best, newx = x_test)
  preds_1se <- predict(model_1se, s = lambda_1se, newx = x_test)

  # Convert predictions to binary outcomes
  preds_bin_min <- ifelse(preds_min > 0.5, 1, 0)
  preds_bin_1se <- ifelse(preds_1se > 0.5, 1, 0)

  # Calculate accuracy
  acc_min <- mean(preds_bin_min == y_test)
  acc_1se <- mean(preds_bin_1se == y_test)

  # Calculate AUC
  auc_min <- pROC::auc(pROC::roc(y_test, preds_min[, 1]))
  auc_1se <- pROC::auc(pROC::roc(y_test, preds_1se[, 1]))

  list(
    model_min = model_min, model_1se = model_1se,
    lambda_min = lambda_best, lambda_1se = lambda_1se,
    acc_min = acc_min, acc_1se = acc_1se,
    auc_min = auc_min, auc_1se = auc_1se
  )
}

y_ind <- ifelse(y == 5, 0, 1)

# fit using the full x (x orginal, x_logit, x_adj, x_mat, x_emp)
x_list <- c(list(x_adj, x_logit), map(eig_results, ~ .x$x_trans))
names(x_list) <- c("org", "logit", "adj", "mat", "emp")

full_fit <- map(x_list, ~ perform_lasso(.x, y_ind))
results <- map(full_fit, ~ c(
  acc_min = .x$acc_min, acc_1se = .x$acc_1se,
  auc_min = .x$auc_min, auc_1se = .x$auc_1se
))
results_df <- as.data.frame(do.call(rbind, results))
print(format(results_df, digits = 3))

#       acc_min acc_1se auc_min auc_1se
# org     0.945   0.950   0.996   0.995
# logit   0.959   0.959   0.998   0.996
# adj     0.964   0.973   0.996   0.995
# mat     0.964   0.964   0.995   0.995
# emp     0.959   0.968   0.997   0.998

# fit using x_adj, x_mat, x_emp with the top n eigenvectors
top_n_eig_vecs <- 10
top_n_fit <- map(
  eig_results, ~ perform_lasso(.x$x_trans[, seq_len(top_n_eig_vecs)], y_ind)
)
top_n_results <- map(top_n_fit, ~ c(
  acc_min = .x$acc_min, acc_1se = .x$acc_1se,
  auc_min = .x$auc_min, auc_1se = .x$auc_1se
))
top_n_results_df <- as.data.frame(do.call(rbind, top_n_results))
print(format(top_n_results_df, digits = 3))

#     acc_min acc_1se auc_min auc_1se
# adj   0.932   0.914   0.986   0.983
# mat   0.927   0.914   0.986   0.984
# emp   0.964   0.959   0.992   0.995
