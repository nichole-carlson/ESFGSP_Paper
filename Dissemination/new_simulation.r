library(MASS)
library(glmnet)
library(pROC)
library(tictoc)
library(table1)
library(ggplot2)
library(reshape2)
library(gridExtra)

current_dir <- getwd()
parent_dir <- dirname(current_dir)
fig_dir <- file.path(parent_dir, "Figures")
source("simulation.r")
source("simWrapper.r")


# ----- Simulation Plan -----
# - X, 1000 * 256, represents the covariate matrix in the pixel space; It was
# generated following multivariate normal distribution with an exponential
# correlation structure, denoted as W. Using eigen decomposition, we can get
# the eigenvectors V from W. Both V and W are 256 * 256.

# - Beta is the corresponding coefficient vector in the pixel space. It should
# be (256, ).

# - X_freq = X %*% V^T is the covariate matrix in the frequency space. It is
# also 1000 * 256.

# - b is the coefficient vector in the frequency space. beta = V^T %*% b.

# We do two simulations:
#   - the first one assumes sparsity in beta. When transforming beta back to
#     16 * 16, we make all areas except the central 8*8 to be zero. The values
#     in the central area will be decided later.
#   - the second one assumes sparsity in b. We randomly assign 10% values to be
#     non-zero.
# We will fit two models for each simulation: using X as covariates, and using
# X_freq as covariates.



# Generate an exponential correlation matrix for a square matrix with n_pixels
# Args:
#   n_pixels: Integer. The number of pixels (or elements) in the matrix.
# Returns:
#   A matrix where each element represents the correlation between two points,
#   based on the exponential decay of the Euclidean distance between them.
gen_exp_corr_matrix <- function(n_pixels) {
  num_col <- sqrt(n_pixels)

  # Function to calculate the Euclidean distance between two points in a 2D grid
  # Args:
  #   i: Integer. Index of the first point in 1D.
  #   j: Integer. Index of the second point in 1D.
  # Returns:
  #   Numeric. The Euclidean distance between the two points.
  calc_dist <- function(i, j) {
    row_i <- (i - 1) %/% num_col
    col_i <- (i - 1) %% num_col
    row_j <- (j - 1) %/% num_col
    col_j <- (j - 1) %% num_col
    sqrt((row_i - row_j)^2 + (col_i - col_j)^2)
  }

  # Generate the correlation matrix using the vectorized distance function
  outer(
    1:n_pixels, 1:n_pixels,
    Vectorize(function(i, j) exp(-calc_dist(i, j)))
  )
}

# Generate a sparse coefficient vector
# Args:
#   len: Integer. The length of the coefficient vector.
#   sparsity: Numeric. Proportion of non-zero elements in the vector. Default is 0.1.
#   effect: Numeric. The effect size assigned to non-zero elements. Default is 0.1.
#   seed: Integer. Random seed for reproducibility.
# Returns:
#   A numeric vector of specified length with a specified proportion of non-zero elements.
create_sparse_vec <- function(len, sparsity = 0.1, effect = 0.1, seed) {
  set.seed(seed)
  vec <- rep(0, len)
  non_zero_idx <- sample(1:len, size = floor(len * sparsity), replace = FALSE)
  vec[non_zero_idx] <- effect
  vec
}

# Define beta for simulation with the center as the black region
# Args:
#   img_size: Integer. The dimension of the square image (number of rows/columns).
#   beta_effect: Numeric. The effect size assigned to the center region.
# Returns:
#   A numeric vector of length img_size^2 with the center region set to beta_effect.
define_center_beta <- function(img_size, beta_effect) {
  p <- img_size^2
  beta <- rep(0, p)
  ctr_start <- (img_size - 8) / 2 + 1
  ctr_end <- (img_size + 8) / 2
  ctr_idx <- as.vector(matrix(1:p, nrow = img_size)[
    ctr_start:ctr_end,
    ctr_start:ctr_end
  ])
  beta[ctr_idx] <- beta_effect
  return(beta)
}


# Generate the X matrix
# Args:
#   n_samples: Integer. The number of samples (rows) to generate.
#   img_size: Integer. The dimension of the square image (number of rows/columns).
#   cov_matrix: Matrix. The covariance matrix used for generating
#   multivariate normal samples.
# Returns:
#   A matrix of generated samples, with n_samples rows and img_size^2 columns.
gen_X_matrix <- function(n_samples, img_size, cov_matrix) {
  mvrnorm(n_samples, mu = rep(0, img_size^2), Sigma = cov_matrix)
}

# Generate probabilities
# Args:
#   x_matrix: Matrix. The design matrix of predictor variables.
#   coeffs: Numeric vector. The coefficient vector.
# Returns:
#   A numeric vector of probabilities for each sample.
calc_probs <- function(x_matrix, coeffs) {
  eta <- x_matrix %*% coeffs
  1 / (1 + exp(-eta))
}

# Generate response variables
# Args:
#   x_matrix: Matrix. The design matrix of predictor variables.
#   coeffs: Numeric vector. The coefficient vector.
# Returns:
#   A numeric vector of binary response variables generated using the probabilities.
gen_responses <- function(x_matrix, coeffs) {
  p <- calc_probs(x_matrix, coeffs)
  rbinom(nrow(x_matrix), 1, p)
}

# Function for Simulation 1
simulate_1 <- function(i, size, num_samples, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  beta <- define_beta(size, effect)
  x <- generate_X(num_samples, size, w)
  y <- generate_response(x, beta)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}

# Function to compare beta effects for Simulation 1
compare_beta_effects_sim1 <- function(effects, size = 16, num_samples = 1000, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))
  set.seed(seed)
  w <- generate_cov_matrix(size)
  for (effect in effects) {
    beta <- define_beta(size, effect)
    x <- generate_X(num_samples, size, w)
    p <- generate_probs(x, beta)
    hist(p,
      breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
      xlab = "Probability p", col = "lightblue", border = "black"
    )
  }
}

png(
  file = file.path(fig_dir, "sim1_p_dist.png"),
  width = 1600, height = 1200, res = 150
)
compare_beta_effects_sim1(c(1, 0.2, 0.1, 0.05, 0.01))
dev.off()

# Run Simulation 1
# For each iteration:
#   - Split into train and test dataset, calculate AUC and accuracy, estimated
#     coefficients under lambda.min and lambda.1se
#   - Perform permutation test. Estimate p-values for pixels

tic()
sim1_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_1(i, size = 16, num_samples = 1000, effect = 0.1, p_train = 0.8, n_perm = 100, seed = 42),
  list_export = c(
    "generate_cov_matrix", "define_beta", "generate_X", "generate_probs",
    "generate_response", "simulate_1", "perform_lasso", "perm_lasso"
  ),
  list_package = c("MASS", "glmnet", "pROC", "foreach", "doParallel")
)
toc()

# # Save Simulation 1 results
# save(
#   sim1_output,
#   file = file.path(parent_dir, "Simulations", paste0("sim1_", format(Sys.time(), "%y%m%d"), ".RData"))
# )
load(file = file.path(current_dir, "Simulations", "sim1_240701.RData"))

# # Generate summary stats for ACC and AUCs
acc_auc_sim1 <- as.data.frame(sim1_output[, 1:4])
colnames(acc_auc_sim1) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = acc_auc_sim1)

# # Calculate significance percentages for Simulation 1
pvals_matrix_sim1 <- sim1_output[, -(1:4)]
sig_perc_sim1 <- calc_sig_perc(pvals_matrix_sim1, method = "bonferroni", alpha = 0.05)

# Create plots for Simulation 1
beta_sim1 <- define_beta(16, 0.1)
plot_beta_sim1 <- plot_heatmap(beta_sim1)
plot_orig_sim1 <- plot_heatmap(sig_perc_sim1$perc_orig, c(0, 100))
plot_adj_sim1 <- plot_heatmap(sig_perc_sim1$perc_adj, c(0, 100))

# Arrange plots for formal presentation for Simulation 1
sim1_heatmap <- grid.arrange(
  plot_beta_sim1, plot_orig_sim1, plot_adj_sim1,
  ncol = 2
)
ggsave(
  filename = file.path(current_dir, "Figures", "sim1_heatmap.png"),
  plot = sim1_heatmap, width = 10, height = 8
)

# Function for Simulation 2
simulate_2 <- function(i, size, num_samples, sparsity, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  v <- eigen_decomp(w)$vectors
  x <- generate_X(num_samples, size, diag(size^2))
  b <- generate_sparse_vector(size^2, sparsity, effect, seed)
  y <- generate_response(x, b)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}

# # Function to compare beta effects for Simulation 2
compare_beta_effects_sim2 <- function(effects, size = 16, num_samples = 1000, sparsity = 0.1, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))
  set.seed(seed)
  w <- generate_cov_matrix(size)
  V <- eigen_decomp(w)$vectors
  for (effect in effects) {
    b <- generate_sparse_vector(size^2, sparsity, effect, seed)
    x <- generate_X(num_samples, size, diag(size^2))
    p <- generate_probs(x, b)
    hist(p, breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1), xlab = "Probability p", col = "lightblue", border = "black")
  }
}

png(
  file = file.path(fig_dir, "sim2_p_dist.png"),
  width = 1600, height = 1200, res = 150
)
compare_beta_effects_sim2(c(1, 0.8, 0.6, 0.4, 0.2, 0.1))
dev.off()

# Run Simulation 2
tic()
sim2_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_2(i, size = 16, num_samples = 1000, sparsity = 0.1, effect = 0.4, p_train = 0.8, n_perm = 100, seed = 42),
  list_export = c(
    "generate_cov_matrix", "eigen_decomp", "generate_sparse_vector", "generate_X", "generate_probs", "generate_response", "perform_lasso",
    "perm_lasso", "simulate_2"
  ),
  list_package = c("MASS", "glmnet", "pROC", "foreach", "doParallel")
)
toc()

# Save Simulation 2 results
# save(
#   sim2_output,
#   file = file.path(parent_dir, "Simulations", paste0("sim2_", format(Sys.time(), "%y%m%d"), ".RData"))
# )
load(file = file.path(current_dir, "Simulations", "sim2_240701.RData"))

# Generate summary stats for ACC and AUCs for Simulation 2
acc_auc_sim2 <- as.data.frame(sim2_output[, 1:4])
colnames(acc_auc_sim2) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = acc_auc_sim2)

# Calculate significance percentages for Simulation 2
pvals_matrix_sim2 <- sim2_output[, -(1:4)]
sig_perc_sim2 <- calc_sig_perc(pvals_matrix_sim2, method = "bonferroni", alpha = 0.05)


# Create plots for Simulation 2
b_sim2 <- generate_sparse_vector(16^2, 0.1, 0.4, 42)
plot_beta_sim2 <- plot_heatmap(b_sim2)
plot_orig_sim2 <- plot_heatmap(sig_perc_sim2$perc_orig, c(0, 100))
plot_adj_sim2 <- plot_heatmap(sig_perc_sim2$perc_adj, c(0, 100))

# Arrange plots for formal presentation for Simulation 2
sim2_heatmap <- grid.arrange(
  plot_beta_sim2, plot_orig_sim2, plot_adj_sim2,
  ncol = 2
)
ggsave(
  filename = file.path(current_dir, "Figures", "sim2_heatmap.png"),
  plot = sim2_heatmap, width = 10, height = 8
)
