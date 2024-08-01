# ----- Simulation Plan -----

# Two simulations: sparsity in pixel and frequency spaces.

# Setup
# X: a column vector representing 256 pixels in the pixel space. Its covariance
# matrix, Sigma, follows an exponential correlation structure. Suppose i and j
# are two pixels in X, Sigma_ij = -exp(dist(i, j)), where dist(i, j) are
# calculated based on their location in the 16*16 matrix.

# Let V be the matrix of eigenvectors of Sigma, with each column representing
# an eigenvector. X_freq = t(V) %*% X, representing X in the frequency space.
# Its covariance matrix is t(V) %*% Sigma %*% V, which is a diagonal matrix.

# Simulations
#   In each iteration, randomly generate X_freq, a vector, from a multivariate
# normal distribution with the covariance matrix equals to t(V) %*% Sigma %*% V.
# We repeat this process 1000 times. Then we calculate X = V %*% X_freq.
#   In simulation 1, we assume sparsity in beta, the coefficient vector for X.
# When converting the (256, ) vector into a 16*16 matrix, it should only have
# non-zero values in the central 8*8 region.
#   In simulation 2, we assume sparsity in b, the coefficient vector for X_freq.
# Most of its 256 entries will be zero, with a randomly 10% be non-zero.

# Model Fitting
# Fit LASSO models using pixel and frequency space covariates.
# Split data: 80% training, 20% test.
# Tune lambda via 10-fold cross-validation.
# Optimal lambda: lowest average binomial deviance.

# Evaluation
# Metrics: accuracy and AUC.
# Permutation test (100 iterations) for p-values.
# Calculate mean, std deviation, and significant p-values percentage.

# Matrix notation
# If X is a 256*1 matrix, V is the matrix of eigenvectors, with each column
# representing an eigenvector, then we know from the content above that
# X_freq = t(V) %*% X, which is also a 256*1 matrix. Suppose know we generate
# 1000 samples of X, and let X_mat to be a 1000*256 matrix, with each row
# as a X, then X_freq_mat = t( t(V) %*% t(X_mat) ) = X_mat %*% V.

# Load required libraries
library(MASS)
library(glmnet)
library(pROC)
library(tictoc)
library(table1)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Set directories
parent_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"
fig_dir <- file.path(parent_dir, "Figures")
source("simulation.r")
source("simWrapper.r")



# ----- Define Functions -----

# Generate an exponential correlation matrix
# Suppose the pixels can construct a squared image. Calculate the correlation
# between all pairs of pixels based on their location in the image.

# Args:
#   n_pixels: Integer. The number of pixels (or elements) in the matrix.
# Returns:
#   A matrix where each element represents the correlation between two points,
#   based on the exponential decay of the Euclidean distance between them.
gen_exp_corr <- function(n_pixels) {
  n_cols <- sqrt(n_pixels)

  # Calculate Euclidean distance between two points in a 2D grid
  # Args:
  #   i: Integer. Index of the first point in 1D.
  #   j: Integer. Index of the second point in 1D.
  # Returns:
  #   Numeric. The Euclidean distance between the two points.
  calc_dist <- function(i, j) {
    row_i <- (i - 1) %/% n_cols
    col_i <- (i - 1) %% n_cols
    row_j <- (j - 1) %/% n_cols
    col_j <- (j - 1) %% n_cols
    sqrt((row_i - row_j)^2 + (col_i - col_j)^2)
  }

  # Generate the correlation matrix using the vectorized distance function
  outer(1:n_pixels, 1:n_pixels, Vectorize(function(i, j) exp(-calc_dist(i, j))))
}

# Eigen decomposition of a covariance matrix. Decreasingly order by eigenvalues
eigen_decomp <- function(mat) {
  eig_res <- eigen(mat, symmetric = TRUE)

  ord_idx <- order(eig_res$values, decreasing = TRUE)
  eig_vecs <- eig_res$vectors[, ord_idx]
  eig_vals <- eig_res$values[ord_idx]

  return(list(vectors = eig_vecs, values = eig_vals))
}


# Generate a sparse coefficient vector in the freq space
# Args:
#   len: Integer. Length of the coefficient vector.
#   sparsity: Numeric. Proportion of non-zero elements.
#   effect_size: Numeric. Effect size of non-zero elements.
#   seed: Integer. Random seed for reproducibility.
# Returns:
#   A numeric vector with specified sparsity and effect size.

gen_b <- function(len, sparsity, effect_size, seed) {
  set.seed(seed)
  vec <- rep(0, len)
  nz_idx <- sample(seq_len(len), size = floor(len * sparsity), replace = FALSE)
  vec[nz_idx] <- effect_size
  return(vec)
}

# Generate a sparse coefficient vector in the pixel space
# Args:
#   img_size: Integer. Dimension of the square image (rows/columns).
#   effect_size: Numeric. Effect size assigned to the center region.
# Returns:
#   A numeric vector with the center region set to the effect size.

gen_beta <- function(img_size, effect_size) {
  len <- img_size^2
  beta <- rep(0, len)
  ctr_start <- (img_size - 8) / 2 + 1
  ctr_end <- (img_size + 8) / 2
  ctr_idx <- as.vector(matrix(1:len, nrow = img_size)[
    ctr_start:ctr_end,
    ctr_start:ctr_end
  ])
  beta[ctr_idx] <- effect_size
  return(beta)
}

# Generate a matrix of samples from a multivariate normal distribution
# Args:
#   n_samples: Integer. Number of samples (rows).
#   cov_matrix: Matrix. Covariance matrix.
# Returns:
#   A matrix of samples with n_samples rows and img_size^2 columns.
gen_x <- function(n_samples, cov_matrix) {
  n_pixels <- nrow(cov_matrix)
  samples <- mvrnorm(n_samples, mu = rep(0, n_pixels), Sigma = cov_matrix)
  if (n_samples == 1) {
    samples <- matrix(samples, nrow = 1)
  }
  return(samples)
}

# Generate binary response variables
# Args:
#   x_mat: Matrix. Design matrix of predictor variables.
#   coefs: Numeric vector. Coefficient vector.
# Returns:
#   A numeric vector of binary response variables.
gen_y <- function(x_mat, coefs) {
  probs <- 1 / (1 + exp(-(x_mat %*% coefs)))
  rbinom(nrow(x_mat), 1, probs)
}



# ----- Choose beta effect size and b effect size -----

# Function to compare beta effects for Simulation 1
compare_beta_effects <- function(effects, n_samples = 1000, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))
  set.seed(seed)

  W <- gen_exp_corr(256)
  V <- eigen_decomp(W)$vectors
  W_freq <- t(V) %*% W %*% V

  for (effect in effects) {
    beta <- gen_beta(img_size = 16, effect_size = effect)
    x_freq <- gen_x(n_samples, W_freq)
    x <- x_freq %*% t(V)
    p <- 1 / (1 + exp(-(x %*% beta)))
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
compare_beta_effects(c(1, 0.2, 0.1, 0.05, 0.01))
dev.off()

# Function to compare b effects for Simulation 2
compare_b_effects <- function(effects, n_samples = 1000, sparsity = 0.1, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))

  W <- gen_exp_corr(256)
  V <- eigen_decomp(W)$vectors
  W_freq <- t(V) %*% W %*% V

  for (effect in effects) {
    b <- gen_b(len = 256, sparsity, effect, seed)
    x_freq <- gen_x(n_samples, W_freq)
    p <- 1 / (1 + exp(-(x_freq %*% b)))
    hist(p, breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1), xlab = "Probability p", col = "lightblue", border = "black")
  }
}

png(
  file = file.path(fig_dir, "sim2_p_dist.png"),
  width = 1600, height = 1200, res = 150
)
compare_b_effects(c(1, 0.8, 0.6, 0.4, 0.2, 0.1))
dev.off()

beta_effect_size <- 0.1
b_effect_size <- 0.4


# ----- Simulation Functions -----



# Perform a single simulation run
# Args:
#   img_size: Integer. Dimension of the square image (rows/columns).
#   n_samples: Integer. Number of samples (rows).
#   seed: Integer. Random seed for reproducibility.
# Returns:
#   A list containing the design matrix, frequency space design matrix,
#   and responses.
run_single_sim <- function(img_size, n_samples, beta_vec, seed) {
  x <- gen_x(n_samples, img_size, cov_matrix)
  x_freq <- x %*% eig_vecs
  responses <- gen_responses(x, beta)
  return(list(x = x, x_freq = x_freq, responses = responses))
}







# Run multiple simulations
# Args:
#   n_sim: Integer. Number of simulations to run.
#   img_size: Integer. Dimension of the square image (rows/columns).
#   n_samples: Integer. Number of samples (rows) for each simulation.
#   beta: Numeric vector. Coefficient vector in pixel space or frequency space.
#   seed: Integer. Seed for reproducibility.
# Returns:
#   A list containing the results of each simulation.
run_all_sims <- function(n_sim, img_size, n_samples, cov_matrix, eig_vecs, beta, seed = 123) {
  set.seed(seed)
  results <- vector("list", n_sim)
  for (i in 1:n_sim) {
    results[[i]] <- run_single_sim(img_size, n_samples, cov_matrix, eig_vecs, beta)
  }
  return(results)
}

# Fit models
# Args:
#   data: List. Contains x_mat, x_freq_mat, and responses.
# Returns:
#   A list containing fitted models using x_mat and x_freq_mat.
fit_models <- function(data) {
  x_mat <- data$x_mat
  x_freq_mat <- data$x_freq_mat
  y <- data$responses

  model_x <- glmnet(x_mat, y, family = "binomial")
  model_x_freq <- glmnet(x_freq_mat, y, family = "binomial")

  return(list(model_x = model_x, model_x_freq = model_x_freq))
}

# Run Simulation 1
run_simulation_1 <- function(img_size, n_samples, beta_effect, n_sim, seed) {
  cov_matrix <- gen_exp_corr(img_size^2)
  eig_vecs <- eigen(cov_matrix)$vectors
  beta <- define_center_beta(img_size, beta_effect)
  b <- t(eig_vecs) %*% beta

  sim_results <- run_all_sims(n_sim, img_size, n_samples, cov_matrix, eig_vecs, beta, seed)

  models <- lapply(sim_results, fit_models)

  return(models)
}


# Run Simulation 2
run_simulation_2 <- function(img_size, n_samples, sparsity, effect_size, n_sim, seed) {
  cov_matrix <- diag(img_size^2)
  eig_vecs <- eigen(cov_matrix)$vectors
  b <- gen_sparse_vec(img_size^2, sparsity, effect_size, seed)

  sim_results <- run_all_sims(n_sim, img_size, n_samples, cov_matrix, eig_vecs, b, seed)

  models <- lapply(sim_results, fit_models)

  return(models)
}

# ----- Run Simulations -----

# Parameters
img_size <- 16
n_samples <- 1000
beta_effect <- 0.1
sparsity <- 0.1
effect_size <- 0.1
n_sim <- 2
seed <- 123

# Run Simulation 1
models_sim1 <- run_simulation_1(img_size, n_samples, beta_effect, n_sim, seed)

# Run Simulation 2
models_sim2 <- run_simulation_2(img_size, n_samples, sparsity, effect_size, n_sim, seed)



# ----- Simulation 1 -----

w <- generate_exp_corr_matrix(256)
v <- perform_eigen_decomp(w)$vectors
beta <- define_center_beta(16, 0.1)
b <- v %*% beta

x <- gen_X_matrix(1000, 16, w)
x_freq <- x %*% t(v)
y <- generate_response(x, beta)

# Function for Simulation 1
simulate_1 <- function(i, size, n_samples, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  beta <- define_beta(size, effect)
  x <- generate_X(n_samples, size, w)
  y <- generate_response(x, beta)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}




# Run Simulation 1
# For each iteration:
#   - Split into train and test dataset, calculate AUC and accuracy, estimated
#     coefficients under lambda.min and lambda.1se
#   - Perform permutation test. Estimate p-values for pixels

tic()
sim1_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_1(i, size = 16, n_samples = 1000, effect = 0.1, p_train = 0.8, n_perm = 100, seed = 42),
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
simulate_2 <- function(i, size, n_samples, sparsity, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  v <- eigen_decomp(w)$vectors
  x <- generate_X(n_samples, size, diag(size^2))
  b <- generate_sparse_vector(size^2, sparsity, effect, seed)
  y <- generate_response(x, b)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}



# Run Simulation 2
tic()
sim2_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_2(i, size = 16, n_samples = 1000, sparsity = 0.1, effect = 0.4, p_train = 0.8, n_perm = 100, seed = 42),
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
