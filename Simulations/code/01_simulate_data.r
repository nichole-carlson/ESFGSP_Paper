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

# Matrix notation
# If X is a 256*1 matrix, V is the matrix of eigenvectors, with each column
# representing an eigenvector, then we know from the content above that
# X_freq = t(V) %*% X, which is also a 256*1 matrix. Suppose know we generate
# 1000 samples of X, and let X_mat to be a 1000*256 matrix, with each row
# as a X, then X_freq_mat = t( t(V) %*% t(X_mat) ) = X_mat %*% V.

# If beta and b are both 256*1 matrix, then t(X) %*% beta = t(X_freq) %*% b.
# t(X_freq) %*% b = t(X) %*% V %*% b. So beta = V %*% b


# Load required libraries
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")



# ----- Functions to simulate data -----

# Generate an exponential correlation matrix: suppose the pixels can
# construct a squared image. Calculate the correlation between all pairs
# of pixels based on their location in the image.
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

# Generate a sparse coefficient vector in the freq space
# Args:
#   len: Integer. Length of the coefficient vector.
#   sparsity: Numeric. Proportion of non-zero elements.
#   effect_size: Numeric. Effect size of non-zero elements.
#   seed: Integer. Random seed for reproducibility.
# Returns:
#   A numeric vector with specified sparsity and effect size.
gen_b <- function(len, sparsity, effect_size) {
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

# Transform between b and beta
b_to_beta <- function(b, V) {
  V %*% b
}

beta_to_b <- function(beta, V) {
  t(V) %*% beta
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

gen_meta_data <- function(img_size) {
  n_pixels <- img_size^2

  W <- gen_exp_corr(n_pixels) # correlation matrix for x
  V <- eigen(W)$vectors
  W_freq <- t(V) %*% W %*% V # correlation matrix for x_freq

  list(x_cov = W, eig_vecs = V, x_freq_cov = W_freq)
}

simulate_xy_freq <- function(n_iter, n_samples, img_size, freq_coefs) {
  meta_data <- gen_meta_data(img_size)
  x_freq_cov <- meta_data$x_freq_cov
  eig_vecs <- meta_data$eig_vecs

  results <- list()

  for (i in seq_len(n_iter)) {
    if (i %% 10 == 0) {
      cat("Simulating", i, "of", n_iter, "iterations \n")
    }
    x_freq <- gen_x(n_samples, x_freq_cov)
    x <- x_freq %*% t(eig_vecs)
    y <- gen_y(x_freq, freq_coefs)
    results[[i]] <- list(x = x, x_freq = x_freq, y = y)
  }

  results
}

run_simulation1 <- function(n_iter, n_samples, img_size, beta_effect_size, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # calculate meta data, W, V, W_freq
  meta_data <- gen_meta_data(img_size)
  x_freq_cov <- meta_data$x_freq_cov
  eig_vecs <- meta_data$eig_vecs

  # generate coefs on pixel space
  beta <- gen_beta(img_size, beta_effect_size)
  b <- beta_to_b(beta, eig_vecs)

  # run simulation on freq space
  simulated_data <- simulate_xy_freq(n_iter, n_samples, img_size, b)

  results <- list(
    meta_data = c(meta_data, list(beta = beta, b = b)),
    data = simulated_data
  )

  return(results)
}

run_simulation2 <- function(n_iter, n_samples, img_size, b_effect_size, sparsity, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # calculate meta data, W, V, W_freq
  meta_data <- gen_meta_data(img_size)
  x_freq_cov <- meta_data$x_freq_cov
  eig_vecs <- meta_data$eig_vecs

  # generate coefs on pixel space
  b <- gen_b(img_size^2, sparsity, b_effect_size)
  beta <- b_to_beta(b, eig_vecs)

  # run simulation on freq space
  simulated_data <- simulate_xy_freq(n_iter, n_samples, img_size, b)

  results <- list(
    meta_data = c(meta_data, list(beta = beta, b = b)),
    data = simulated_data
  )

  return(results)
}


# ----- Choose beta effect size and b effect size -----

# Function to compare beta effects for Simulation 1
compare_beta_effects <- function(effects, n_samples = 1000, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))
  set.seed(seed)

  meta <- gen_meta_data(img_size = 16)
  x_cov <- meta$x_cov
  x <- gen_x(n_samples, x_cov)

  for (effect in effects) {
    beta <- gen_beta(img_size = 16, effect_size = effect)
    p <- 1 / (1 + exp(-(x %*% beta)))
    hist(p,
      breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
      xlab = "Probability p", col = "lightblue", border = "black"
    )
  }
}

png(
   file = file.path(fig_dir, "sim1_p_dist.png"),
   width = 1200, height = 800, res = 150
)
compare_beta_effects(c(1, 0.2, 0.1, 0.05, 0.01))
dev.off()

# Function to compare b effects for Simulation 2
compare_b_effects <- function(effects, n_samples = 1000, sparsity = 0.1, seed = 42) {
  set.seed(seed)
  par(mfrow = c(ceiling(length(effects) / 2), 2))

  meta <- gen_meta_data(img_size = 16)
  x_freq_cov <- meta$x_freq_cov
  x_freq <- gen_x(n_samples, x_freq_cov)

  for (effect in effects) {
    b <- gen_b(len = 256, sparsity, effect)
    p <- 1 / (1 + exp(-(x_freq %*% b)))
    hist(p, breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1), xlab = "Probability p", col = "lightblue", border = "black")
  }
}

png(
   file = file.path(fig_dir, "sim2_p_dist.png"),
   width = 1200, height = 800, res = 150
)
compare_b_effects(c(1, 0.8, 0.6, 0.4, 0.2, 0.1))
dev.off()



# ----- Run Simulations -----

n_sim <- 500
n_samples <- 1000
img_size <- 16
beta_effect <- 0.1
b_effect <- 0.4
b_sparsity <- 0.1
seed <- 42

cat("Simulating Data for", n_sim, "iterations of Simulation 1 ... \n")
tic()
sim1_data <- run_simulation1(
  n_sim, n_samples, img_size, beta_effect, seed
)
toc()

cat("Simulating Data for", n_sim, "iterations of Simulation 2 ... \n")
tic()
sim2_data <- run_simulation2(
  n_sim, n_samples, img_size, b_effect, b_sparsity, seed
)
toc()

filename <- paste0("simulated_data_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  n_sim, n_samples, img_size, beta_effect, b_effect, b_sparsity, seed,
  sim1_data, sim2_data,
  file = file.path(results_data_dir, filename)
)
