# ----- Simulation Plan -----

# Two simulations: sparsity in pixel and frequency spaces.

# Setup
# Suppose x is a column vector representing 256 pixels in the pixel space. It's
# covariance matrix, Sigma, follows an exponential correlation structure:
# suppose i and j are two pixels in x, sigma_ij = exp(-dist(i, j)), where
# dist(i, j) is calculated based on their 2-D location in the 16*16 matrix.
# 
# Let V be the matrix of eigenvectors of Sigma, calculated using the method
# brought in the Spatial filtering paper, with each column representing an
# eigenvector. x_freq = t(V) %*% x, representing x in the freq space. It's
# covariance matrix is t(V) %*% Sigma %*% V, which should be a diagonal matrix.
# 
# Matrix notation
# Suppose X reperesents a n*256 matrix, with each row as t(x) in the above
# context, then we can calculate X_freq = X %*% V, which is also a n*256 matrix. 
# 
# We further assume beta and b as 256*1, representing the coefficient vector
# in the pixel space and freq space, respectively. Since X_freq %*% b = X %*%
# beta, we can get V %*% b = beta.
# 
# Simulations
# We do two simulations. The first one assume sparsity in beta. When converting
# beta into a 16*16 matrix, it should only have non-zero values in the central
# 8*8 region. We will simulate the covariates and outcome in the pixel space. The second
# one assumes sparsity in b. We will simulate the covariates and outcome in the
# freq space this time.


# Load required libraries
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")



# ----- Functions to simulate data -----

eigen_decomp <- function(mat) {
  n_cols <- ncol(mat)
  cent_mat <- diag(n_cols) - matrix(1, n_cols, n_cols) / n_cols
  eig_res <- eigen(cent_mat %*% mat %*% cent_mat, symmetric = TRUE)

  ord_idx <- order(eig_res$values, decreasing = TRUE)
  eig_vecs <- eig_res$vectors[, ord_idx]
  eig_vals <- eig_res$values[ord_idx]

  return(list(vectors = eig_vecs, values = eig_vals))
}

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

# Generates a diagonal matrix with values decay faster at the beginning and 
# slow down towards the end.
# d_i = start_value * exp(-k * log(1 + b * i))
# k is the decay rate
# b is a constant that controls the rate of initial decay
# gen_diag_corr <- function(n, start_value = 6, end_value = 0.5, b = 0.1 ) {
#   k <- log(start_value / end_value) / log(1 + b * (n - 1))
#   diag_vals <- start_value * exp(-k * log(1 + b * seq(0, n - 1)))
#   diag_mat <- diag(diag_vals)
# 
#   return(diag_mat)
# }


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

run_sim1 <- function(n_iter, n_samples, img_size, beta_effect_size, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # generate correlation structure in the pixel space
  n_pixels <- img_size^2
  x_cov <- gen_exp_corr(n_pixels)
  eig_res <- eigen_decomp(x_cov)
  eig_vecs <- eig_res$vectors
  eig_vals <- eig_res$values

  # generate coefs in the pixel space
  beta <- gen_beta(img_size, beta_effect_size)
  b <- beta_to_b(beta, eig_vecs)

  # generate covariates in the pixel space (X)
  x <- gen_x(n_samples * n_iter, x_cov)
  x_freq <- x %*% eig_vecs
  y <- gen_y(x, beta)

  # return
  res <- list(
    meta_data = list(
      n_iter = n_iter, n_samples = n_samples, 
      x_cov = x_cov, eig_vecs = eig_vecs, eig_vals = eig_vals,
      beta = beta, b = b
    ),
    data = list(x = x, y = y, x_freq = x_freq)
  )
  return(res)
}

run_sim2 <- function(n_iter, n_samples, img_size, b_effect_size, sparsity, seed = NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
 
  # generate correlation structure in the freq space
  n_pixels <- img_size^2
  x_freq_cov <- diag(seq(6, 0.5, length.out = n_pixels))
  x_cov <- gen_exp_corr(n_pixels) # same assumption of exp corr structre for X
  eig_res <- eigen_decomp(x_cov)
  eig_vecs <- eig_res$vectors
  eig_vals <- eig_res$values

  # generate coefs in the freq space
  b <- gen_b(n_pixels, sparsity, b_effect_size)
  beta <- b_to_beta(b, eig_vecs)

  # generate covariates in the freq space (X_freq)
  x_freq <- gen_x(n_samples * n_iter, x_freq_cov)
  x <- x_freq %*% t(eig_vecs)
  y <- gen_y(x_freq, b)

  # return
  res <- list(
    meta_data = list(
      n_iter = n_iter, n_samples = n_samples,
      x_freq_cov = x_freq_cov, eig_vecs = eig_vecs, eig_vals = eig_vals,
      beta = beta, b = b
    ),
    data = list(x = x, y = y, x_freq = x_freq)
  )
  return(res)
}


# ----- Choose beta effect size and b effect size -----

# Choose beta effect size for sim1
png(
   file = file.path(fig_dir, "sim1_p_dist.png"),
   width = 1200, height = 800, res = 150
)

beta_effects <- c(1, 0.2, 0.1, 0.05, 0.01)
par(mfrow = c(ceiling(length(beta_effects) / 2), 2))
for (effect in beta_effects){
  sim1_1iter <- run_sim1(
    n_iter = 1, n_samples = 1000, img_size = 16, 
    beta_effect_size = effect, seed = 42
  )  
  x <- sim1_1iter$data$x 
  beta <- sim1_1iter$meta_data$beta
  p <- 1 / (1 + exp(-(x %*% beta)))
  hist(
    p,
    breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
    xlab = "Probability p", col = "lightblue", border = "black"
  )
}

dev.off()

# Choose b effect for sim2
png(
 file = file.path(fig_dir, "sim2_p_dist.png"),
 width = 1200, height = 800, res = 150
)

b_effects <- seq(0.3, 0.1, by = -0.05)
par(mfrow = c(ceiling(length(b_effects) / 2), 2))
for (effect in b_effects){
  sim2_1iter <- run_sim2(
    n_iter = 1, n_samples = 1000, img_size = 16,
    b_effect_size = effect, sparsity = 0.1, seed = 42
  )
  x_freq <- sim2_1iter$data$x_freq
  b <- sim2_1iter$meta_data$b
  p <- 1 / (1 + exp(-(x_freq %*% b)))
  hist(
    p,
    breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
    xlab = "Probability p", col = "lightblue", border = "black"
  )
}

dev.off()


# ----- Run Simulations -----

n_sim <- 500
n_samples <- 1000
img_size <- 16
beta_effect <- 0.1
b_effect <- 0.2
b_sparsity <- 0.1
seed <- 42

cat("Simulating Data for", n_sim, "iterations of Simulation 1 ... \n")
tic()
sim1_data <- run_sim1(n_sim, n_samples, img_size, beta_effect, seed)
toc()

cat("Simulating Data for", n_sim, "iterations of Simulation 2 ... \n")
tic()
sim2_data <- run_sim2(n_sim, n_samples, img_size, b_effect, b_sparsity, seed)
toc()

sim1_1iter <- run_sim1(1, n_samples, img_size, beta_effect, seed)
sim2_1iter <- run_sim2(1, n_samples, img_size, b_effect, b_sparsity, seed)

filename <- paste0("simulated_data_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  n_sim, n_samples, img_size, beta_effect, b_effect, b_sparsity, seed,
  sim1_data, sim2_data, sim1_1iter, sim2_1iter,
  file = file.path(results_data_dir, filename)
)

filename2 <- paste0("simulated_data_1iter_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  n_sim, n_samples, img_size, beta_effect, b_effect, b_sparsity, seed,
  sim1_1iter, sim2_1iter,
  file = file.path(results_data_dir, filename2)
)


