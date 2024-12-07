# ----- Simulation Plan -----



# Load required libraries
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")



# ----- Functions to simulate data -----

# Given the spatial adjacency matrix C, eigendecompose MCM, where M is a
# centering matrix. Returns the ordered eigenvalues and corresponding
# eigenvectors.
#
# Args:
#   mat: The spatial adjacency matrix C.
# Returns:
#   a list contains vectors and values.
eigen_decomp <- function(mat) {
  n_cols <- ncol(mat)
  cent_mat <- diag(n_cols) - matrix(1, n_cols, n_cols) / n_cols
  eig_res <- eigen(cent_mat %*% mat %*% cent_mat, symmetric = TRUE)

  ord_idx <- order(eig_res$values, decreasing = TRUE)
  eig_vecs <- eig_res$vectors[, ord_idx]
  eig_vals <- eig_res$values[ord_idx]

  return(list(vectors = eig_vecs, values = eig_vals))
}

# Generate a 2-neighbor adjacency matrix.
# For any two pixels in a matrix, they are considered adjacent if they are
# connected directly through edges, or indirectly through a shared neighboring
# pixel.
#
# Args:
#  n_row: int, number of rows in a matrix of pixels (an image).
#  n_col: int, number of columns in an image.
#  d_max: float. Maximum distance between pixels to be considered as adjacent.
#         For 2-neighbor adj matrix, d_max = 2.
# Return:
#  A p * p matrix encodes pixel adjacent: 1 for connected, 0 for not.
# p = n_row * n_col.
generate_n_neighbor_matrix <- function(n_row, n_col, d_max) {
  # Generate all coordinates
  coords <- expand.grid(row = 1:n_row, col = 1:n_col)
  # Compute Euclidean distance between every pair of coordinates
  distances <- as.matrix(dist(coords))
  # Decided adjacency by threshold (d_max)
  adj_mat <- (distances <= d_max) * 1
  # Change diagonal values to be 0
  diag(adj_mat) <- 0

  return(adj_mat)
}

# Generate a correlation matrix with exponential decay. For any two pixels
# i and j, those their 2D Euclidean distnace is d, then their correlation
# is exp(-d).
#
# Args:
#  n_row: int, number of rows in a matrix of pixels (consider an image).
#  n_col, int, number of columns in a matrix of pixels.
generate_exp_corr_matrix <- function(n_row, n_col) {
  # Generate all coorediates
  coords <- expand.grid(row = 1:n_row, col = 1:n_col)
  # Compute Euclidean distance between every pair of coordinates
  distances <- as.matrix(dist(coords))
  # Calculate each pair's correlation
  exp_corr_mat <- exp(-distances)

  return(exp_corr_mat)
}

# Generate a vector with a degree of sparsity.
#
# Args:
#  vec_len: int, the length of the vector.
#  sparse_level: [0, 1], the percentage of non-zero elements.
#  effect_size: float, the non-zero value.
#  seed.
#
# Returns:
#  Vector.
generate_1d_sparse_vector <- function(vec_len, sparse_level, effect_size, seed = NULL) {
  vec <- rep(0, vec_len)
  nz_index <- sample(
    seq_len(vec_len),
    size = floor(vec_len * sparse_level), # the number of items to choose
    replace = FALSE
  )
  vec[nz_index] <- effect_size
  return(vec)
}

# Generate a vector with a degree of sparsity in the 2D space.
# Image this vector is the coefficent vector for the pixels of a 2D image. The
# coefficients for a specific sub-area of the image has non-zero values.
#
# Args:
#  n_row: int, the dimension of the 2D image.
#  n_col: int, the dimension of the 2D image.
#  effect_size: float, the non-zero value.
#  active_area: vector, the indices for row_start, row_end, col_start, col_end.
generate_2d_sparse_vector <- function(n_row, n_col, effect_size, active_area) {
  # Initilize vector with zeros
  vec <- rep(0, n_row * n_col)
  # Calculae indices of the active_area
  indices <- as.vector(
    matrix(seq_along(vec), nrow = n_row)[
      active_area[1]:active_area[2], # row start and end
      active_area[3]:active_area[4] # column start and end
    ]
  )
  vec[indices] <- effect_size

  return(vec)
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

run_sim2 <- function(n_iter, n_samples, img_size, b_effect_size, sparsity, seed = NULL) {
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
for (effect in beta_effects) {
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
for (effect in b_effects) {
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
