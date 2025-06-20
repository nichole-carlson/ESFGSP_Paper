# simulate_data.R

library(MASS)


# --------------------
# Generate exponential correlation matrix for 16x16 image grid
generate_exp_corr <- function(n_row, n_col, rate) {
  # Generate all coorediates
  coords <- expand.grid(row = 1:n_row, col = 1:n_col)
  # Compute Euclidean distance between every pair of coordinates
  distances <- as.matrix(dist(coords))

  corr_mat <- exp(-rate * distances)

  return(corr_mat)
}


# Generate a 2-neighbor adjacency matrix for a grid of pixels.
gen_n_neighbor_adj_mat <- function(n_row, n_col, d_max) {
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


# Eigendecomposition of M*C*M where M is centering matrix, C is adjacency mat
eigen_decomp_mcm <- function(adj_mat) {
  n <- ncol(adj_mat)
  # Center matrix: M = I - (1/n) * 11'
  cent_adj_mat <- diag(n) - adj_matrix(1 / n, n, n)
  # mcm
  mcm <- cent_adj_mat %*% adj_mat %*% cent_adj_mat
  mcm <- (mcm + t(mcm)) / 2 # for numerical stability

  eigen(mcm, symmetric = TRUE)
}


# Generate multivariate normal data
generate_image_data <- function(n, cov_matrix, mu = NULL) {
  p <- nrow(cov_matrix)

  if (is.null(mu)) {
    mu <- rep(0, p)
  }

  image_data <- MASS::mvrnorm(n, mu = mu, Sigma = cov_matrix)

  return(image_data)
}


# Generate block coefficients (for simulation 1)
# grid_dims accepts c(nrow, ncol) for image dimensions
# block_dims the same for block dimensions
generate_block_coefs <- function(grid_dims, block_dims, beta_value) {
  if (length(grid_dims) != 2 || length(block_dims) != 2) {
    stop("grid_dims and block_dims must be length 2 vectors")
  }

  if (any(block_dims > grid_dims)) {
    stop("block dimensions cannot exceed grid dimensions")
  }

  grid_nrow <- grid_dims[1]
  grid_ncol <- grid_dims[2]
  block_nrow <- block_dims[1]
  block_ncol <- block_dims[2]

  beta_2d <- matrix(0, nrow = grid_nrow, ncol = grid_ncol)

  # Calculate central block position
  start_row <- (grid_nrow - block_nrow) %/% 2 + 1
  start_col <- (grid_ncol - block_ncol) %/% 2 + 1
  end_row <- start_row + block_nrow - 1
  end_col <- start_col + block_ncol - 1

  beta_2d[start_row:end_row, start_col:end_col] <- beta_value

  # Convert to vector (rowwise)
  beta_1d <- as.vector(t(beta_2d))

  return(beta_1d)
}


# Generate sparse coefficients (for simulation 2)
generate_sparse_coefs <- function(n, p_nonzero, beta_value) {
  beta <- rep(0, n)
  indicies <- sample(1:n, size = floor(n * p_nonzero))
  beta[indicies] <- beta_value

  return(beta)
}


# Generate binary outcome from logistic model
generate_outcomes <- function(x, beta, beta0 = 0) {
  if (!is.matrix(beta)) {
    beta <- matrix(beta, ncol = 1)
  }

  if (ncol(x) != nrow(beta)) {
    stop("Dimension mismatch: ncol(x) must equal length(beta)")
  }

  n <- nrow(x)

  eta <- beta0 + x %*% beta
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n, 1, p)

  return(y)
}


# Transform data between spaces
transform_data <- function(x, e, to_freq = TRUE) {
  if (to_freq) {
    x %*% e
  } else {
    x %*% t(e)
  }
}


# Transform coefs between spaces
transform_coef <- function(coefs, e, to_freq = TRUE) {
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


# ---------- Function to use ----------
simulate_data <- function(n_sample, cov_matrix, coef_vec, adj_matrix, on_freq) {
  # Simulate items
  x_orig <- generate_image_data(n_sample, cov_matrix)
  transform_mat <- eigen_decomp_mcm(adj_matrix)$vectors
  x_trans <- transform_data(x_orig, transform_mat, !on_freq)
  coef_trans <- transform_coef(coef_vec, transform_mat, !on_freq)
  y <- generate_outcomes(x_orig, coef_vec)

  # Prepare for return
  x_arr <- cbind(x_orig, x_trans)
  colnames(x_arr) <- c("orig", "trans")

  coef_arr <- cbind(coef_vec, coef_trans)
  colnames(coef_arr) <- c("orig", "trans")

  return(list(
    x = x_arr,
    coef = coef_arr,
    outcome = y,
    on_freq = on_freq,
    transform_mat = transform_mat
  ))
}
