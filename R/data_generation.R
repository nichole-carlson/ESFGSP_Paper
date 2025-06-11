# data_generation.R
# Functions for generating data, correlation matrix, coefficients, and outcomes

library(MASS)


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


# Generate multivariate normal data
generate_image_data <- function(n, cov_matrix, mu = NULL) {
  p <- nrow(cov_matrix)

  if (is.null(mu)) {
    mu <- rep(0, p)
  }

  image_data <- MASS::mvrnorm(n, mu = mu, Sigma = cov_matrix)

  return(iamge_data)
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
