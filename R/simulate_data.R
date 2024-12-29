library(MASS)


# ----- Functions for simulating data -----

# Perform eigendecomposition of MCM, where M is a centering matrix.
#
# Args:
#   mat: Matrix. The spatial adjacency matrix (C) to be decomposed.
#
# Returns:
#   A list containing:
#     - vectors: Matrix of eigenvectors ordered by descending eigenvalues.
#     - values: Vector of eigenvalues in descending order.
#
# Notes:
#   The centering matrix (M) is computed as (I - 1/n * 1 * 1^T), where n is
#   the number of columns in the matrix. This function applies the
#   transformation MCM and returns the eigenvalues and eigenvectors of the
#   resulting matrix.
eigen_decomp <- function(mat) {
  n_cols <- ncol(mat)
  cent_mat <- diag(n_cols) - matrix(1, n_cols, n_cols) / n_cols
  eig_res <- eigen(cent_mat %*% mat %*% cent_mat, symmetric = TRUE)

  ord_idx <- order(eig_res$values, decreasing = TRUE)
  eig_vecs <- eig_res$vectors[, ord_idx]
  eig_vals <- eig_res$values[ord_idx]

  return(list(vectors = eig_vecs, values = eig_vals))
}


# Generate a 2-neighbor adjacency matrix for a grid of pixels.
#
# Args:
#   n_row: Integer. Number of rows in the pixel grid (e.g., an image).
#   n_col: Integer. Number of columns in the pixel grid.
#   d_max: Numeric. Maximum distance for pixels to be considered adjacent.
#          For a 2-neighbor adjacency matrix, set d_max = 2.
#
# Returns:
#   A p x p matrix (where p = n_row * n_col) encoding pixel adjacency:
#     - 1 indicates that two pixels are connected.
#     - 0 indicates that two pixels are not connected.
#
# Notes:
#   Pixels are considered adjacent if they are directly connected via edges or
#   indirectly connected through a shared neighbor within the given distance.
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


# Generate a correlation matrix with exponential decay.
#
# For any two pixels i and j, their correlation is exp(-d), where d is their
# 2D Euclidean distance.
#
# Args:
#   n_row: Integer. Number of rows in the pixel grid (e.g., an image).
#   n_col: Integer. Number of columns in the pixel grid.
#
# Returns:
#   A p x p correlation matrix (where p = n_row * n_col) with values determined
#   by the exponential decay of the 2D Euclidean distance between pixels.
generate_exp_corr_matrix <- function(n_row, n_col) {
  # Generate all coorediates
  coords <- expand.grid(row = 1:n_row, col = 1:n_col)
  # Compute Euclidean distance between every pair of coordinates
  distances <- as.matrix(dist(coords))
  # Calculate each pair's correlation
  exp_corr_mat <- exp(-distances)

  return(exp_corr_mat)
}


# Generate a diagonal matrix with values that decay non-linearly.
#
# The diagonal values decay rapidly at the beginning and slow down towards the
# end, following the formula:
#   d_i = start_value * exp(-k * log(1 + b * i))
# where:
#   - k is the decay rate, computed based on start_value and end_value.
#   - b is a constant controlling the rate of initial decay.
#
# Args:
#   n_row: Integer. Number of rows in the pixel grid (e.g., an image).
#   n_col: Integer. Number of columns in the pixel grid.
#   start_value: Numeric. Initial value for the diagonal (default: 6).
#   end_value: Numeric. Final value for the diagonal (default: 0.5).
#   b: Numeric. Controls the rate of initial decay (default: 0.1).
#
# Returns:
#   A numeric diagonal matrix of size p x p (where p = n_row * n_col), with
#   the specified decay pattern.
#
# Notes:
#   The decay factor k is calculated to ensure the diagonal values smoothly
#   transition from start_value to end_value across n elements.
generate_diag_corr_matrix <- function(
    n_row, n_col, start_value = 6, end_value = 0.5, b = 0.1) {
  p <- n_row * n_col
  # Compute the decay rate k
  k <- log(start_value / end_value) / log(1 + b * (p - 1))
  # Generate the diagonal values with non-linear decay
  diag_vals <- start_value * exp(-k * log(1 + b * seq(0, p - 1)))
  # Create a diagonal matrix with the computed values
  diag_mat <- diag(diag_vals)
  return(diag_mat)
}


# Calculate the indices for the center area of a matrix.
#
# The center area is automatically calculated as half the dimensions of the
# input matrix (rounded down for odd dimensions). Returns the row and column
# indices for the center area as a single vector.
#
# Args:
#   n_row: Integer. Number of rows in the matrix.
#   n_col: Integer. Number of columns in the matrix.
#
# Returns:
#   A numeric vector containing the start and end indices for the rows and
#   columns of the center area in the format:
#   c(row_start, row_end, col_start, col_end).
#
# Notes:
#   For matrices with odd dimensions, the center area will be slightly smaller
#   as it is calculated by rounding down the dimensions.
calculate_center_indices <- function(n_row, n_col) {
  # Calculate half dimensions for the center area
  center_rows <- floor(n_row / 2)
  center_cols <- floor(n_col / 2)

  # Calculate start and end indices for rows
  row_start <- floor((n_row - center_rows) / 2) + 1
  row_end <- row_start + center_rows - 1

  # Calculate start and end indices for columns
  col_start <- floor((n_col - center_cols) / 2) + 1
  col_end <- col_start + center_cols - 1

  # Return indices as a vector
  return(c(row_start, row_end, col_start, col_end))
}


# Simulate a sparse vector with a specified level of sparsity.
#
# Args:
#   vec_len: Integer. Length of the vector to be simulated.
#   sparse_level: Numeric (0 to 1). Proportion of non-zero elements.
#   effect_size: Numeric. Value assigned to non-zero elements in the vector.
#   seed: Optional. Integer seed for reproducibility.
#
# Returns:
#   A numeric vector of length vec_len, with the specified level of sparsity.
#
# Notes:
#   Non-zero elements are randomly assigned within the vector based on
#   sparse_level.
simulate_1d_sparse_vector <- function(
    vec_len, sparse_level, effect_size, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  vec <- rep(0, vec_len)
  nz_index <- sample(
    seq_len(vec_len),
    size = floor(vec_len * sparse_level), # the number of items to choose
    replace = FALSE
  )
  vec[nz_index] <- effect_size

  return(vec)
}

# Generate a sparse vector representing coefficients in a 2D pixel space.
#
# This vector represents a 2D image where coefficients in a specified sub-area
# are assigned non-zero values.
#
# Args:
#   n_row: Integer. Number of rows in the 2D image.
#   n_col: Integer. Number of columns in the 2D image.
#   effect_size: Numeric. Value assigned to non-zero items in the active area.
#   active_area: Numeric vector. Indices specifying the sub-area where non-zero
#                values are applied.
#                Format: c(row_start, row_end, col_start, col_end).
#
# Returns:
#   A numeric vector of length n_row * n_col, with sparsity determined by the
#   specified active area.
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


# Generate a matrix of samples from a multivariate normal distribution.
#
# Args:
#   n_samples: Integer. Number of samples (rows) to generate.
#   cov_matrix: Matrix. Covariance matrix defining the multivariate normal
#               distribution.
#   seed: Optional. Integer seed for reproducibility.
#
# Returns:
#   A numeric matrix of samples with n_samples rows and columns equal to the
#   dimension of the covariance matrix (img_size^2 for a 2D pixel grid).
simulate_mvn_samples <- function(n_samples, cov_matrix, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_features <- nrow(cov_matrix)
  samples <- MASS::mvrnorm(
    n_samples,
    mu = rep(0, n_features),
    Sigma = cov_matrix
  )
  samples <- matrix(samples, nrow = n_samples) # make sure output is matrix
  return(samples)
}


# Simulate binary response variables based on a logistic model.
#
# Args:
#   x: Matrix or data frame. Predictor variables (n_samples x n_features).
#   beta: Numeric vector. Coefficients for the logistic model.
#   seed: Optional. Integer seed for reproducibility.
#
# Returns:
#   A numeric vector of binary response variables (0 or 1) generated using the
#   logistic model: P(y = 1) = 1 / (1 + exp(-x * beta)).
simulate_binary_response <- function(x, beta, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  probs <- 1 / (1 + exp(-as.matrix(x) %*% beta))
  binary_response <- rbinom(n = length(probs), size = 1, prob = probs)
  return(binary_response)
}


# This function generates datasets where predictors are simulated in the pixel
# space and transformed into the frequency space using eigendecomposition of a
# user-provided spatial adjacency matrix. Binary outcomes are generated based
# on sparse coefficients in the pixel space.
#
# Args:
#   n_samples: Integer. Number of samples (observations) per simulation.
#   n_row: Integer. Number of rows in the spatial grid (image dimensions).
#   n_col: Integer. Number of columns in the spatial grid.
#   effect_size: Numeric. Non-zero value assigned to coefficients in the active
#                region of the pixel space.
#   c_adj: Matrix. Spatial adjacency matrix used for eigendecomposition to
#          transform data between pixel and frequency spaces.
#   seed: Optional. Integer seed for reproducibility (default: NULL).
#
# Returns:
#   A list containing:
#     - beta: Numeric vector. Sparse coefficient vector in the pixel space.
#     - b: Numeric vector. Coefficient vector transformed to the frequency
#          space.
#     - data: List containing:
#         - x: Numeric matrix (n_samples x (n_row * n_col)). Predictors in the
#              pixel space.
#         - x_freq: Numeric matrix (n_samples x (n_row * n_col)). Predictors
#                   transformed to the frequency space.
#         - y: Numeric vector. Binary outcomes (0 or 1) generated from x and
#              beta.
run_pixel_to_freq_simulation <- function(
    n_samples, n_row, n_col, effect_size, c_adj, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # generate the exponential decay covariance matrix for X
  c_cov <- generate_exp_corr_matrix(n_row, n_col)

  # generate coefficient vector beta
  center_indices <- calculate_center_indices(n_row, n_col)
  beta <- generate_2d_sparse_vector(n_row, n_col, effect_size, center_indices)

  # perform eigendecomposition based on the spatial adjacency matrix given.
  eigen_decomposition <- eigen_decomp(c_adj)
  eigen_vectors <- eigen_decomposition$vectors

  # calculate the coefficient vector in the frequency space
  b <- t(eigen_vectors) %*% beta

  # simulate x and y
  x <- simulate_mvn_samples(n_samples, c_cov)

  # calculate corresponding predictors in the frequency space
  x_freq <- x %*% eigen_vectors

  # simulate outcome y
  y <- simulate_binary_response(x, beta)

  # all arguments
  hparams <- as.list(match.call())[-1]

  simulations <- list(
    beta = beta,
    b = b,
    x = x,
    x_freq = x_freq,
    y = y,
    hparams = hparams
  )

  return(simulations)
}


# This function simulates a dataset where the predictors are generated in the
# freq space and transformed back to the pixel space using eigendecomposition
# of a spatial adjacency matrix. Binary outcomes are generated based on sparse
# coefficients in the frequency space.
#
# Args:
#   n_samples: Integer. Number of samples (observations) per simulation.
#   n_row: Integer. Number of rows in the spatial grid (image dimensions).
#   n_col: Integer. Number of columns in the spatial grid.
#   sparse_level: Numeric (0 to 1). Proportion of non-zero coefficients in the
#                 sparse coefficient vector (b).
#   effect_size: Numeric. Value assigned to non-zero coefficients in b.
#   seed: Optional. Integer seed for reproducibility (default: NULL).
#
# Returns:
#   A list containing:
#     - beta: Numeric vector. Coefficient vector in the pixel space, derived
#             from the frequency space coefficients.
#     - b: Numeric vector. Sparse coefficient vector in the frequency space.
#     - data: List contains:
#         - x: Numeric matrix (n_samples x (n_row * n_col)). Predictors in the
#              pixel space.
#         - x_freq: Numeric matrix (n_samples x (n_row * n_col)). Predictors
#                   in the frequency space.
#         - y: Numeric vector. Binary outcomes (0 or 1) generated from x_freq
#              and b.
run_simulation_1b <- function(
    n_samples, n_row, n_col, sparse_level, effect_size, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # generate the diagonal covariance matrix for x
  c_cov <- diag(seq(6, 0.5, length.out = n_row * n_col))

  # generate coefficient vector b
  b <- simulate_1d_sparse_vector(
    vec_len = ncol(c_cov),
    sparse_level = sparse_level,
    effect_size = effect_size
  )

  # choose the spatial adjacency matrix and perform eigendecomposition
  c_adj <- c_cov
  eigen_decomposition <- eigen_decomp(c_adj)
  eigen_vectors <- eigen_decomposition$vectors

  # calculate the coefficient vector in the pixel space
  beta <- eigen_vectors %*% b

  # generate predictors in the frequency space
  x_freq <- simulate_mvn_samples(n_samples, c_cov)

  # transform it back to the pixel space
  x <- x_freq %*% t(eigen_vectors)

  # simulate outcome y
  y <- simulate_binary_response(x_freq, b)

  # save all arguments
  hparams <- as.list(match.call())[-1]

  simulations <- list(
    beta = beta,
    b = b,
    x = x,
    x_freq = x_freq,
    y = y,
    hparams = hparams
  )

  return(simulations)
}
