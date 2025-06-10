# List of required packages
packages <- c("MASS")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load packages
sapply(packages, require, character.only = TRUE)



# ----- Functions for simulating data -----

# Generate a matrix of samples from a multivariate normal distribution.
simulate_mvn_samples <- function(n_sample, cov_matrix, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_features <- nrow(cov_matrix)
  samples <- MASS::mvrnorm(
    n_sample,
    mu = rep(0, n_features),
    Sigma = cov_matrix
  )
  samples <- matrix(samples, nrow = n_sample) # make sure output is matrix
  return(samples)
}


# Simulates y ~ bin(p) in pixel space.
#
# Args:
#   data_mat is n x p matrix. param_vec is either p x 1 matrix or length-p vec.
#
#   from_pixel_space = TRUE if data_mat and param_vec are in pixel space.
#
#   e is p x p matrix. Required if from_pixel_space is FALSE.
simulate_pixel_outcome <- function(data_mat,
                                   param_vec,
                                   from_pixel_space = TRUE,
                                   e = NULL,
                                   seed = NULL) {
  if (is.null(dim(param_vec))) {
    param_vec <- matrix(param_vec, ncol = 1)
  } else if (ncol(param_vec) != 1) {
    stop("param_vec must be a column vector (ncol = 1).")
  }

  if (!from_pixel_space) {
    if (is.null(e)) stop("e must be provided when not from pixel space.")
    data_mat <- data_mat %*% t(e) # map freq to pixel
    param_vec <- t(e) %*% param_vec # map freq to pixel
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  eta <- drop(data_mat %*% param_vec)
  prob <- plogis(eta)
  y <- rbinom(length(prob), 1, prob)
  return(y)
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
  hparams$eigen_vectors <- eigen_vectors
  hparams$eigen_values <- eigen_decomposition$values

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
  hparams$eigen_vectors <- eigen_vectors
  hparams$eigen_values <- eigen_decomposition$values

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
