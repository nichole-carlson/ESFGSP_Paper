# transformations.R
# Functions for E matrix generation and space transformations


# Eigendecomposition of M*C*M where M is centering matrix
eigen_decomp_mcm <- function(mat) {
  n <- ncol(mat)
  # Center matrix: M = I - (1/n) * 11'
  cent_mat <- diag(n) - matrix(1 / n, n, n)
  # mcm
  mcm <- cent_mat %*% mat %*% cent_mat
  mcm <- (mcm + t(mcm)) / 2 # for numerical stability

  eigen(mcm, symmetric = TRUE)
}

# Generate E matrix from correlation structure
generate_e_matrix <- function(corr_matrix) {
  # Adjacency matrix
  adj_matrix <- corr_matrix - diag(nrow(corr_matrix))

  eigen_decomp_mcm(adj_matrix)$vectors
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
