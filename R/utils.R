# Perform eigendecomposition of MCM, where M is a centering matrix.
#
# Args:
#   mat: Matrix. The spatial adjacency matrix (C) to be decomposed.
#
# Returns:
#   A list containing:
#     - vectors: Matrix of eigenvectors ordered by decreasing eigenvalues.
#     - values: Vector of eigenvalues in decreasing order.
eigen_decomp_mcm <- function(mat) {
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
#   A p x p correlation matrix (where p = n_row * n_col).
gen_exp_corr_mat <- function(n_row, n_col) {
  # Generate all coorediates
  coords <- expand.grid(row = 1:n_row, col = 1:n_col)
  # Compute Euclidean distance between every pair of coordinates
  distances <- as.matrix(dist(coords))
  # Calculate each pair's correlation
  exp_corr_mat <- exp(-distances)

  return(exp_corr_mat)
}
