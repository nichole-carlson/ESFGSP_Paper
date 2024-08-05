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
library(gridExtra)

fig_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations/results/figures"


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

# Function to plot a 1D vector as a 16x16 heatmap
#
# Args:
#   vec: A numeric vector of length 256 representing a 1D 16x16 image.
#   value_limits: An optional numeric vector of length 2 specifying the value limits for the color scale (default is NULL).
#   title: A string specifying the title of the plot (default is "Heatmap").
#
# Returns:
#   A ggplot2 object representing the heatmap.
plot_heatmap <- function(vec, value_limits = NULL) {
  # Convert 1D vector to 2D matrix
  len <- length(vec)
  img_size <- as.integer(sqrt(len))
  mat <- matrix(vec, img_size, img_size, byrow = TRUE)

  # Convert the 2D matrix to a long data frame: x, y, value
  data_long <- reshape2::melt(mat)

  # Ensure value_limits is of length 2 if provided
  if (!is.null(value_limits) && length(value_limits) != 2) {
    stop("value_limits must be a vector of length 2, like c(low, high)")
  }

  # Create the plot using ggplot2
  plot <- ggplot(data_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(
      low = "white", high = "black",
      limits = value_limits,
      oob = scales::squish
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(x = "", y = "")

  return(plot)
}

generate_x_and_x_freq <- function(n_samples, img_size) {
  n_pixels <- img_size^2
  W <- gen_exp_corr(n_pixels) # correlation matrix for x
  V <- eigen(W)$vectors
  W_freq <- t(V) %*% W %*% V # correlation matrix for x_freq

  x_freq <- gen_x(n_samples, W_freq)
  x <- x_freq %*% t(V)

  list(x = x, x_freq = x_freq, V = V)
}

run_simulation1 <- function(n_iter, n_samples, img_size, beta_effect_size, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # run generate_x_and_x_freq function one time to get V, beta and b
  res <- generate_x_and_x_freq(n_samples, img_size)
  V <- res$V
  beta <- gen_beta(img_size, beta_effect_size)
  b <- beta_to_b(beta, V)

  result <- list(
    coefs = list(beta = beta, b = b, V = V),
    data = list()
  )

  for (i in seq_len(n_iter)) {
    res <- generate_x_and_x_freq(n_samples, img_size)
    x <- res$x
    x_freq <- res$x_freq
    y <- gen_y(x_freq, b)

    result$data[[i]] <- list(x = x, x_freq = x_freq, y = y)
  }
}

run_simulation2 <- function(n_iter, n_samples, img_size, b_effect_size, sparsity, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  res <- generate_x_and_x_freq(n_samples, img_size)
  V <- res$V
  b <- gen_b(img_size^2, sparsity, b_effect_size)
  beta <- b_to_beta(b, V)

  result <- list(
    coefs = list(beta = beta, b = b, V = V),
    data = list()
  )

  for (i in seq_len(n_iter)) {
    res <- generate_x_and_x_freq(n_samples, img_size)
    x <- res$x
    x_freq <- res$x_freq
    y <- gen_y(x_freq, b)

    result$data[[i]] <- list(x = x, x_freq = x_freq, y = y)
  }
}


# ----- Choose beta effect size and b effect size -----

# Function to compare beta effects for Simulation 1
compare_beta_effects <- function(effects, n_samples = 1000, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))
  set.seed(seed)

  res <- generate_x_and_x_freq(n_samples, img_size = 256)
  x <- res$x

  for (effect in effects) {
    beta <- gen_beta(img_size = 16, effect_size = effect)
    p <- 1 / (1 + exp(-(x %*% beta)))
    hist(p,
      breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
      xlab = "Probability p", col = "lightblue", border = "black"
    )
  }
}

# png(
#   file = file.path(fig_dir, "sim1_p_dist.png"),
#   width = 1600, height = 1200, res = 150
# )
# compare_beta_effects(c(1, 0.2, 0.1, 0.05, 0.01))
# dev.off()

# Function to compare b effects for Simulation 2
compare_b_effects <- function(effects, n_samples = 1000, sparsity = 0.1, seed = 42) {
  par(mfrow = c(ceiling(length(effects) / 2), 2))

  res <- generate_x_and_x_freq(n_samples, img_size = 256)
  x_freq <- res$x_freq

  for (effect in effects) {
    b <- gen_b(len = 256, sparsity, effect, seed)
    p <- 1 / (1 + exp(-(x_freq %*% b)))
    hist(p, breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1), xlab = "Probability p", col = "lightblue", border = "black")
  }
}

# png(
#   file = file.path(fig_dir, "sim2_p_dist.png"),
#   width = 1600, height = 1200, res = 150
# )
# compare_b_effects(c(1, 0.8, 0.6, 0.4, 0.2, 0.1))
# dev.off()



# ----- Run Simulations -----

# Global settings
# n_sim <- 500
# n_samples <- 1000
# img_size <- 16
# beta_effect <- 0.1
# b_effect <- 0.4
# b_sparsity <- 0.1
# seed <- 42


# Visualize beta and b using selected effect size
visualize_coefs <- function(img_size, beta_effect, b_effect, b_sparsity, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate global variables
  W <- gen_exp_corr(img_size^2)
  V <- eigen(W)$vectors
  W_freq <- t(V) %*% W %*% V

  beta <- gen_beta(img_size, beta_effect)
  b <- gen_b(img_size^2, b_sparsity, b_effect)

  # Visualize coefs
  # When assigning sparsity to beta
  p1 <- plot_heatmap(beta)
  p2 <- plot_heatmap(beta_to_b(beta, V))
  # When assigning sparsity to b
  p3 <- plot_heatmap(b)
  p4 <- plot_heatmap(b_to_beta(b, V))

  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

# coef_plots <- visualize_coefs(
#   img_size = img_size,
#   beta_effect = beta_effect,
#   b_effect = b_effect,
#   b_sparsity = b_sparsity,
#   seed = seed
# )
# ggsave(
#   file.path(fig_dir, paste0("coef_heatmaps_", format(Sys.Date(), "%y%m%d"), ".png")),
#   coef_plots
# )

# Simulate x and y for two sparsity patterns
# tic()
# cat("Simulating Data for", n_sim, "iterations ... \n")
# simulated_data <- simulate_data(
#   n_sim = n_sim,
#   n_samples = n_samples,
#   img_size = img_size,
#   beta_effect = beta_effect,
#   b_effect = b_effect,
#   b_sparsity = b_sparsity,
#   seed = seed
# )
# toc()

# Function to visualize the simulated data
visualize_simulated_data <- function(sim_data) {
  n_sim <- length(sim_data)
  plots <- list()
  for (i in seq_len(n_sim)) {
    name <- names(sim_data)[[i]]
    data <- sim_data[[name]]

    x <- data$x
    x_freq <- data$x_freq
    y <- data$y

    # Calculate mean differences based on y assignment
    mean_diff_x <- colMeans(x[y == 1, ]) - colMeans(x[y == 0, ])
    mean_diff_x_freq <- colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ])

    # Define the common range for both heatmaps
    common_range <- range(mean_diff_x, mean_diff_x_freq)

    # Generate heatmaps
    p1 <- plot_heatmap(mean_diff_x, common_range)
    p2 <- plot_heatmap(mean_diff_x_freq, common_range)

    plots[[length(plots) + 1]] <- p1
    plots[[length(plots) + 1]] <- p2
  }
  do.call(grid.arrange, c(plots, ncol = 2))
}

# # Extract data from the first simulation and visualize
# ggsave(
#   file.path(fig_dir, paste0("group_mean_image_", format(Sys.Date(), "%y%m%d"), ".png")),
#   visualize_simulated_data(simulated_data)
# )
