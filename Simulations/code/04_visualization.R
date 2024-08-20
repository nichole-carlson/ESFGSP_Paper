library(ggplot2)
library(gridExtra)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
fig_dir <- file.path(simulation_dir, "results", "figures")
res_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load Simulated Data for a Single Iteration -----
load(file = file.path(res_data_dir, "sim1_1iter_240815.RData"))
load(file = file.path(res_data_dir, "sim2_1iter_240815.RData"))

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
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
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


# ----- Visualize group mean difference in X/X_freq between y=0 and y=1 -----
visual_group_diff <- function(sim_data) {
  data <- sim_data$data
  x <- data$x
  x_freq <- data$x_freq
  y <- data$y

  # Calculate mean differences based on y assignment
  mean_diff_x <- colMeans(x[y == 1, ]) - colMeans(x[y == 0, ])
  mean_diff_x_freq <- colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ])

  # Define the range of X and X_freq
  # x_range <- c(floor(min(x)), ceiling(max(x)))
  # x_freq_range <- c(floor(min(x_freq)), ceiling(max(x_freq)))

  # Generate heatmaps
  p1 <- plot_heatmap(mean_diff_x)
  p2 <- plot_heatmap(mean_diff_x_freq)

  grid.arrange(p1, p2, ncol = 2)
}

png(
  file = file.path(fig_dir, "group_diff_sim1.png"),
  width = 1200, height = 500, res = 150
)
visual_group_diff(sim1_1iter)
dev.off()


png(
  file = file.path(fig_dir, "group_diff_sim2.png"),
  width = 1200, height = 500, res = 150
)
visual_group_diff(sim2_1iter)
dev.off()


# ----- Visualize beta and b using selected effect size -----
visualize_coefs <- function(sim_data) {
  meta_data <- sim_data$meta_data
  b <- meta_data$b
  beta <- meta_data$beta

  p1 <- plot_heatmap(beta)
  p2 <- plot_heatmap(b)

  grid.arrange(p1, p2, ncol = 2)
}

png(
  file = file.path(fig_dir, "coefs_sim1.png"),
  width = 1200, height = 500, res = 150
)
visualize_coefs(sim1_1iter)
dev.off()

png(
  file = file.path(fig_dir, "coefs_sim2.png"),
  width = 1200, height = 500, res = 150
)
visualize_coefs(sim2_1iter)
dev.off()


# ----- Visualize the top three and the bottom eigenvectors -----
# the eigenvectors are the same across sim1 and sim2
visualize_eigvecs <- function(sim_data) {
  x_cov <- sim_data$meta_data$x_cov
  eig_res <- eigen(x_cov)
  eig_vecs <- eig_res$vectors[, order(eig_res$values, decreasing = TRUE)]
  p1 <- plot_heatmap(eig_vecs[, 1])
  p2 <- plot_heatmap(eig_vecs[, 2])
  p3 <- plot_heatmap(eig_vecs[, 3])
  p4 <- plot_heatmap(eig_vecs[, 256])

  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

png(
  file = file.path(fig_dir, "top_bottom_eigvecs.png"),
  width = 1200, height = 800, res = 150
)
visualize_eigvecs(sim1_1iter)
dev.off()
