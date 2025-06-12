# List of required packages
packages <- c("ggplot2", "reshape2", "patchwork")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load packages
sapply(packages, require, character.only = TRUE)


# Convert a vector to a heatmap
#
# This function reshapes a 1D vector into a 2D matrix (assuming the length is
# a perfect square) and plots it as a heatmap.
vector_to_heatmap <- function(vec, value_limits = NULL) {
  # Ensure the vector length is a perfect square
  len <- length(vec)
  img_size <- as.integer(sqrt(len))
  if (img_size^2 != len) {
    stop("The length of 'vec' must be a perfect square.")
  }

  # Reshape vector to matrix
  mat <- matrix(vec, nrow = img_size, byrow = TRUE)

  # Convert matrix to long format for ggplot
  data_long <- reshape2::melt(mat)

  # Validate value_limits
  if (!is.null(value_limits) && length(value_limits) != 2) {
    stop("value_limits must be a numeric vector of length 2.")
  }

  # Create the heatmap plot
  ggplot(data_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = value_limits, oob = scales::squish
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(x = NULL, y = NULL)
}

# Plot a scatterplot of two numeric variables
plot_scatter <- function(x, y, xlab = NULL, ylab = NULL) {
  # Calculate y-axis limits, ensuring no NAs
  y_range <- range(y, na.rm = TRUE)
  y_limits <- c(floor(y_range[1]), ceiling(y_range[2]))

  # Create a data frame for plotting
  dat <- data.frame(x = x, y = y)

  # Generate the scatterplot
  ggplot(dat, aes(x = x, y = y)) +
    geom_point(color = "blue", size = 2) + # Size for better visibility
    labs(x = xlab, y = ylab) +
    scale_y_continuous(limits = y_limits) +
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank() # Cleaner plot
    )
}


# Plot each column of input matrix as a heatmap
create_heatmap_plots <- function(data_matrix, shared_limits = TRUE) {
  lim <- if (shared_limits) range(data_matrix) else NULL

  plot_list <- lapply(colnames(data_matrix), function(l) {
    .vec <- data_matrix[, l]
    vector_to_heatmap(.vec, lim) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}


# Plot each column of input matrix as a scatterplot
create_scatter_plots <- function(data_matrix, x_vals) {
  plot_list <- lapply(colnames(data_matrix), function(l) {
    .vec <- data_matrix[, l]
    plot_scatter(x_vals, .vec) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}
