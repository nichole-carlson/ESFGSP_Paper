library(table1)
library(ggplot2)
library(gridExtra)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load model fitting results -----
load(file = file.path(results_data_dir, "model_metrics_240815.RData"))
load(file = file.path(results_data_dir, "sim1_1iter_240815.RData"))

my.render.cont <- function(x) {
  with(stats.default(x), sprintf("%.3f (%0.3f)", MEAN, SD))
}


# ----- Create table for AUC and accuracy -----
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim1_auc_acc$pixel,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim1_auc_acc$freq,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim2_auc_acc$pixel,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim2_auc_acc$freq,
  render.continuous = my.render.cont
)


# ----- Visualize estimated coefficients -----

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
      low = "white",
      high = "black",
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

# For coefficients in the pixel space (beta), visualize by heatmap
png(
  file = file.path(fig_dir, "perc_sign_pvals_pixel_sim1.png"),
  width = 1600, height = 1200, res = 150
)
plot_heatmap(
  colMeans(sim1_pvals$pixel < 0.05),
  c(0, 1)
)
dev.off()

png(
  file = file.path(fig_dir, "perc_sign_pvals_pixel_sim2.png"),
  width = 1600, height = 1200, res = 150
)
plot_heatmap(
  colMeans(sim2_pvals$pixel < 0.05),
  c(0, 1)
)
dev.off()

# For coefficients in the frequency space (b), use scatterplot against order
# of eigenvalues
x_cov <- sim1_1iter$meta_data$x_cov
eig_vals <- eigen(x_cov)$values

png(
  file = file.path(fig_dir, "perc_sign_pvals_freq_sim1.png"),
  width = 1600, height = 1200, res = 150
)
plot(
  order(eig_vals) / length(eig_vals),
  colMeans(sim1_pvals$freq) < 0.05,
  xlab = "Order of Eigenvalues",
  ylab = "Percentage of Significant P-values"
)
dev.off()

png(
  file = file.path(fig_dir, "perc_sign_pvals_freq_sim2.png"),
  width = 1600, height = 1200, res = 150
)
plot(
  order(eig_vals) / length(eig_vals),
  colMeans(sim2_pvals$freq) < 0.05,
  xlab = "Order of Eigenvalues",
  ylab = "Percentage of Significant P-values"
)
dev.off()
