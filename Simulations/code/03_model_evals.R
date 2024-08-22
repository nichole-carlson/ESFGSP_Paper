library(table1)
library(ggplot2)
library(gridExtra)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load model fitting results -----
load(file = file.path(results_data_dir, "model_metrics_240821.RData"))
load(file = file.path(results_data_dir, "sim1_1iter_240815.RData"))

my.render.cont <- function(x) {
  with(stats.default(x), sprintf("%.3f (%0.3f)", MEAN, SD))
}


# ----- Create table for AUC and accuracy -----
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim1_auc_acc$pixel),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim1_auc_acc$freq),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim2_auc_acc$pixel),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim2_auc_acc$freq),
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

# For beta (coef in the pixel space), use heatmaps
# mean of estimated coefs
p_sim1_beta_est_min <- plot_heatmap(colMeans(sim1_coefs$pixel$lambda_min))
p_sim1_beta_est_1se <- plot_heatmap(colMeans(sim1_coefs$pixel$lambda_1se))
p_sim2_beta_est_min <- plot_heatmap(colMeans(sim2_coefs$pixel$lambda_min))
p_sim2_beta_est_1se <- plot_heatmap(colMeans(sim2_coefs$pixel$lambda_1se))

png(
  file = file.path(fig_dir, "beta_estimates.png"),
  width = 1200, height = 800, res = 150
)
grid.arrange(
  p_sim1_beta_est_min, p_sim1_beta_est_1se,
  p_sim2_beta_est_min, p_sim2_beta_est_1se,
  ncol = 2
)
dev.off()

# percentage of significant p-values
p_beta_pvals_sim1 <- plot_heatmap(colMeans(sim1_pvals$pixel < 0.05), c(0, 1))
p_beta_pvals_sim2 <- plot_heatmap(colMeans(sim2_pvals$pixel < 0.05), c(0, 1))

png(
  file = file.path(fig_dir, "perc_sign_pvals_beta.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(p_beta_pvals_sim1, p_beta_pvals_sim2, ncol = 2)
dev.off()


# For b (coef in freq space), use scatterplot against order of eigenvalues
x_cov <- sim1_1iter$meta_data$x_cov
eig_vals <- eigen(x_cov)$values
eig_vals_order <- order(eig_vals) / length(eig_vals)

# mean of estimated coefs
b_coef_data <- rbind(
  data.frame(
    e = eig_vals_order, b = colMeans(sim1_coefs$freq$lambda_min),
    sim = "sim1", l = "min"
  ),
  data.frame(
    e = eig_vals_order, b = colMeans(sim1_coefs$freq$lambda_1se),
    sim = "sim1", l = "1se"
  ),
  data.frame(
    e = eig_vals_order, b = colMeans(sim2_coefs$freq$lambda_min),
    sim = "sim2", l = "min"
  ),
  data.frame(
    e = eig_vals_order, b = colMeans(sim2_coefs$freq$lambda_1se),
    sim = "sim2", l = "1se"
  )
)

p_b_est <- ggplot(b_coef_data, aes(x = e, y = b)) +
  geom_point(color = "blue") +
  ylim(-0.1, 0.5) +
  labs(
    x = "order of eigenvalues",
    y = "mean estimated coefs"
  ) +
  theme_minimal() +
  facet_wrap(~ sim + l, ncol = 2) +
  theme(strip.text = element_blank())

png(
  file = file.path(fig_dir, "b_estimates.png"),
  width = 1200, height = 800, res = 150
)
p_b_est
dev.off()

# percentage of significant p-values
b_pval_data <- rbind(
  data.frame(
    e = eig_vals_order, p = colMeans(sim1_pvals$freq < 0.05), sim = "sim1"
  ),
  data.frame(
    e = eig_vals_order, p = colMeans(sim2_pvals$freq < 0.05), sim = "sim2"
  )
)

p_b_pvals <- ggplot(b_pval_data, aes(x = e, y = p)) +
  geom_point(color = "blue") +
  ylim(0.5, 1) +
  labs(
    x = "order of eigenvalues",
    y = "percentage of significant p-values"
  ) +
  theme_minimal() +
  facet_wrap(~sim, ncol = 2) +
  theme(strip.text = element_blank())

png(
  file = file.path(fig_dir, "perc_sign_pvals_b.png"),
  width = 1200, height = 500, res = 150
)
p
dev.off()
