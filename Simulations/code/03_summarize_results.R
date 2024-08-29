library(ggplot2)
library(gridExtra)
library(table1)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
fig_dir <- file.path(simulation_dir, "results", "figures")
results_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load model fitting results and simulated data -----
load(file = file.path(results_data_dir, "simulated_data_1iter_240826.RData"))
load(file = file.path(results_data_dir, "model_metrics_240826.RData"))

sim_data_list <- list(sim1_1iter, sim2_1iter)
sim_coefs_list <- list(sim1_coefs, sim2_coefs)
sim_pvals_list <- list(sim1_pvals, sim2_pvals)
lambda_options <- c("lambda_min", "lambda_1se")
eig_vals_ordered <- {
  eig_vals <- sim1_1iter$meta_data$eig_vals
  order(eig_vals) / length(eig_vals)
}


# ----- Functions for visualization -----
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

# specifically defined for plotting eigenvalues vs other values
plot_scatterplot <- function(x, y, xlab = NULL, ylab = NULL) {
  ymin <- floor(min(y, na.rm = TRUE))
  ymax <- ceiling(max(y, na.rm = TRUE))

  dat <- data.frame(x = x, y = y)
  ggplot(dat, aes(x = x, y = y)) +
    geom_point(color = "blue") +
    labs(x = xlab, y = ylab) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    theme_minimal() +
    theme(strip.text = element_blank())
}



# ----- Visualize actual coefs and group means -----
p_group_mean_diff <- list()
for (i in seq_len(length(sim_data_list))) {
  dat <- sim_data_list[[i]]$data
  x <- dat$x
  x_freq <- dat$x_freq
  y <- dat$y

  mean_diff_x <- colMeans(x[y == 1, ]) - colMeans(x[y == 0, ])
  mean_diff_x_freq <- colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ])

  p_group_mean_diff[[paste0("sim", i)]][["beta"]] <- plot_heatmap(mean_diff_x)
  p_group_mean_diff[[paste0("sim", i)]][["b"]] <- plot_scatterplot(
    eig_vals_ordered, mean_diff_x_freq,
    "Order of eigenvalues", "Group mean diff"
  )
}

p_actual_coefs <- list()
for (i in seq_len(length(sim_data_list))) {
  meta_data <- sim_data_list[[i]]$meta_data
  beta <- meta_data$beta
  b <- meta_data$b

  p_actual_coefs[[paste0("sim", i)]][["beta"]] <- plot_heatmap(beta)
  p_actual_coefs[[paste0("sim", i)]][["b"]] <- plot_scatterplot(
    eig_vals_ordered, b,
    "Order of eigenvalues", "Coefficient values"
  )
}

for (i in seq_len(length(sim_data_list))) {
  filename1 <- paste0("actual_coefs_sim", i, ".png")
  png(file = file.path(fig_dir, filename1), width = 1200, height = 500)
  grid.arrange(grobs = p_actual_coefs[[paste0("sim", i)]], ncol = 2)
  dev.off()

  filename2 <- paste0("group_mean_diff_sim", i, ".png")
  png(file = file.path(fig_dir, filename2), width = 1200, height = 500)
  grid.arrange(grobs = p_group_mean_diff[[paste0("sim", i)]], ncol = 2)
  dev.off()
}


# ----- Create table for AUC and accuracy -----
my.render.continuous <- function(x) {
  with(stats.default(x), sprintf("%.3f (%0.3f)", MEAN, SD))
}

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
# For beta (coef in the pixel space), use heatmaps
p_est_coefs <- list()
for (i in seq_len(length(sim_coefs_list))) {
  dat <- sim_coefs_list[[i]]
  for (j in lambda_options) {
    p_est_coefs[[paste0("sim", i)]]$beta[[j]] <- plot_heatmap(colMeans(dat$pixel[[j]]))
    p_est_coefs[[paste0("sim", i)]]$b[[j]] <- plot_scatterplot(eig_vals_ordered, colMeans(dat$freq[[j]]))
  }
}

png(
  file = file.path(fig_dir, "beta_estimates.png"),
  width = 1200, height = 800, res = 150
)
grid.arrange(grobs = c(p_est_coefs$sim1$beta, p_est_coefs$sim2$beta), ncol = 2)
dev.off()

# For b, use scatterplot
png(
  file = file.path(fig_dir, "b_estimates.png"),
  width = 1200, height = 800, res = 150
)
grid.arrange(grobs = c(p_est_coefs$sim1$b, p_est_coefs$sim2$b), ncol = 2)
dev.off()


# ----- Visualize significant p-values -----
# percentage of significant p-values
p_pvals <- list()
for (i in seq_len(length(sim_pvals_list))) {
  dat <- sim_pvals_list[[i]]
  p_pvals[[paste0("sim", i)]]$beta <- plot_heatmap(colMeans(dat$pixel < 0.05), c(0, 1))
  p_pvals[[paste0("sim", i)]]$b <- plot_scatterplot(eig_vals_ordered, colMeans(dat$freq < 0.05), c(0, 1))
}

png(
  file = file.path(fig_dir, "perc_sign_pvals_beta.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals$sim1$beta, p_pvals$sim2$beta), ncol = 2)
dev.off()

png(
  file = file.path(fig_dir, "perc_sign_pvals_b.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals$sim1$b, p_pvals$sim2$b), ncol = 2)
dev.off()
