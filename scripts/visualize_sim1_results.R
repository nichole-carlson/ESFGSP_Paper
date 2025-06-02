# Visualization plan for simulation results.
#
# Simulation 1 is performed in pixel space.
#
# For the true simulated data, we will show:
#   - Group mean difference at each pixel (heatmap)
#   - Group mean difference at each frequency (scatterplot)
#   - True coefs in pixel space (heatmap)
#   - True coefs transformed to freq space (scatterplot vs. eigenvalue order)
#
# For the fitted models, each iteration includes four models:
#   - Pixel space: fitted with lambda.min and lambda.1se
#   - Frequency space: fitted with lambda.min and lambda.1se
#
# For each model, we will generate:
#   - A heatmap of coefficients in pixel space
#   - A scatterplot of the corresponding coefficients in freq space
#   - A heatmap of p<0.05 in pixel space
#   - A scatterplot of the corresponding p<0.05 in freq space
#
# Note: The .rds file (~1GB) must be downloaded before use.
# A copy is available at: P:/ESFGSPproject/DataLibrary

library(rprojroot)
library(ggplot2)
library(patchwork)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- "" # temporary .rds location

source(file.path(proj_dir, "R", "summarize_results.R"))

# Read the .rds file saving all iterations
res <- readRDS(file.path(data_dir, "sim1_combined_results.rds"))

x_arr <- res$x
y_arr <- res$y
beta_vec <- res$beta
e <- res$e
hparams <- res$hparams
auc_acc_df <- res$auc_acc
coef_arr <- res$coefs
p_arr <- res$pvals


# ---------- Group Mean Difference ----------
# For 1 iteration
x_1 <- x_arr[, , 1]
y_1 <- y_arr[1, ]
# In the pixel space
p_group_mean_pixel <- vector_to_heatmap(
  colMeans(x_1[y_1 == 1, ]) - colMeans(x_1[y_1 == 0, ])
)
# In the freq space
x_index <- seq_along(beta_vec) / length(beta_vec)
x_freq_1 <- x_1 %*% e
p_group_mean_freq <- plot_scatter(
  x = x_index,
  y = colMeans(x_freq_1[y_1 == 1, ]) - colMeans(x_freq_1[y_1 == 0, ])
)


# ---------- Ture Coefficients Visualization ----------
# True coefs transformed to freq space: coef_freq = t(e) %*% coef_pixel
b_vec <- as.vector(t(e) %*% matrix(beta_vec, ncol = 1))

# Heatmap for true beta
p_true_beta <- vector_to_heatmap(beta_vec)

# Scatter plot of b_vec with x-axis as scaled index: i/256 for i = 1 to 256
x_index <- seq_along(b_vec) / length(b_vec)
p_true_b <- plot_scatter(x = x_index, y = b_vec)


# ---------- Average Estimated Coefs ----------
# Average estimated beta, fit on pixel space
beta_pixel <- coef_arr[, , "pixel", ]
beta_pixel_avg <- apply(beta_pixel, MARGIN = c(1, 2), FUN = mean)
p_est_beta_pixel <- {
  lim <- range(beta_pixel_avg) # lower and upper

  plot_list <- lapply(colnames(beta_pixel_avg), function(l) {
    .vec <- beta_pixel_avg[, l]
    vector_to_heatmap(.vec, lim) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}

# Average estimated b, fit in pixel space
b_pixel_avg <- t(e) %*% beta_pixel_avg
p_est_b_pixel <- {
  plot_list <- lapply(colnames(b_pixel_avg), function(l) {
    .vec <- b_pixel_avg[, l]
    plot_scatter(x_index, .vec) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}

# Average estimated b, fit in freq space
b_freq <- coef_arr[, , "freq", ]
b_freq_avg <- apply(b_freq, MARGIN = c(1, 2), FUN = mean)
p_est_b_freq <- {
  plot_list <- lapply(colnames(b_freq_avg), function(l) {
    .vec <- b_freq_avg[, l]
    plot_scatter(x_index, .vec) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}

# Average estimated beta, fit in freq space
beta_freq_avg <- e %*% b_freq_avg
p_est_beta_freq <- {
  lim <- range(beta_freq_avg)

  plot_list <- lapply(colnames(beta_freq_avg), function(l) {
    .vec <- beta_freq_avg[, l]
    vector_to_heatmap(.vec, lim) +
      labs(subtitle = l) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  })
  patchwork::wrap_plots(plot_list, ncol = 2)
}
