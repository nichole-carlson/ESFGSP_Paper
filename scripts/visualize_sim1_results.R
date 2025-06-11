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
# A copy of data is available at: P:/ESFGSPproject/DataLibrary

library(rprojroot)
library(ggplot2)
library(patchwork)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- "/Volumes/alzheimersdisease/ESFGSPproject/DataLibrary"

# vector_to_heatmap, plot_scatter
source(file.path(proj_dir, "R", "summarize_results.R"))
# transform_data, transform_coef
source(file.path(proj_dir, "R", "transformations.R"))

# Read in .rds file saves all iteration fits, and data from iter1
iter1 <- readRDS(file.path(data_dir, "data_015.rds"))
res <- readRDS(file.path(data_dir, "sim1_combined_results.rds"))

x <- iter1$x
y <- iter1$y
beta <- iter1$beta
e <- iter1$e
hparams <- iter1$hparams

auc_acc_df <- res$auc_acc
coef_arr <- res$coefs
p_arr <- res$pvals


# ---------- Group Mean Difference ----------
# In the pixel space
p_group_mean_pixel <- vector_to_heatmap(
  colMeans(x[y == 1, ]) - colMeans(x[y == 0, ])
)
# In the freq space
index <- seq_along(beta) / length(beta)
x_freq <- transform_data(x, e, to_freq = TRUE)
p_group_mean_freq <- plot_scatter(
  x = index,
  y = colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ])
)


# ---------- Ture Coefficients Visualization ----------
b <- transform_coef(beta, e, to_freq = TRUE)

# Heatmap for true beta
p_true_beta <- vector_to_heatmap(beta)

# Scatter plot of b with x-axis as scaled index: i/256 for i = 1 to 256
index <- seq_along(b) / length(b)
p_true_b <- plot_scatter(x = x_index, y = b)


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
b_pixel_avg <- transform_coef(beta_pixel_avg, e)
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
beta_freq_avg <- transform_coef(b_freq_avg, e)
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
