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
data_dir <- "/Volumes/alzheimersdisease/ESFGSPproject/DataLibrary/sim1"
fig_dir <- file.path(proj_dir, "results", "figures", "sim1")

# vector_to_heatmap, plot_scatter, create_heatmap_plots, create_scatter_plots
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

p <- length(beta)
index <- seq_len(p) / p
x_freq <- transform_data(x, e, to_freq = TRUE)
b <- transform_coef(beta, e, to_freq = TRUE)


# ---------- Group Mean Difference ----------
# In the pixel space
p_group_mean_pixel <- vector_to_heatmap(
  colMeans(x[y == 1, ]) - colMeans(x[y == 0, ])
)
# In the freq space
p_group_mean_freq <- plot_scatter(
  x = index,
  y = colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ])
)


# ---------- Ture Coefficients Visualization ----------
# Heatmap for true beta
p_true_beta <- vector_to_heatmap(beta)

# Scatter plot of b with x-axis as scaled index: i/256 for i = 1 to 256
p_true_b <- plot_scatter(x = index, y = b)


# ---------- Average Estimated Coefs ----------
# Average est coefs fit in pixel space
beta_pixel_avg <- apply(coef_arr[, , "pixel", ], MARGIN = c(1, 2), FUN = mean)
# Average est coefs fit in freq space
b_freq_avg <- apply(coef_arr[, , "freq", ], MARGIN = c(1, 2), FUN = mean)

p_est_beta_pixel <- create_heatmap_plots(beta_pixel_avg)
p_est_b_pixel <- create_scatter_plots(transform_coef(beta_pixel_avg, e), index)
p_est_b_freq <- create_scatter_plots(b_freq_avg, index)
p_est_beta_freq <- create_heatmap_plots(transform_coef(b_freq_avg, e))



# ---------- Percentage of p<0.05 Visualization ----------
# Percentage of p<0.05 for the model fit on pixel space
pixel_pval <- apply(
  p_arr[, , "pixel", ],
  MARGIN = c(1, 2),
  FUN = function(x) mean(x < 0.05)
)
freq_pval <- apply(
  p_arr[, , "freq", ],
  MARGIN = c(1, 2),
  FUN = function(x) mean(x < 0.05)
)

p_pixel_pval <- create_heatmap_plots(pixel_pval)
p_freq_pval <- create_scatter_plots(freq_pval, index)


# ---------- Save figures ----------
ggplot2::ggsave(
  file.path(fig_dir, "group_mean_diff_pixel.pdf"),
  p_group_mean_pixel,
  width = 4, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "group_mean_diff_freq.pdf"),
  p_group_mean_freq,
  width = 4, height = 4
)

ggplot2::ggsave(
  file.path(fig_dir, "true_coef_pixel.pdf"),
  p_true_beta,
  width = 4, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "true_coef_freq.pdf"),
  p_true_b,
  width = 4, height = 4
)

ggplot2::ggsave(
  file.path(fig_dir, "est_coef_beta_pixel.pdf"),
  p_est_beta_pixel,
  width = 8, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_b_pixel.pdf"),
  p_est_b_pixel,
  width = 8, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_beta_freq.pdf"),
  p_est_beta_freq,
  width = 8, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_b_freq.pdf"),
  p_est_b_freq,
  width = 8, height = 4
)

ggplot2::ggsave(
  file.path(fig_dir, "pval_pixel.pdf"),
  p_pixel_pval,
  width = 8, height = 4
)
ggplot2::ggsave(
  file.path(fig_dir, "pval_freq.pdf"),
  p_freq_pval,
  width = 8, height = 4
)
