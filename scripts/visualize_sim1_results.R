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
source(file.path(proj_dir, "R", "simulate_data.R"))

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


# ---------- AUC/ACC ----------
auc_acc_df |>
  dplyr::group_by(space, lambda) |>
  dplyr::summarize(
    mean_auc = mean(auc),
    se_auc = sd(auc),
    mean_acc = mean(accuracy),
    se_acc = sd(accuracy)
  )


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
coefs_avg <- apply(coef_arr, 1:4, mean)
# Average est coefs fit in pixel space
p_est_beta_pixel <- create_heatmap_plots(t(coefs_avg["pixel", , "orig", ]))
p_est_b_pixel <- create_scatter_plots(t(coefs_avg["pixel", , "trans", ]), index)
# Average est coefs fit in freq space
p_est_b_freq <- create_scatter_plots(t(coefs_avg["freq", , "orig", ]), index)
p_est_beta_freq <- create_heatmap_plots(t(coefs_avg["freq", , "trans", ]))


# Goal: select the top b values under lambda.min and lambda.1se, project
# back to beta, visualize
est_b <- t(coefs_avg["freq", , "orig", ]) # (256, 2)

# Select the top n values by absolute value, set others to zero
select_top_sparse <- function(coefs, top_n) {
  # Get indices of top n values by absolute value
  top_indices <- order(abs(coefs), decreasing = TRUE)[1:top_n]

  # Create sparse version with only top n coefficients
  sparse_coefs <- rep(0, length(coefs))
  sparse_coefs[top_indices] <- coefs[top_indices]

  return(sparse_coefs)
}

lambda_min_plots <- purrr::map(seq(1, 20, 2), ~ {
  transformed_vec <- est_b[, 1] |>
    select_top_sparse(top_n = .x) |>
    transform_coef(e = e, to_freq = FALSE)

  vector_to_heatmap(transformed_vec) +
    ggplot2::ggtitle(paste("Top", .x, "coefficients"))
}) |>
  patchwork::wrap_plots(ncol = 5)

lambda_1se_plots <- purrr::map(seq(1, 20, 2), ~ {
  transformed_vec <- est_b[, 2] |>
    select_top_sparse(top_n = .x) |>
    transform_coef(e = e, to_freq = FALSE)

  vector_to_heatmap(transformed_vec) +
    ggplot2::ggtitle(paste("Top", .x, "coefficients"))
}) |>
  patchwork::wrap_plots(ncol = 5)

ggplot2::ggsave(
  file.path(fig_dir, "est_beta_freq_lmin_top20.pdf"),
  lambda_min_plots,
  width = 28, height = 16
)
ggplot2::ggsave(
  file.path(fig_dir, "est_beta_freq_l1se_top20.pdf"),
  lambda_1se_plots,
  width = 28, height = 16
)


# ---------- Percentage of p<0.05 Visualization ----------
# Percentage of p<0.05 for the model fit on pixel space
pvals_signif <- apply(p_arr, 1:4, function(x) mean(x < 0.05))

p_pval_pixel_orig <- create_heatmap_plots(t(pvals_signif["pixel", , "orig", ]))
p_pval_pixel_trans <- create_scatter_plots(
  t(pvals_signif["pixel", , "trans", ]), index
)
p_pval_freq_orig <- create_scatter_plots(
  t(pvals_signif["freq", , "orig", ]), index
)
p_pval_freq_trans <- create_heatmap_plots(t(pvals_signif["freq", , "trans", ]))



# ---------- Save figures ----------
ggplot2::ggsave(
  file.path(fig_dir, "group_mean_diff_pixel.pdf"),
  p_group_mean_pixel,
  width = 8, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "group_mean_diff_freq.pdf"),
  p_group_mean_freq,
  width = 8, height = 8
)

ggplot2::ggsave(
  file.path(fig_dir, "true_coef_pixel.pdf"),
  p_true_beta,
  width = 8, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "true_coef_freq.pdf"),
  p_true_b,
  width = 8, height = 8
)

ggplot2::ggsave(
  file.path(fig_dir, "est_coef_beta_pixel.pdf"),
  p_est_beta_pixel,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_b_pixel.pdf"),
  p_est_b_pixel,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_beta_freq.pdf"),
  p_est_beta_freq,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "est_coef_b_freq.pdf"),
  p_est_b_freq,
  width = 12, height = 8
)

ggplot2::ggsave(
  file.path(fig_dir, "pval_pixel_orig.pdf"),
  p_pval_pixel_orig,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "pval_pixel_trans.pdf"),
  p_pval_pixel_trans,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "pval_freq_orig.pdf"),
  p_pval_freq_orig,
  width = 12, height = 8
)
ggplot2::ggsave(
  file.path(fig_dir, "pval_freq_trans.pdf"),
  p_pval_freq_trans,
  width = 12, height = 8
)
