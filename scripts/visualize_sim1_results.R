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

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- "" # temporary .rds location

source(file.path(proj_dir, "R", "summarize_results.R"))

# Read the .rds file saving all iterations
res <- readRDS(file.path(data_dir, "sim1_combined_results.rds"))

x_list <- res$x
y_list <- res$y
beta_vec <- res$beta
e <- res$e
hparams <- res$hparams
auc_acc_df <- res$auc_acc
coefs_pvals_df <- res$coefs_pvals


# ---------- Group Mean Difference ----------
# For 1 iteration
x_1 <- x_list[[1]]
y_1 <- y_list[[1]]
# In the pixel space
p_group_mean_pixel <- vector_to_heatmap(
  colMeans(x_1[y_1 == 1, ]) - colMeans(x_1[y_1 == 0, ])
)
# In the freq space
x_index <- seq_along(y_1) / length(y_1)
x_freq_1 <- x_1 %*% e
p_group_mean_freq <- plot_scatter(
  x = x_index,
  y = colMeans(x_freq_1[y_1 == 1, ]) - colMeans(x_freq_1[y_1 == 0, ])
)


# ---------- Real Coefficients Visualization ----------
# true coefs transformed to freq space: coef_freq = t(e) %*% coef_pixel
b_vec <- as.vector(t(e) %*% matrix(beta_vec, ncol = 1))

# heatmap for true beta
p_beta_heat <- vector_to_heatmap(beta_vec)

# Scatter plot of b_vec with x-axis as scaled index: i/256 for i = 1 to 256
x_index <- seq_along(b_vec) / length(b_vec)
p_b_scatter <- plot_scatter(x = x_index, y = b_vec)
