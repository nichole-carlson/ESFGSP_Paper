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

library(rprojroot)
library(ggplot2)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- "/scratch/alpine/sren@xsede.org/esfgsp"

source(file.path(proj_dir, "R", "summarize_results.R"))

# Read the .rds file saving all iterations
res <- readRDS(file.path(data_dir, "sim1_combined_results.rds"))

# Matrix of eigenvectors (decreasingly ordered)
e <- res$data[[1]]$e

# Simulation 1 true coefs in pixel space
beta_matrix <- matrix(0, nrow = 16, ncol = 16)
beta_matrix[5:12, 5:12] <- 0.1
beta_vec <- as.vector(beta_matrix)

# true coefs transformed to freq space: coef_freq = t(e) %*% coef_pixel
b_vec <- t(e) %*% matrix(beta_vec, ncol = 1) |>
  as.vector()
