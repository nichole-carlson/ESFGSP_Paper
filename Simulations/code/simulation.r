# Simulate imaging data to see how models work
library(ggplot2)
library(reshape2)

# Defaultly do parallel
# number of cores used = cores detected - 1

# Step 0: Set Up Functions

# - 0.1 Functions for data generation:
#   - Generate group indicator: 1000 in group A; 1000 in group B
#   - Generate beta as a 16x16 matrix with:
#       5 in the center 8x8 block, 0 elsewhere
#   - Generate an exponential correlation matrix:
#       For pixels i, j, correlation = exp(-dist(i, j))
#       where dist(i, j) is their Euclidean distance
#   - Use the correlation matrix to generate epsilon as:
#       multivariate normal (dims: 100, 2000, 256)
#   - Broadcast and generate y

# - 0.2 Functions for matrix visualization:
#   - Heatmap: given a (256, ) vector, visualize it as a 16*16 matrix
#   - Boxplot: given the p-values for outer area, show the boxplots

# - 0.3 Indicies for center space and outer space
# - 0.4 Functions for train/test split
# - 0.5 Functions for performing LASSO (predictions) and perm_lasso (p-values)
# - 0.6 Functions for eigen decomposition
# - 0.7 Functions for p-value adjustment and summary over iterations


# Step 1: VBM Model
# - Fit linear model on each pixel, with the pixel value as the outcome and
#   group indicator as the covariate
# - Adjust p-values for multitesting

# Step 4: spVBM Model

# Step 5: LASSO Model
# - This time using the group indicator as the outcome and all pixel values
#   as covariates
# - Calculate p-values by either the LASSO projection method,
#   or using permutation tests

##########################################################################










# # 0.6 Functions for eigen decomposition
# Perform eigen decomposition on a matrix
# Args:
#   mat: Matrix. The covariance matrix to be decomposed.
# Returns:
#   A list containing the eigenvectors and eigenvalues of the matrix.
eigen_decomp <- function(mat) {
  n_cols <- ncol(mat)
  cent_mat <- diag(n_cols) - matrix(1, n_cols, n_cols) / n_cols
  eig_res <- eigen(cent_mat %*% mat %*% cent_mat, symmetric = TRUE)

  ord_idx <- order(eig_res$values, decreasing = TRUE)
  eig_vecs <- eig_res$vectors[, ord_idx]
  eig_vals <- eig_res$values[ord_idx]

  return(list(vectors = eig_vecs, values = eig_vals))
}



# Function to calculate the percentages of significant p-values
#
# Args:
#   pvals: A matrix of p-values with dimensions (n_sim, n_pixels).
#   method: A string specifying the p-value adjustment method (default is "bonferroni").
#   alpha: A numeric value specifying the significance threshold (default is 0.05).
#
# Returns:
#   A list containing:
#     - pvals: The original matrix of p-values.
#     - pvals_adj: The matrix of adjusted p-values.
#     - perc_orig: A vector of percentages of significant p-values (original) for each pixel.
#     - perc_adj: A vector of percentages of significant p-values (adjusted) for each pixel.
calc_sig_perc <- function(pvals, method = "bonferroni", alpha = 0.05) {
  # Apply the specified p-value adjustment method across each row
  pvals_adj <- t(apply(pvals, 1, p.adjust, method = method))

  # Calculate the percentage of p-values < alpha for each column
  perc_orig <- colSums(pvals < alpha) / nrow(pvals) * 100
  perc_adj <- colSums(pvals_adj < alpha) / nrow(pvals) * 100

  # Return a list with both percentages and the adjusted p-values
  list(
    pvals = pvals,
    pvals_adj = pvals_adj,
    perc_orig = perc_orig,
    perc_adj = perc_adj
  )
}


# # Model 1: VBM --------------------------------------------------------

# # This function is designed exclusively for parallel execution.
# # In each iteration, the function generates 2000 images. It then utilizes a
# # grouping indicator to estimate the pixel values for each image. Each image
# # consists of a total of 256 pixels.
# # The output is a matrix with dimensions n_iters (number of iterations) by
# # n_pixels (number of pixels), containing the p-values for each pixel.

# vbm_fsim <- function(i) {
#   # simulate data
#   simulated_data <- generate_data(center_effect = center_effect)

#   n_pixel <- ncol(simulated_data) - 1
#   pvals <- rep(NA, n_pixel)

#   for (j in seq_len(n_pixel)) {
#     pixel_value <- simulated_data[, j + 1]
#     group_ind <- simulated_data[, 1]
#     model <- lm(pixel_value ~ group_ind)
#     pvals[j] <- summary(model)$coefficients[2, 4]
#   }

#   return(pvals)
# }

# vbm_pkgs <- c()
# vbm_objs <- c()
# list_package <- c(gen_data_pkgs, vbm_pkgs)

# # Run the function
# set.seed(42)
# tic()
# vbm_pvals <- simWrapper(
#   n_sim = 1000,
#   f_sim = vbm_fsim,
#   list_export = c(gen_data_objs, vbm_objs, "list_package"),
#   list_package = list_package
# )
# toc()

# # Adjust p-values for multiple testing
# # Calculate the perc of p-values < 0.05 for each pixel

# vbm_results <- calc_pval_adj(vbm_pvals)

# # Visualize the p-values as an image. Each pixel is the number of significant
# # p-values across iterations.
# image_vbm_pvals <- plot_matrix(vbm_results$perc, c(0, 100))
# image_vbm_pvals_corr <- plot_matrix(vbm_results$perc_corr, c(0, 100))

# ggsave(file.path(image_path, "vbm_pvals.png"), plot = image_vbm_pvals)
# ggsave(file.path(image_path, "vbm_pvals_corr.png"), plot = image_vbm_pvals_corr)

# # Visualize the p-values from the outer area with histogram
# image_vbm_boxplot <- plot_boxplot(
#   vbm_results$perc[e_indices],
#   vbm_results$perc_corr[e_indices]
# )

# ggsave(
#   file.path(image_path, "vbm_boxplots.png"),
#   plot = image_vbm_boxplot,
#   width = 8, height = 4, dpi = 300
# )



# # spVBM ------------------------------------------------------------
# # In this simulation, there is no subject-level non-spatial random effects
# # because each subject has only a single slice.

# # get the eigenvectors, space=1 means no approximation for the spacing
# # between the coordinates
# # eigen_vecs <- mymeigen2D(
# #   coords = df_long[, c("x", "y")], id = df_long$image_id, space = 1
# # )
# # tic()
# # spvbm_fit <- myresf_vc(
# #   y = df_long[["pixel_value"]],
# #   x = df_long[, "group_ind"],
# #   xgroup = factor(df_long$image_id),
# #   meig = eigen_vecs
# # )
# # toc()
# # # saveRDS(spvbm_fit, file = "spvbm_fit.rds")
# # # readRDS(file = "spvbm_fit.rds")

# # spvbm_coefs <- spvbm_fit$b_vc
# # spvbm_coefs_2d <- matrix(rowSums(spvbm_coefs), 16, 16)
# # image(
# #   t(spvbm_coefs_2d),
# #   col = gray(seq(1, 0, length = 256)),
# #   main = "2D Visualization of Estimated Parameters"
# # )
# # plot.new()


# # Model 3: LASSO -------------------------------------------------------------
# # For each iteration,
# # group_indicator as the outcome;
# # y (2000, 256) as the predictors;

# # Generate a sample dataset. See whether using all pixels will cause perfect
# # separation.
# set.seed(121)
# df <- generate_data(center_effect = center_effect)
# x <- df[, -1]
# y <- df[, 1]
# perform_lasso(x, y, p_train = 0.8) # use all pixels cause perfect separation

# # search how many pixels in the LASSO model will cause perfect separation

# i <- 1 # number of pixels from each side
# m <- 0 # maximum AUC/accuracy
# res <- matrix(NA, nrow = 0, ncol = 4)
# while (m < 1) {
#   ia <- sample(c_indices, i) + 1
#   ib <- sample(e_indices, i) + 1
#   x <- df[, c(ia, ib)]
#   y <- df[, 1]
#   res <- rbind(res, perform_lasso(x, y, p_train = 0.8))
#   m <- max(res[i, ])
#   i <- i + 1
# }
# # when use two points from each side, AUC/accuracy equals 1.
# rm(df, x, y, ia, ib, i, m, res)

# # This function is designed exclusively for parallel execution.

# lasso_fsim <- function(i) {
#   # simulate data
#   simulated_data <- generate_data(center_effect = center_effect)
#   x <- simulated_data[, -1]
#   y <- simulated_data[, 1]

#   # fit lasso
#   pvals <- perm_lasso(x, y, n_perm = 100)

#   return(pvals)
# }

# lasso_pkgs <- c("glmnet", "pROC")
# lasso_objs <- c("perm_lasso")
# list_package <- c(gen_data_pkgs, lasso_pkgs)

# # Run the function
# set.seed(42)
# tic()
# lasso_pvals <- simWrapper(
#   n_sim = 1000,
#   f_sim = lasso_fsim,
#   list_export = c(gen_data_objs, lasso_objs, "list_package"),
#   list_package = list_package
# )
# toc()
# # save(
# #   freq_pvals,
# #   file = paste0("lasso_pvals_", format(Sys.time(), "%y%m%d"), ".RData")
# # )
# load(file = "lasso_pvals_240522.RData")

# # Adjust p-values for multiple testing
# # Calculate the perc of p-values < 0.05 for each pixel
# lasso_results <- calc_pval_adj(lasso_pvals)

# # Visualize the p-values as an image. Each pixel is the number of significant
# # p-values across iterations.
# image_lasso_pvals <- plot_matrix(lasso_results$perc, c(0, 100))
# image_lasso_pvals_corr <- plot_matrix(lasso_results$perc_corr, c(0, 100))

# ggsave(file.path(image_path, "lasso_pvals.png"), plot = image_lasso_pvals)
# ggsave(file.path(image_path, "lasso_pvals_corr.png"), plot = image_lasso_pvals_corr)

# # Visualize the p-values from the outer area with histogram
# image_lasso_boxplot <- plot_boxplot(
#   lasso_results$perc[e_indices],
#   lasso_results$perc_corr[e_indices]
# )

# ggsave(
#   file.path(image_path, "lasso_boxplots.png"),
#   plot = image_lasso_boxplot,
#   width = 8, height = 4, dpi = 300
# )



# # Model 4: Frequency ------------------------------------------------------

# # Predict group_ind using image projected on the frequency domain
# # Find a fix correlation matrix, such as exp corr mat
# # Do the transformation use
# #   1. positive eigenvalues only
# #   2. all eigenvalues

# # Similar to LASSO model, evaluate model performance by perdiction accuracy
# # and p-values.

# # Generate a sample dataset. See how many eigenvectors will cause perfect sep.
# set.seed(121)
# df <- generate_data(center_effect = center_effect)
# x <- df[, -1]
# y <- df[, 1]
# corr_mat <- exp_corr_mat(ncol(x))
# eig_comp <- eig_decomp(corr_mat)

# i <- 2 # of positive eigenvectors to use
# m <- 0 # max AUC/accuracy
# res <- matrix(NA, nrow = 0, ncol = 4)
# while (m < 1) {
#   x_trans <- x %*% eig_comp$eigenvectors[, seq_len(i)]
#   res <- rbind(res, perform_lasso(x_trans, y, p_train = 0.8))
#   m <- max(res[i - 1, ])
#   i <- i + 1
# }
# #      min_acc min_auc 1se_acc 1se_auc
# # [1,]  0.4875     0.5  0.4875     0.5
# # [2,]  0.4875     0.5  0.4875     0.5   use 3 eigenvectors
# # [3,]  1.0000     1.0  1.0000     1.0

# rm(df, x, y, x_trans, i, m, res)

# freq_fsim <- function(i) {
#   # simulate data
#   simulated_data <- generate_data(center_effect = center_effect)
#   x <- simulated_data[, -1]
#   y <- simulated_data[, 1]

#   # build exponential correlation matrix
#   corr_mat <- exp_corr_mat(ncol(x))

#   # eigendecompose
#   eig_comp <- eig_decomp(corr_mat)

#   # eigen transpose in two ways
#   pos_eig_vals <- eig_comp$eigenvalues > 0
#   x_trans_pos <- x %*% eig_comp$eigenvectors[, pos_eig_vals]
#   x_trans_all <- x %*% eig_comp$eigenvectors

#   # fit lasso
#   pvals_pos <- perm_lasso(x_trans_pos, y, n_perm = 100)
#   pvals_all <- perm_lasso(x_trans_all, y, n_perm = 100)

#   # Combine results
#   results <- cbind(pvals_all, pvals_pos)

#   return(results)
# }

# freq_pkgs <- c("glmnet")
# freq_objs <- c("exp_corr_mat", "eig_decomp", "perm_lasso")
# list_package <- c(gen_data_pkgs, freq_pkgs)

# # Run the function
# set.seed(42)
# tic()
# freq_pvals <- simWrapper(
#   n_sim = 100,
#   f_sim = freq_fsim,
#   list_export = c(gen_data_objs, freq_objs, "list_package"),
#   list_package = list_package
# )
# toc()
# # save(
# #   freq_pvals,
# #   file = paste0("freq_pvals_", format(Sys.time(), "%y%m%d"), ".RData")
# # )
# load(file = "freq_pvals_240520.RData")

# freq_pvals <- freq_pvals[, seq_len(256)]

# # Adjust p-values for multiple testing
# # Calculate the perc of p-values < 0.05 for each pixel
# freq_results <- calc_pval_adj(freq_pvals)

# # Visualize the p-values as an image. Each pixel is the number of significant
# # p-values across iterations.
# image_freq_pvals <- plot_matrix(freq_results$perc, c(0, 100))
# image_freq_pvals_corr <- plot_matrix(freq_results$perc_corr, c(0, 100))

# ggsave(file.path(image_path, "freq_pvals.png"), plot = image_freq_pvals)
# ggsave(file.path(image_path, "freq_pvals_corr.png"), plot = image_freq_pvals_corr)

# # Visualize the p-values from the outer area with histogram
# image_freq_boxplot <- plot_boxplot(
#   freq_results$perc[e_indices],
#   freq_results$perc_corr[e_indices]
# )

# ggsave(
#   file.path(image_path, "freq_boxplots.png"),
#   plot = image_freq_boxplot,
#   width = 8, height = 4, dpi = 300
# )
