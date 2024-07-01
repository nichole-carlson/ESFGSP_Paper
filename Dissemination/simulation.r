# Simulate imaging data to see how models work

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

# setup_packages <- function(packages) {
#   for (package in packages) {
#     if (!require(package, character.only = TRUE, quietly = TRUE)) {
#       cat("Installing package:", package, "\n")
#       install.packages(package)
#       if (!require(package, character.only = TRUE, quietly = TRUE)) {
#         cat("Failed to install package:", package, "\n")
#       }
#     }
#   }
# }

# List of packages to check and install if necessary
# packages_to_install <- c(
#   "MASS", "stats", "tidyverse", "hdi", "tictoc", "viridis", "reshape2",
#   "parallel", "glmnet", "pROC"
# )

# setup_packages(packages_to_install)

# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")

parent_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Dissemination"
# image_path <- file.path(proj_path, "Figures")
source(file.path(parent_dir, "simWrapper.r"))

# 0.1 Functions for data generation

# # Generate an exponential correlation matrix for a square matrix with n_pixels
# exp_corr_mat <- function(n_pixel) {
#   n_col <- sqrt(n_pixel)

#   # given any index from 1 to n^2, calculate its 2D indices
#   convert_1d_to_2d <- function(index, n_col) {
#     # Adjust for zero-based indexing
#     index_zero_based <- index - 1

#     # Calculate row and column
#     row_index <- index_zero_based %/% n_col
#     column_index <- index_zero_based %% n_col

#     # Return the 2D coordinates (row, column)
#     return(c(row_index, column_index))
#   }

#   # Function to calculate the Euclidean distance between two points in a 2D grid
#   dist <- function(i, j, n_col) {
#     # Convert 1D indices to 2D coordinates
#     coord1 <- convert_1d_to_2d(i, n_col)
#     coord2 <- convert_1d_to_2d(j, n_col)

#     # Calculate Euclidean distance
#     distance <- sqrt((coord1[1] - coord2[1])^2 + (coord1[2] - coord2[2])^2)

#     return(distance)
#   }

#   vec_dist <- Vectorize(dist, c("i", "j"))

#   # Generate the correlation matrix
#   corr_mat <- outer(seq_len(n_pixel), seq_len(n_pixel), function(i, j) {
#     dist_ij <- vec_dist(i, j, n_col)
#     exp(-dist_ij)
#   })

#   return(corr_mat)
# }

# # The output is a 2000 * 257 matrix. The first column is the group indicator.
# # Each column of the remaining represents a pixel.

# generate_data <- function(
#     na = 1000, # n obs from group A
#     nb = 1000,
#     image_size = c(16, 16),
#     center_size = 8,
#     center_effect) {
#   n_row <- image_size[1]
#   n_col <- image_size[2]
#   n_pixel <- n_row * n_col
#   # generate group indicator
#   group_ind <- c(rep(1, na), rep(0, nb)) # (2000, )

#   # group effect by pixel
#   beta <- matrix(0, nrow = n_row, ncol = n_col)
#   index_st <- (n_row - center_size) %/% 2 + 1
#   index_end <- index_st + center_size - 1
#   beta[index_st:index_end, index_st:index_end] <- center_effect
#   beta <- as.vector(beta) # (256, )

#   # multivariate normal error
#   corr_mat <- exp_corr_mat(n_pixel)

#   # generate errors (2000, 256)
#   epsilon <- MASS::mvrnorm(
#     n = na + nb, mu = rep(0, n_pixel), Sigma = corr_mat
#   )

#   # generate y
#   y <- outer(group_ind, beta) + epsilon

#   return(cbind(group_ind, y))
# }

# gen_data_objs <- c("generate_data", "exp_corr_mat")
# gen_data_pkgs <- c("MASS")

# 0.2 Functions for matrix visualization
# Input should be a (256, ) vector, representing a 1D 16*16 image.
plot_matrix <- function(vec, value_limits = c()) {
  # Convert 1D vector to 2D matrix
  mat <- matrix(vec, 16, 16, byrow = TRUE)

  # Convert the 2D matrix to a long data frame: x, y, value
  data_long <- reshape2::melt(mat)

  # Ensure value_limits is of length 2
  if (length(value_limits) != 2) {
    stop("value_limits must be a vector of length 2, like c(low, high)")
  }

  # Create the plot using ggplot2
  plot <- ggplot(data_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(
      low = "white", high = "black",
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

# # # Decide the strength of center effect by visualization
# # set.seed(121)
# # df1 <- generate_data(center_effect = 4)
# # plot_matrix(df1[78, -1], range(df1))

# # df2 <- generate_data(center_effect = 5)
# # plot_matrix(df2[78, -1], range(df2))

# # image_4 <- plot_matrix(df1[78, -1], range(df1, df2))
# # image_5c <- plot_matrix(df2[1001, -1], range(df1, df2)) # w/o center effect
# # image_5 <- plot_matrix(df2[78, -1], range(df1, df2)) # w/ center effect

# # ggsave(file.path(image_path, "ex_image_4.png"), plot = image_4)
# # ggsave(file.path(image_path, "ex_image_5.png"), plot = image_5)
# # ggsave(file.path(image_path, "ex_image_5c.png"), plot = image_5c)

# # Set up global center effect
# center_effect <- 5
# gen_data_objs <- c(gen_data_objs, "center_effect")

# # vec1 should be a vector of p-values before multiple testing adjustment,
# # vec2 should be after adjustment
# plot_boxplot <- function(vec1, vec2, y_range = c()) {
#   # Create a data frame for plotting
#   data <- data.frame(
#     perc = c(vec1, vec2),
#     grp = factor(
#       c(
#         rep("Unadjusted", length(vec1)),
#         rep("Adjusted", length(vec2))
#       ),
#       levels = c("Unadjusted", "Adjusted")
#     )
#   )

#   # Create boxplots
#   p <- ggplot(data, aes(x = grp, y = perc)) +
#     geom_boxplot() +
#     facet_wrap(~grp, scales = "free_x") +
#     labs(y = "Percentage of Significant P-Values", x = "") +
#     theme_minimal() +
#     theme(
#       plot.background = element_rect(fill = "white"),
#       panel.background = element_rect(fill = "white"),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.x = element_blank()
#     )

#   # Apply y-axis range if specified
#   if (length(y_range) == 2) {
#     p <- p + scale_y_continuous(
#       limits = y_range,
#       breaks = scales::pretty_breaks(n = 10)
#     )
#   } else {
#     p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
#   }

#   return(p)
# }

# # 0.3 Indicies for center space and outer space
# c_indices <- {
#   c_rows <- 5:12
#   as.vector(outer(c_rows, c_rows, FUN = function(i, j) (i - 1) * 16 + j))
# }
# e_indices <- {
#   all_indices <- 1:256
#   setdiff(all_indices, c_indices)
# }


# 0.5 Functions for performing LASSO (predictions) and perm_lasso (p-values)
# select optimal lambda via cross-validation
# evaluate performance via accuracy and AUC
perform_lasso <- function(x, y, p_train, seed = 42) {
  # Split the dataset
  set.seed(seed)
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]

  # Perform cross-validation to find the optimal lambda values
  cv_model <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
  lambda_min <- cv_model$lambda.min
  lambda_1se <- cv_model$lambda.1se

  # Function to evaluate model performance
  eval_perf <- function(l) {
    mod <- glmnet(x_train, y_train, alpha = 1, lambda = l, family = "binomial")
    preds <- predict(mod, newx = x_test, type = "response")[, 1]
    preds_bin <- ifelse(preds > 0.5, 1, 0)

    acc <- mean(preds_bin == y_test)
    auc <- pROC::auc(pROC::roc(y_test, preds))

    c(acc = acc, AUC = auc)
  }

  # Get performance metrics for both lambda.min and lambda.1se
  perf_min <- eval_perf(lambda_min)
  perf_1se <- eval_perf(lambda_1se)

  # Combine results into a coherent matrix
  results <- c(perf_min, perf_1se)
  results <- matrix(results, nrow = 1)
  colnames(results) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")

  return(results)
}

# Estimate p-value for each pixel using permutation test
# The output is a 1 * n_pixel matrix with the value as p-values
perm_lasso <- function(x, y, n_perm, seed = 42) {
  set.seed(seed)

  cv_lasso <- function(x, y) {
    cv_fit <- cv.glmnet(x, y, alpha = 1)
    best_lambda <- cv_fit$lambda.min
    coefs <- coef(cv_fit, s = "lambda.min")[-1, 1]
    return(coefs)
  }

  # get original coef estimates
  orig_coefs <- cv_lasso(x, y)

  # perform permutations
  perm_coefs <- matrix(NA, n_perm, length(orig_coefs))

  for (i in seq_len(n_perm)) {
    y_perm <- sample(y)
    perm_coefs[i, ] <- cv_lasso(x, y_perm)
  }

  # calculate p-values
  p_vals <- matrix(colMeans(abs(perm_coefs) > abs(orig_coefs)), nrow = 1)
  return(p_vals)
}

# # 0.6 Functions for eigen decomposition
# eig_decomp <- function(C) {
#   p <- ncol(C)
#   M <- diag(p) - matrix(1, p, p) / p
#   eig_data <- eigen(M %*% C %*% M, symmetric = TRUE)

#   order_idx <- order(eig_data$values, decreasing = TRUE)
#   eig_vecs <- eig_data$vectors[, order_idx]
#   eig_vals <- eig_data$values[order_idx]

#   return(list(
#     eigenvectors = eig_vecs,
#     eigenvalues = eig_vals
#   ))
# }

# 0.7 Functions for p-value adjustment and summary over iterations
# Input should be a matrix, # of rows equals # of simulations
calc_pval_adj <- function(pvals) {
  # Apply Bonferroni correction across each row
  pvals_corr <- t(apply(pvals, 1, p.adjust, method = "bonferroni"))

  # Calculate the percentage of p-values < 0.05 for each column
  perc_orig <- colSums(pvals < 0.05) / nrow(pvals) * 100
  perc_corr <- colSums(pvals_corr < 0.05) / nrow(pvals) * 100

  # Return a list with both percentages
  list(
    pvals = pvals,
    pvals_corr = pvals_corr,
    perc = perc_orig,
    perc_corr = perc_corr
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
