##########################################################################
# Simulate imaging data to see how models work

# Defaultly do parallel
# number of cores used = cores detected - 1

# Step 1: Simulate Data
# - Generate group indicator, 1000 in group A; 1000 in group B
# - Generate the beta as 16*16 matrix, it has 1 in the center 8*8 and 0 in
#    other areas.
# - Generate the exponential correlation matrix
#   - For two pixels i, j in the picture, their correlation are calculated
#     as exp(-dist(i, j)). dist(i, j) are defined as their Euclidean dist.
# - Use the corr matrix to generate the epsilon as a multivariate normal.
#    Epsilon should be (100, 2000, 256)
# - Broadcast and generate y

# Step 2: Visualization
# - Pick one image from the group A, another from group B
# - Visualize them by heatmap, check how obvious the center effect is

# Step 3: VBM Model
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

setup_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", package, "\n")
      install.packages(package)
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        cat("Failed to install package:", package, "\n")
      }
    }
  }
}

# List of packages to check and install if necessary
packages_to_install <- c(
  "MASS", "stats", "tidyverse", "hdi", "tictoc", "viridis", "reshape2",
  "parallel", "glmnet", "pROC"
)

setup_packages(packages_to_install)

# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")

source("/Users/siyangren/Documents/ESFGSP/simWrapper.r")



# 1 Define function to simulate data ------------------------------------
# The output is a 2000 * 257 matrix. The first column is the group indicator.
# Each column of the remaining represents a pixel.

generate_data <- function(
    na = 1000, # n obs from group A
    nb = 1000,
    image_size = c(16, 16),
    center_size = 8,
    center_effect) {
  n_row <- image_size[1]
  n_col <- image_size[2]
  n_pixel <- n_row * n_col
  # generate group indicator
  group_ind <- c(rep(1, na), rep(0, nb)) # (2000, )

  # group effect by pixel
  beta <- matrix(0, nrow = n_row, ncol = n_col)
  index_st <- (n_row - center_size) %/% 2 + 1
  index_end <- index_st + center_size - 1
  beta[index_st:index_end, index_st:index_end] <- center_effect
  beta <- as.vector(beta) # (256, )

  # multivariate normal error
  convert_1d_to_2d <- function(index, ncol) {
    # Adjust for zero-based indexing
    index_zero_based <- index - 1

    # Calculate row and column
    row_index <- index_zero_based %/% ncol
    column_index <- index_zero_based %% ncol

    # Return the 2D coordinates (row, column)
    return(c(row_index, column_index))
  }

  # Function to calculate the Euclidean distance between two points in a 2D grid
  dist <- function(i, j, ncol) {
    # Convert 1D indices to 2D coordinates
    coord1 <- convert_1d_to_2d(i, ncol)
    coord2 <- convert_1d_to_2d(j, ncol)

    # Calculate Euclidean distance
    distance <- sqrt((coord1[1] - coord2[1])^2 + (coord1[2] - coord2[2])^2)

    return(distance)
  }

  vec_dist <- Vectorize(dist, c("i", "j"))

  exp_corr_mat <- outer(seq_len(n_pixel), seq_len(n_pixel), function(i, j) {
    dist_ij <- vec_dist(i, j, n_col)
    exp(-dist_ij)
  })

  # generate errors (2000, 256)
  epsilon <- MASS::mvrnorm(
    n = na + nb, mu = rep(0, n_pixel), Sigma = exp_corr_mat
  )

  # generate y
  y <- outer(group_ind, beta) + epsilon

  return(cbind(group_ind, y))
}

gen_data_objs <- c("generate_data")
gen_data_pkgs <- c("MASS")


# 2 Visualization --------------------------------------------------------
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

# Decide the strength of center effect by visualization
set.seed(121)
df1 <- generate_data(center_effect = 4)
plot_matrix(df1[78, -1], range(df1))

df2 <- generate_data(center_effect = 5)
plot_matrix(df2[78, -1], range(df2))

image_4 <- plot_matrix(df1[78, -1], range(df1, df2))
image_5c <- plot_matrix(df2[1001, -1], range(df1, df2)) # w/o center effect
image_5 <- plot_matrix(df2[78, -1], range(df1, df2)) # w/ center effect

image_path <- file.path(getwd(), "Figures")

ggsave(file.path(image_path, "ex_image_4.png"), plot = image_4)
ggsave(file.path(image_path, "ex_image_5.png"), plot = image_5)
ggsave(file.path(image_path, "ex_image_5c.png"), plot = image_5c)

# Set up global center effect
center_effect <- 5
gen_data_objs <- c(gen_data_objs, "center_effect")


# 3. VBM --------------------------------------------------------------

n_iter <- 1000

# This function is designed exclusively for parallel execution.
# In each iteration, the function generates 2000 images. It then utilizes a
# grouping indicator to estimate the pixel values for each image. Each image
# consists of a total of 256 pixels.
# The output is a matrix with dimensions n_iters (number of iterations) by
# n_pixels (number of pixels), containing the p-values for each pixel.

vbm_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(center_effect = center_effect)

  n_pixel <- ncol(simulated_data) - 1
  pvals <- rep(NA, n_pixel)

  for (j in seq_len(n_pixel)) {
    pixel_value <- simulated_data[, j + 1]
    group_ind <- simulated_data[, 1]
    model <- lm(pixel_value ~ group_ind)
    pvals[j] <- summary(model)$coefficients[2, 4]
  }

  return(pvals)
}

vbm_pkgs <- c()
vbm_objs <- c()
list_package <- c(gen_data_pkgs, vbm_pkgs)

set.seed(42)
tic()
vbm_pvals <- simWrapper(
  n_sim = n_iter,
  f_sim = vbm_fsim,
  list_export = c(gen_data_objs, vbm_objs, "list_package"),
  list_package = list_package
)
toc()

vbm_pvals_corr <- t(apply(vbm_pvals, 1, p.adjust, method = "bonferroni"))

# calculate the perc of p-values < 0.05 for each pixel
vbm_pvals_perc <- colSums(vbm_pvals < 0.05) / n_iter * 100
vbm_pvals_corr_perc <- colSums(vbm_pvals_corr < 0.05) / n_iter * 100

image_vbm_pvals <- plot_matrix(vbm_pvals_perc, c(0, 100))
image_vbm_pvals_corr <- plot_matrix(vbm_pvals_corr_perc, c(0, 100))

image_path <- file.path(getwd(), "Figures")

ggsave(file.path(image_path, "vbm_pvals.png"), plot = image_vbm_pvals)
ggsave(file.path(image_path, "vbm_pvals_corr.png"), plot = image_vbm_pvals_corr)

# Build


# spVBM ------------------------------------------------------------
# In this simulation, there is no subject-level non-spatial random effects
# because each subject has only a single slice.

# get the eigenvectors, space=1 means no approximation for the spacing
# between the coordinates
# eigen_vecs <- mymeigen2D(
#   coords = df_long[, c("x", "y")], id = df_long$image_id, space = 1
# )
# tic()
# spvbm_fit <- myresf_vc(
#   y = df_long[["pixel_value"]],
#   x = df_long[, "group_ind"],
#   xgroup = factor(df_long$image_id),
#   meig = eigen_vecs
# )
# toc()
# # saveRDS(spvbm_fit, file = "spvbm_fit.rds")
# # readRDS(file = "spvbm_fit.rds")

# spvbm_coefs <- spvbm_fit$b_vc
# spvbm_coefs_2d <- matrix(rowSums(spvbm_coefs), 16, 16)
# image(
#   t(spvbm_coefs_2d),
#   col = gray(seq(1, 0, length = 256)),
#   main = "2D Visualization of Estimated Parameters"
# )
# plot.new()


# LASSO -------------------------------------------------------------
# For each iteration,
# group_indicator as the outcome;
# y (2000, 256) as the predictors;

train_test_split <- function(x, y, p_train) {
  set.seed(123)
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]

  return(list(
    x_train = x_train,
    x_test = x_test,
    y_train = y_train,
    y_test = y_test
  ))
}

# select optimal lambda via cross-validation
# evaluate performance via accuracy and AUC
perform_lasso <- function(x, y, p_train) {
  # Split the dataset
  split <- train_test_split(x, y, p_train)

  # Extract training and testing data
  x_tr <- split$x_train
  x_te <- split$x_test
  y_tr <- split$y_train
  y_te <- split$y_test

  # Perform cross-validation to find the optimal lambda values
  cv_mod <- cv.glmnet(x_tr, y_tr, alpha = 1, family = "binomial")
  l_min <- cv_mod$lambda.min
  l_1se <- cv_mod$lambda.1se

  # Function to evaluate model performance
  eval_perf <- function(l) {
    mod <- glmnet(x_tr, y_tr, alpha = 1, lambda = l, family = "binomial")
    preds <- predict(mod, newx = x_te, type = "response")[, 1]
    preds_bin <- ifelse(preds > 0.5, 1, 0)

    acc <- mean(preds_bin == y_te)
    auc <- pROC::auc(pROC::roc(y_te, preds))

    c(acc = acc, AUC = auc)
  }

  # Get performance metrics for both lambda.min and lambda.1se
  perf_min <- eval_perf(l_min)
  perf_1se <- eval_perf(l_1se)

  # Combine results into a coherent matrix
  results <- c(perf_min, perf_1se)
  results <- matrix(results, nrow = 1)
  colnames(results) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")

  return(results)
}

rm(df1, df2)
set.seed(121)
df <- generate_data(center_effect = 5)
x <- df[, -1]
y <- df[, 1]
perform_lasso(x, y, p_train = 0.8) # use all pixels cause perfect separation

# search how many pixels in the LASSO model will cause perfect separation

# indices for center pixels
c_indices <- {
  c_rows <- 5:12
  as.vector(outer(c_rows, c_rows, FUN = function(i, j) (i - 1) * 16 + j))
}
# indices for edge pixels
e_indices <- {
  all_indices <- 1:256
  setdiff(all_indices, c_indices)
}

i <- 1 # number of pixels from each side
m <- 0 # maximum AUC/accuracy
res <- matrix(NA, nrow = 0, ncol = 4)
while (m < 1) {
  ia <- sample(c_indices, i) + 1
  ib <- sample(e_indices, i) + 1
  x <- df[, c(ia, ib)]
  y <- df[, 1]
  res <- rbind(res, perform_lasso(x, y, p_train = 0.8))
  m <- max(res[i, ])
  i <- i + 1
}
# when use two points from each side, AUC/accuracy equals 1.
rm(x, y, ia, ib, i, m, res)

# Estimate p-value for each pixel using permutation test
perm_lasso <- function(x, y, n_perm) {
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

lasso_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(center_effect = 5)
  x <- simulated_data[, -1]
  y <- simulated_data[, 1]

  # fit lasso
  pvals <- perm_lasso(x, y, n_perm = 100)

  return(pvals)
}

lasso_pkgs <- c("glmnet", "pROC")
lasso_objs <- c("perm_lasso")
list_package <- c(gen_data_pkgs, lasso_pkgs)

set.seed(42)
tic()
lasso_pvals <- simWrapper(
  n_sim = 100,
  f_sim = lasso_fsim,
  list_package = list_package,
  list_export = c(gen_data_objs, lasso_objs, "list_package")
)
toc()

lasso_pvals_corr <- t(apply(lasso_pvals, 1, p.adjust, method = "bonferroni"))

# calculate the perc of p-values < 0.05 for each pixel
lasso_pvals_perc <- colSums(lasso_pvals < 0.05) / nrow(lasso_pvals) * 100
lasso_pvals_corr_perc <- colSums(lasso_pvals_corr < 0.05) / nrow(lasso_pvals) * 100

image_lasso_pvals <- plot_matrix(lasso_pvals_perc, c(0, 100))
image_lasso_pvals_corr <- plot_matrix(lasso_pvals_corr_perc, c(0, 100))

image_path <- file.path(getwd(), "Figures")

ggsave(file.path(image_path, "lasso_pvals.png"), plot = image_lasso_pvals)
ggsave(file.path(image_path, "lasso_pvals_corr.png"), plot = image_lasso_pvals_corr)


# Frequency --------------------------------------------------------
# Predict group_ind using image projected on the frequency domain
# Find a fix correlation matrix, such as exp corr mat
# Do the transformation use
#   1. positive eigenvalues only
#   2. all eigenvalues

exp_corr_mat <- function(n) {
  dist_mat <- outer(seq_len(n), seq_len(n), function(x, y) abs(x - y))
  corr_mat <- exp(-dist_mat / max(dist_mat))
  return(corr_mat)
}

eig_decomp <- function(C) {
  p <- ncol(C)
  M <- diag(p) - matrix(1, p, p) / p
  eig_data <- eigen(M %*% C %*% M, symmetric = TRUE)

  order_idx <- order(eig_data$values, decreasing = TRUE)
  eig_vecs <- eig_data$vectors[, order_idx]
  eig_vals <- eig_data$values[order_idx]

  return(list(
    eigenvectors = eig_vecs,
    eigenvalues = eig_vals
  ))
}



# parallel data generation and freq model fitting above
n_iter <- 100
n_perm <- 1000

freq_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(
    na = na,
    nb = nb,
    n_pixel = n_pixel,
    center_size = center_size
  )
  x <- simulated_data[, -1]
  y <- simulated_data[, 1]

  # build exponential correlation matrix
  corr_mat <- exp_corr_mat(ncol(x))

  # eigendecompose
  eig_comp <- eig_decomp(corr_mat)

  # eigen transpose in two ways
  pos_eigenvalues <- eig_comp$eigenvalues > 0
  x_trans_pos <- x %*% eig_comp$eigenvectors[, pos_eigenvalues]
  x_trans <- x %*% eig_comp$eigenvectors

  # Perform lasso regression on the transformed datasets
  evals_pos <- perform_lasso(x_trans_pos[, 1:2], y, p_train)
  evals <- perform_lasso(x_trans, y, p_train)

  # Combine results, assume evals and evals_pos are matrices or data frames
  results <- cbind(evals_pos, evals)

  return(results)
}

freq_objs <- c("exp_corr_mat", "eig_decomp", "perm_lasso", "n_perm")
freq_pkgs <- c("glmnet")
list_package <- c(gen_data_pkgs, freq_pkgs)

set.seed(42)
tic()
freq_pvals <- call_simWrapper(
  n_sim = n_iter,
  f_sim = freq_fsim,
  list_export = c(gen_data_objs, freq_objs, "list_package"),
  list_package = list_package
)
toc()




generate_data <- function(beta_eff) {
  na + nb <- 2000
  square_size <- sqrt(256)

  # generate group indicator
  group_ind <- c(rep(1, 1000), rep(0, 1000)) # (2000, )

  # group effect by pixel
  beta <- matrix(0, square_size, square_size)
  index_st <- (square_size - center_size) %/% 2 + 1
  index_end <- index_st + center_size - 1
  beta[index_st:index_end, index_st:index_end] <- beta_eff
  beta <- as.vector(beta) # (256, )

  # multivariate normal error
  exp_corr_mat <- function(n) {
    dist_mat <- outer(seq_len(n), seq_len(n), function(x, y) abs(x - y))
    corr_mat <- exp(-dist_mat / max(dist_mat))
    return(corr_mat)
  }
  corr_mat <- exp_corr_mat(n_pixel)

  # generate errors (2000, 256)
  epsilon <- mvrnorm(n = na + nb, mu = rep(0, n_pixel), Sigma = corr_mat)

  # generate y
  y <- outer(group_ind, beta) + epsilon

  return(cbind(group_ind, y))
}

beta_eff <- 1.5
df_sample <- generate_data(beta_eff)
plot_matrix(df_sample[78, -1], range(df_sample))
