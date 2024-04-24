setup_packages <- function(packages) {
  for (package in packages) {
    # Check if the package is installed
    if (!require(package, character.only = TRUE)) {
      cat(sprintf("Installing package: %s\n", package))
      install.packages(package)

      # Attempt to load the package again after installing
      if (!require(package, character.only = TRUE)) {
        cat(sprintf("Failed to install package: %s\n", package))
      } else {
        cat(sprintf("Package loaded successfully: %s\n", package))
      }
    } else {
      cat(sprintf("Package already installed and loaded: %s\n", package))
    }
  }
}

# List of packages to check and install if necessary
packages_to_install <- c(
  "MASS", "stats", "tidyverse", "hdi", "tictoc", "viridis", "reshape2",
  "parallel", "glmnet", "pROC"
)

setup_packages(packages_to_install)

source("/Users/siyangren/Documents/ESFGSP/simWrapper.r")
# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
# source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")



# Parallel global settings --------------------------------------------
TF_parallel <- TRUE
n_cores <- detectCores() - 1

call_simWrapper <- function(
    n_sim, f_sim, list_package, list_export, params = list()) {
  # Establish default parameters for the simWrapper function
  default_params <- list(
    n_sim = n_sim,
    f_sim = f_sim,
    TF_parallel = TF_parallel,
    n_cores = n_cores,
    list_package = list_package,
    list_export = list_export
  )

  # Merge user-specified parameters with defaults
  # Any additional parameters in 'params' will override the defaults if provided
  args <- modifyList(default_params, params)

  # Use do.call to execute simWrapper with these arguments
  do.call("simWrapper", args)
}


# Simulate data ------------------------------------------------------
# 1. Generate group indicator, 1000 in group A; 1000 in group B
# 2. Generate the beta as 16*16 matrix, it has 1 in the center 8*8 and 0 in
#    other areas.
# 3. Generate the exponential correlation matrix
# 4. Use the corr matrix to generate the epsilon as a multivariate normal.
#    Epsilon should be (100, 2000, 256)
# 5. Broadcast and generate y

n_a <- 1000
n_b <- 1000
n_pixel <- 256
center_size <- 8

generate_data <- function(n_a, n_b, n_pixel, center_size) {
  n_image <- n_a + n_b
  square_size <- sqrt(n_pixel)

  # generate group indicator
  group_ind <- c(rep(1, n_a), rep(0, n_b)) # (2000, )

  # group effect by pixel
  beta <- matrix(0, square_size, square_size)
  index_st <- (square_size - center_size) %/% 2 + 1
  index_end <- index_st + center_size - 1
  beta[index_st:index_end, index_st:index_end] <- 8
  beta <- as.vector(beta) # (256, )

  # multivariate normal error
  exp_corr_mat <- function(n) {
    dist_mat <- outer(seq_len(n), seq_len(n), function(x, y) abs(x - y))
    corr_mat <- exp(-dist_mat / max(dist_mat))
    return(corr_mat)
  }
  corr_mat <- exp_corr_mat(n_pixel)

  # generate errors (2000, 256)
  epsilon <- mvrnorm(n = n_image, mu = rep(0, n_pixel), Sigma = corr_mat)

  # generate y
  y <- outer(group_ind, beta) + epsilon

  return(cbind(group_ind, y))
}

gen_data_objs <- c(
  "n_a", "n_b", "n_pixel", "center_size", "generate_data"
)
gen_data_pkgs <- c("MASS")


# Visualization --------------------------------------------------------
set.seed(121)
df_sample <- generate_data(n_a, n_b, n_pixel, center_size)

plot_matrix <- function(vec, value_limits = c()) {
  # Convert the matrix to a long format data frame
  mat <- matrix(vec, 16, 16, byrow = TRUE)

  # Convert the matrix to a long format data frame
  data_long <- melt(mat)

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

image_a <- plot_matrix(df_sample[78, -1], range(df_sample))
image_b <- plot_matrix(df_sample[1024, -1], range(df_sample))

image_path <- file.path(getwd(), "Figures")

ggsave(file.path(image_path, "image_ex1.png"), plot = image_a)
ggsave(file.path(image_path, "image_ex2.png"), plot = image_b)



# VBM ----------------------------------------------------------------

n_iter <- 1000

vbm_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(
    n_a = n_a,
    n_b = n_b,
    n_pixel = n_pixel,
    center_size = center_size
  )

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

tic()
vbm_pvals <- call_simWrapper(
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

# n_iter <- 1000

# lasso_fsim <- function(i) {
#   # simulate data
#   simulated_data <- generate_data(
#     n_a = n_a,
#     n_b = n_b,
#     n_pixel = n_pixel,
#     center_size = center_size
#   )
#   x <- simulated_data[, -1]
#   y <- simulated_data[, 1]

#   # fit lasso model
#   model <- lasso.proj(x[, c(1, 3, 123, 124)], y, family = "binomial")

#   return(model$pval.corr)

#   # model <- cv.glmnet(x[, c(1, 3, 123, 124)], y, family = "binomial", alpha = 1)
#   # best_lambda <- model$lambda.min
#   # coefs <- coef(model, s = best_lambda)[-1, 1]

#   # return(coefs)
# }

# lasso_pkgs <- c("hdi", "glmnet")
# lasso_objs <- c()
# list_package <- c(gen_data_pkgs, lasso_pkgs)

# set.seed(42)
# tic()
# lasso_coefs <- call_simWrapper(
#   n_sim = n_iter,
#   f_sim = lasso_fsim,
#   list_export = c(gen_data_objs, lasso_objs, "list_package"),
#   list_package = list_package
# )
# toc()

# lasso_pvals_mat <- colSums(lasso_pvals < 0.05) / n_iter * 100

# no pixel is significant even before correction

n_iter <- 1000
p_train <- 0.8

perform_lasso <- function(x, y, p_train) {
  # Split data into training and test sets
  set.seed(123)
  n <- nrow(x)
  train_idx <- sample(seq_len(n), size = floor(n * p_train))

  x_train <- as.matrix(x[train_idx, ])
  y_train <- y[train_idx]
  x_test <- as.matrix(x[-train_idx, ])
  y_test <- y[-train_idx]

  # Cross-validation to find the optimal lambda
  cv_model <- cv.glmnet(x_train, y_train, alpha = 1)
  lambda_best <- cv_model$lambda.min
  lambda_1se <- cv_model$lambda.1se

  # Fit Lasso models on training set with both lambda.min and lambda.1se
  model_min <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_best)
  model_1se <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_1se)

  # Evaluate model performance on the test set
  preds_min <- predict(model_min, s = lambda_best, newx = x_test)
  preds_1se <- predict(model_1se, s = lambda_1se, newx = x_test)

  # Calculate AUC
  auc_min <- pROC::auc(pROC::roc(y_test, preds_min[, 1]))
  auc_1se <- pROC::auc(pROC::roc(y_test, preds_1se[, 1]))

  res <- matrix(c(auc_min, auc_1se), nrow = 1)
  return(res)
}

lasso_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(
    n_a = n_a,
    n_b = n_b,
    n_pixel = n_pixel,
    center_size = center_size
  )
  x <- simulated_data[, -1][, c(1, 122)]
  y <- simulated_data[, 1]

  # fit lasso
  aucs <- perform_lasso(x, y, p_train)

  return(aucs)
}

lasso_pkgs <- c("glmnet", "pROC")
lasso_objs <- c("p_train", "perform_lasso")
list_package <- c(gen_data_pkgs, lasso_pkgs)

set.seed(42)
tic()
lasso_aucs <- call_simWrapper(
  n_sim = n_iter,
  f_sim = lasso_fsim,
  list_export = c(gen_data_objs, lasso_objs, "list_package"),
  list_package = list_package
)
toc()


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

perm_lasso <- function(x, y, n_perm) {
  # get original coef estimates
  cv_fit <- cv.glmnet(x, y, alpha = 1)
  best_lambda <- cv_fit$lambda.min
  orig_coef <- coef(cv_fit, s = "lambda.min")[-1, 1]

  perm_coefs <- matrix(NA, n_perm, length(orig_coef))

  # perform permutations
  for (i in seq_len(n_perm)) {
    y_perm <- sample(y)
    perm_fit <- glmnet(x, y_perm, alpha = 1)
    perm_coefs[i, ] <- coef(perm_fit, s = best_lambda)[-1, 1]
  }

  # calculate p-values
  p_vals <- matrix(colMeans(abs(perm_coefs) >= abs(orig_coef)), nrow = 1)
  return(p_vals)
}



# parallel data generation and freq model fitting above
n_iter <- 100
n_perm <- 1000

freq_fsim <- function(i) {
  # simulate data
  simulated_data <- generate_data(
    n_a = n_a,
    n_b = n_b,
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
  x_trans_pos <- x %*% eig_comp$eigenvectors[, eig_comp$eigenvalues > 0]
  x_trans <- x %*% eig_comp$eigenvectors

  # perform lasso with permutation test
  p_vals <- perm_lasso(x_trans, y, n_perm)
  p_vals_pos <- perm_lasso(x_trans_pos, y, n_perm)

  return(cbind(p_vals, p_vals_pos))
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

saveRDS(freq_pvals, file = "freq_pvals.rds")

# p-values when using all eigenvectors for transformation
freq_pvals_total <- freq_pvals[, seq_len(n_pixel)]
# p-values when using eigenvectors with postive eigenvalues only
freq_pvals_pos <- freq_pvals[, (n_pixel + 1):ncol(freq_pvals)]
