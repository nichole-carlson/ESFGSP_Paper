library(MASS)
library(stats)
library(tidyverse)
library(hdi)
library(tictoc)
library(viridis)
library(reshape2)

source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")


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
n_image <- n_a + n_b
n_pixel <- 256
n_iter <- 100
image_size <- sqrt(n_pixel)
center_size <- 8

# generate z, group indicator
group_ind <- c(rep(1, n_a), rep(0, n_b))
group_ind_long <- rep(rep(group_ind, each = n_pixel), n_iter)

# group beta, effect by pixel
beta <- matrix(0, image_size, image_size)
index_st <- (image_size - center_size) %/% 2 + 1
index_end <- index_st + center_size - 1
beta[index_st:index_end, index_st:index_end] <- 1
beta <- as.vector(beta)
beta_long <- rep(beta, n_iter * n_image)

# multivariate normal error
exp_corr_mat <- function(n, rate) {
  dist_mat <- outer(seq_len(n), seq_len(n), function(x, y) abs(x - y))
  corr_mat <- exp(-rate * dist_mat)
  return(corr_mat)
}
corr_mat <- exp_corr_mat(n_pixel, rate = 1)

epsilon <- mvrnorm(n = n_iter * n_image, mu = rep(0, n_pixel), Sigma = corr_mat)

# y
y <- beta_long * group_ind_long + epsilon
y_3d <- array(NA, dim = c(n_iter, n_image, n_pixel))
s <- 1
for (i in seq_len(n_iter)) {
  for (j in seq_len(n_image)) {
    y_ij <- y[s:(s + n_pixel - 1)]
    y_3d[i, j, ] <- y_ij
    s <- s + n_pixel
  }
}
rm(s, i, j)


# Visualization --------------------------------------------------------
plot_matrix <- function(vec) {
  # Convert the matrix to a long format data frame
  mat <- matrix(vec, 16, 16, byrow = TRUE)
  data_long <- melt(mat)

  # Create the plot using ggplot2
  plot <- ggplot(data_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "Value") +
    theme_minimal() +
    labs(x = "", y = "")

  return(plot)
}

plot_matrix(y_3d[1, 62, ])
plot_matrix(y_3d[35, 78, ])


# VBM ----------------------------------------------------------------

vbm_pvals <- matrix(0, n_iter, n_pixel)

for (i in seq_len(n_iter)) {
  for (j in seq_len(n_pixel)) {
    y_pixel <- y_3d[i, , j]
    model <- lm(y_pixel ~ group_ind)
    vbm_pvals[i, j] <- summary(model)$coefficients[2, 4]
  }
  vbm_pvals[i, ] <- p.adjust(vbm_pvals[i, ], method = "holm")
}

# calculate the perc of p-values < 0.05 for each pixel
vbm_pvals_mat <- colSums(vbm_pvals < 0.05) / n_iter * 100 %>%
  matrix(., image_size, image_size)

ggplot(melt(vbm_pvals_mat), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  labs(title = "Corrected p-values") +
  theme_minimal()


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
# randomly divide the data into two parts, fit LASSO in one part, obtain pvals
# using OLS on the other part

lasso_pvals <- matrix(NA, 100, 256)
tic()
for (i in seq_len(n_iter)) {
  y_i <- y_3d[i, , ]
  model <- lasso.proj(
    y_i,
    group_ind,
    family = "binomial",
    parallel = TRUE,
    ncores = 4
  )
  lasso_pvals[i, ] <- model$pval.corr
}
toc()

lasso_pvals_mat <- colSums(lasso_pvals < 0.05) / n_iter * 100 %>%
  matrix(., image_size, image_size)

# no pixel is significant even before correction


# Frequency --------------------------------------------------------
# Predict z using y projected on the frequency domain
# For each iteration,
# 1. Calculate the correlation matrix C;
# 2. Calculate the eigenvalues and eigenvectors;
# 3. Project the pixel values (y) by different eigenvectors

emp_corr <- function(mat) {
  c_emp <- cor(mat)
  diag(c_emp) <- 0
  return(c_emp)
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

weighted_contribs_mat <- matrix(NA, n_iter, n_pixel)

for (i in 1:n_iter) {
  y_i <- y_3d[i, , ]
  C <- emp_corr(y_i)
  eig_comp <- eig_decomp(C)

  # Subset the eigenvectors with positive eigenvalues
  pos_eig_vecs <- eig_comp$eigenvectors[, eig_comp$eigenvalues > 0]

  # Transform y using the positive eigenvectors
  y_trans <- y_i %*% pos_eig_vecs

  # Fit LASSO model
  model <- lasso.proj(
    y_trans,
    group_ind,
    family = "binomial",
    parallel = TRUE,
    ncores = 4
  )

  # Find significant predictors
  pvals <- model$pval.corr
  sig <- (pvals < 0.05)
  if (sum(sig) == 0) {
    weighted_contribs_mat[i, ] <- rep(0, n_pixel)
  } else {
    sig_eig_vecs <- pos_eig_vecs[, sig, drop = FALSE]
    contribs <- abs(sig_eig_vecs)
    sig_vals <- -log(pvals[sig])
    weighted_contribs_mat[i, ] <- contribs %*% matrix(sig_vals, ncol = 1)[, 1]
  }
}

weighted_contribs %>%
  matrix(., image_size, image_size) %>%
  melt() %>%
  ggplot(., aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  labs(title = "Heatmap of Variable Significance") +
  theme_minimal()
