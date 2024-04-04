library(MASS)
library(stats)
library(tidyverse)
library(tictoc)

source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")


# Simulate data --------------------------------------------------------
n_a <- 1000 # of images in group A
n_b <- 1000
n <- n_a + n_b
n_pixels <- 256
image_size <- sqrt(n_pixels)
center_size <- 8

group_ind <- c(rep(1, n_pixels * n_a), rep(0, n_pixels * n_b)) # (n*n_pixels, 1)

beta <- matrix(0, image_size, image_size)
center_start <- floor((image_size - center_size) / 2) + 1
center_end <- center_start + center_size - 1
beta[center_start:center_end, center_start:center_end] <- 1
beta <- as.vector(beta) # (n_pixels, 1)

exp_corr_mat <- function(n, rate) {
  indices <- 0:(n - 1)
  dist_mat <- abs(outer(indices, indices, "-"))
  corr_mat <- exp(-rate * dist_mat)
  return(corr_mat)
}

corr_mat <- exp_corr_mat(n_pixels, rate = 1)

# epsilon is (n*n_pixels, 1)
epsilon <- as.vector(mvrnorm(n, mu = rep(0, n_pixels), Sigma = corr_mat))

pixel_value <- beta * group_ind + epsilon

df_long <- {
  coord_grid <- expand.grid(x = 1:image_size, y = 1:image_size)
  data.frame(
    image_id = rep(1:n, each = n_pixels),
    pixel_id = rep(1:n_pixels, n),
    x = rep(coord_grid$x, n),
    y = rep(coord_grid$y, n),
    pixel_value = pixel_value,
    group_ind = group_ind
  )
}


# VBM ----------------------------------------------------------------
results <- list()
for (pixel in seq_len(ncol(pixel_value))) {
  y_pixel <- pixel_value[, pixel]
  model <- lm(y_pixel ~ z - 1)
  results[[pixel]] <- summary(model)
}
vbm_pvals <- sapply(results, function(res) coef(res)[1, 4])
vbm_coefs <- sapply(results, function(res) coef(res)[1, 1])

vbm_coefs_2d <- matrix(vbm_coefs, nrow = 16, ncol = 16)
image(
  t(vbm_coefs_2d),
  col = gray.colors(256),
  main = "2D Visualization of Estimated Parameters"
)
plot.new()

vbm_adjust_pvals <- p.adjust(vbm_pvals, method = "bonferroni")
vbm_adjust_pvals_2d <- matrix(vbm_adjust_pvals, nrow = 16, ncol = 16)
image(
  t(vbm_adjust_pvals_2d),
  col = gray.colors(256), main = "Corrected p-values"
)


# spVBM ------------------------------------------------------------
# In this simulation, there is no subject-level non-spatial random effects
# because each subject has only a single slice.

# get the eigenvectors, space=1 means no approximation for the spacing
# between the coordinates
eigen_vecs <- mymeigen2D(
  coords = df_long[, c("x", "y")], id = df_long$image_id, space = 1
)
tic()
spvbm_fit <- myresf_vc(
  y = df_long[["pixel_value"]],
  x = df_long[, "group_ind"],
  xgroup = factor(df_long$image_id),
  meig = eigen_vecs
)
toc()
# saveRDS(spvbm_fit, file = "spvbm_fit.rds")
# readRDS(file = "spvbm_fit.rds")

spvbm_coefs <- spvbm_fit$b_vc
spvbm_coefs_2d <- matrix(rowSums(spvbm_coefs), 16, 16)
image(
  t(spvbm_coefs_2d),
  col = gray(seq(1, 0, length = 256)),
  main = "2D Visualization of Estimated Parameters"
)
plot.new()
