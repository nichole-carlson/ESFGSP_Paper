library(MASS)
library(stats)
library(tidyverse)

source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/resf_vc.R")
source("/Users/siyangren/Documents/ra-cida/spatial-filtering/code/mymeigen2D.R")


# Simulate data --------------------------------------------------------
n_a <- 1000
n_b <- 1000
s <- 256

z <- c(rep(1, n_a), rep(0, n_b))

beta <- matrix(0, 16, 16)
central_start <- (16 - 8) %/% 2
central_end <- central_start + 8
beta[(central_start + 1):(central_end), (central_start + 1):(central_end)] <- 1
beta <- as.vector(beta)

exp_corr_mat <- function(n, rate) {
  indices <- 0:(n - 1)
  dist_mat <- abs(outer(indices, indices, "-"))
  corr_mat <- exp(-rate * dist_mat)
  return(corr_mat)
}

corr_mat <- exp_corr_mat(s, rate = 1)

y <- matrix(0, nrow = n_a + n_b, ncol = s)
for (i in 1:(n_a + n_b)) {
  y[i, ] <- mvrnorm(1, mu = z[i] * beta, Sigma = corr_mat)
}


# VBM ----------------------------------------------------------------
results <- list()
for (pixel in seq_len(ncol(y))) {
  y_pixel <- y[, pixel]
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

spvbm_df <- {
  coords <- expand.grid(x = 1:16, y = 1:16)
  coords$pixel_id <- 1:256
  y_long <- reshape2::melt(y)
  colnames(y_long) <- c("image_id", "pixel_id", "pixel_value")
  y_long$z <- z[y_long$image_id]
  merge(y_long, coords, by = "pixel_id")
}
# space=1 means no approximation for the spacing between the coordinates
eigen_vecs <- mymeigen2D(
  coords = spvbm_df[, c("x", "y")], id = spvbm_df$image_id, space = 1
)
spvbm_fit <- myresf_vc(
  y = spvbm_df[["pixel_value"]],
  x = spvbm_df[, "z"],
  xgroup = factor(spvbm_df$image_id),
  meig = eigen_vecs
)
