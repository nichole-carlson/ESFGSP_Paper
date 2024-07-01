library(MASS)
library(glmnet)
library(pROC)
library(tictoc)

parent_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/"
source(file.path(parent_dir, "Dissemination", "simulation.r"))
fig_dir <- file.path(parent_dir, "Figures")


# Define the covariance matrix generation function
gen_cov_matrix <- function(img_size) {
  coords <- expand.grid(1:img_size, 1:img_size)
  dist_matrix <- as.matrix(dist(coords))
  w <- exp(-dist_matrix)
  return(w)
}

# Define the beta vector generation function
define_beta <- function(img_size, beta_effect) {
  p <- img_size^2
  beta <- rep(0, p)
  center_idx <- as.vector(matrix(1:p, nrow = img_size)[5:12, 5:12])
  beta[center_idx] <- beta_effect
  return(beta)
}

# Define the X matrix generation function
gen_x <- function(num_samples, img_size, w) {
  p <- img_size^2
  x <- mvrnorm(num_samples, mu = rep(0, p), Sigma = w)
  return(x)
}

# Define the response variable generation function
gen_response <- function(x, beta) {
  eta <- x %*% beta
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(nrow(x), 1, p)
  return(y)
}


# compare_beta_effects <- function(beta_effects = seq(0.1, 1.0, by = 0.1)) {
#   par(mfrow = c(ceiling(length(beta_effects) / 2), 2))
#   for (beta_effect in beta_effects) {
#     sim1_obj <- sim1$new(beta_effect = beta_effect)
#     p <- sim1_obj$check_p()
#     title <- paste("beta_effect =", beta_effect)
#     hist(p,
#       breaks = 30, main = title, , xlim = c(0, 1),
#       xlab = "Probability p", col = "lightblue", border = "black"
#     )
#   }
# }

# # Run the function to compare different beta effects
# png(
#   file = file.path(fig_dir, "sim1_p_dist.png"),
#   width = 1600, height = 1200, res = 150
# )
# compare_beta_effects(beta_effects = c(1, 0.2, 0.1, 0.05, 0.01))
# dev.off()


# Run simulation 1 for 100 times.
# For each iteration:
#   - Split into train and test dataset, calculate AUC and accuracy, estimated
#     coefficients under lambda.min and lambda.1se
#   - Perform permutation test. Estimate p-values for pixels


# Main simulation function
simulate_lasso <- function(i, img_size, num_samples, beta_effect, p_train, n_perm, seed) {
  perf_metrics <- matrix(NA, nrow = 1, ncol = 4)
  perm_pvals <- matrix(NA, nrow = 1, ncol = img_size^2)

  set.seed(seed + i)
  w <- gen_cov_matrix(img_size)
  beta <- define_beta(img_size, beta_effect)
  x <- gen_x(num_samples, img_size, w)
  y <- gen_response(x, beta)
  lasso_results <- perform_lasso(x, y, p_train, seed = seed + i)
  perf_metrics[1, ] <- lasso_results

  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  perm_pvals[1, ] <- p_vals

  return(cbind(perf_metrics, perm_pvals))
}


tic()
sim_output_1 <- simWrapper(
  n_sim = 100,
  f_sim = function(i) simulate_lasso(i, img_size = 16, num_samples = 1000, beta_effect = 0.1, p_train = 0.8, n_perm = 100, seed = 42),
  list_export = c("gen_cov_matrix", "define_beta", "gen_x", "gen_response", "simulate_lasso", "perform_lasso", "perm_lasso"),
  list_package = c("MASS", "glmnet", "pROC", "foreach", "doParallel")
)
toc()

save(
  sim_output_1,
  file = file.path(dirname(parent_dir), "Simulations", "sim_1_output.RData")
)
