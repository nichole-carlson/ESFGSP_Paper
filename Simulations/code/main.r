# Function for Simulation 1
simulate_1 <- function(i, size, n_samples, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  beta <- define_beta(size, effect)
  x <- generate_X(n_samples, size, w)
  y <- generate_response(x, beta)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}


# Run Simulation 1
# For each iteration:
#   - Split into train and test dataset, calculate AUC and accuracy, estimated
#     coefficients under lambda.min and lambda.1se
#   - Perform permutation test. Estimate p-values for pixels

tic()
sim1_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_1(i, size = 16, n_samples = 1000, effect = 0.1, p_train = 0.8, n_perm = 100, seed = 42),
  list_export = c(
    "generate_cov_matrix", "define_beta", "generate_X", "generate_probs",
    "generate_response", "simulate_1", "perform_lasso", "perm_lasso"
  ),
  list_package = c("MASS", "glmnet", "pROC", "foreach", "doParallel")
)
toc()






# # Save Simulation 1 results
# save(
#   sim1_output,
#   file = file.path(parent_dir, "Simulations", paste0("sim1_", format(Sys.time(), "%y%m%d"), ".RData"))
# )
load(file = file.path(current_dir, "Simulations", "sim1_240701.RData"))

# # Generate summary stats for ACC and AUCs
acc_auc_sim1 <- as.data.frame(sim1_output[, 1:4])
colnames(acc_auc_sim1) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = acc_auc_sim1)

# # Calculate significance percentages for Simulation 1
pvals_matrix_sim1 <- sim1_output[, -(1:4)]
sig_perc_sim1 <- calc_sig_perc(pvals_matrix_sim1, method = "bonferroni", alpha = 0.05)

# Create plots for Simulation 1
beta_sim1 <- define_beta(16, 0.1)
plot_beta_sim1 <- plot_heatmap(beta_sim1)
plot_orig_sim1 <- plot_heatmap(sig_perc_sim1$perc_orig, c(0, 100))
plot_adj_sim1 <- plot_heatmap(sig_perc_sim1$perc_adj, c(0, 100))

# Arrange plots for formal presentation for Simulation 1
sim1_heatmap <- grid.arrange(
  plot_beta_sim1, plot_orig_sim1, plot_adj_sim1,
  ncol = 2
)
ggsave(
  filename = file.path(current_dir, "Figures", "sim1_heatmap.png"),
  plot = sim1_heatmap, width = 10, height = 8
)

# Function for Simulation 2
simulate_2 <- function(i, size, n_samples, sparsity, effect, p_train, n_perm, seed) {
  set.seed(seed + i)
  w <- generate_cov_matrix(size)
  v <- eigen_decomp(w)$vectors
  x <- generate_X(n_samples, size, diag(size^2))
  b <- generate_sparse_vector(size^2, sparsity, effect, seed)
  y <- generate_response(x, b)
  perform_metrics <- perform_lasso(x, y, p_train, seed = seed + i)
  p_vals <- perm_lasso(x, y, n_perm, seed = seed + i)
  cbind(perform_metrics, p_vals)
}



# Run Simulation 2
tic()
sim2_output <- simWrapper(
  n_sim = 500,
  f_sim = function(i) simulate_2(i, size = 16, n_samples = 1000, sparsity = 0.1, effect = 0.4, p_train = 0.8, n_perm = 100, seed = 42),
  list_export = c(
    "generate_cov_matrix", "eigen_decomp", "generate_sparse_vector", "generate_X", "generate_probs", "generate_response", "perform_lasso",
    "perm_lasso", "simulate_2"
  ),
  list_package = c("MASS", "glmnet", "pROC", "foreach", "doParallel")
)
toc()

# Save Simulation 2 results
# save(
#   sim2_output,
#   file = file.path(parent_dir, "Simulations", paste0("sim2_", format(Sys.time(), "%y%m%d"), ".RData"))
# )
load(file = file.path(current_dir, "Simulations", "sim2_240701.RData"))

# Generate summary stats for ACC and AUCs for Simulation 2
acc_auc_sim2 <- as.data.frame(sim2_output[, 1:4])
colnames(acc_auc_sim2) <- c("min_acc", "min_auc", "1se_acc", "1se_auc")
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = acc_auc_sim2)

# Calculate significance percentages for Simulation 2
pvals_matrix_sim2 <- sim2_output[, -(1:4)]
sig_perc_sim2 <- calc_sig_perc(pvals_matrix_sim2, method = "bonferroni", alpha = 0.05)


# Create plots for Simulation 2
b_sim2 <- generate_sparse_vector(16^2, 0.1, 0.4, 42)
plot_beta_sim2 <- plot_heatmap(b_sim2)
plot_orig_sim2 <- plot_heatmap(sig_perc_sim2$perc_orig, c(0, 100))
plot_adj_sim2 <- plot_heatmap(sig_perc_sim2$perc_adj, c(0, 100))

# Arrange plots for formal presentation for Simulation 2
sim2_heatmap <- grid.arrange(
  plot_beta_sim2, plot_orig_sim2, plot_adj_sim2,
  ncol = 2
)
ggsave(
  filename = file.path(current_dir, "Figures", "sim2_heatmap.png"),
  plot = sim2_heatmap, width = 10, height = 8
)
