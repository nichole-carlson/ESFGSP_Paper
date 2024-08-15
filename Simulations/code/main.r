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
