# ----- Create table for AUC and accuracy -----
my.render.continuous <- function(x) {
  with(stats.default(x), sprintf("%.3f (%0.3f)", MEAN, SD))
}

table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim1_auc_acc$pixel),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim1_auc_acc$freq),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim2_auc_acc$pixel),
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = as.data.frame(sim2_auc_acc$freq),
  render.continuous = my.render.cont
)


# ----- Visualize significant p-values -----
# percentage of significant p-values
p_pvals <- list()
for (i in seq_len(length(sim_pvals_list))) {
  dat <- sim_pvals_list[[i]]
  p_pvals[[paste0("sim", i)]]$beta <- plot_heatmap(colMeans(dat$pixel < 0.05), c(0, 1))
  p_pvals[[paste0("sim", i)]]$b <- plot_scatterplot(eig_vals_ordered, colMeans(dat$freq < 0.05))
}

png(
  file = file.path(fig_dir, "perc_sign_pvals_beta.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals$sim1$beta, p_pvals$sim2$beta), ncol = 2)
dev.off()

png(
  file = file.path(fig_dir, "perc_sign_pvals_b.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals$sim1$b, p_pvals$sim2$b), ncol = 2)
dev.off()

p_pvals_hdi <- list()
p_pvals_hdi$sim1$beta <- plot_heatmap(colMeans(sim1_pvals_hdi$pixel < 0.05), c(0, 1))
p_pvals_hdi$sim1$b <- plot_scatterplot(eig_vals_ordered, colMeans(sim1_pvals_hdi$freq < 0.05))
p_pvals_hdi$sim2$beta <- plot_heatmap(colMeans(sim2_pvals_hdi$pixel < 0.05), c(0, 1))
p_pvals_hdi$sim2$b <- plot_scatterplot(eig_vals_ordered, colMeans(sim2_pvals_hdi$freq < 0.05))

png(
  file = file.path(fig_dir, "perc_sign_pvals_hdi_beta.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals_hdi$sim1$beta, p_pvals_hdi$sim2$beta), ncol = 2)
dev.off()

png(
  file = file.path(fig_dir, "perc_sign_pvals_hdi_b.png"),
  width = 1200, height = 500, res = 150
)
grid.arrange(grobs = list(p_pvals_hdi$sim1$b, p_pvals_hdi$sim2$b), ncol = 2)
dev.off()
