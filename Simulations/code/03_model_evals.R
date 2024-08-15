library(table1)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
results_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load model fitting results -----
load(file = file.path(results_data_dir, "model_metrics_240807.RData"))

my.render.cont <- function(x) {
  with(stats.default(x), sprintf("%.3f (%0.3f)", MEAN, SD))
}


table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim1_model_perf$pixel,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim1_model_perf$freq,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim2_model_perf$pixel,
  render.continuous = my.render.cont
)
table1(
  ~ min_acc + min_auc + `1se_acc` + `1se_auc`,
  data = sim2_model_perf$freq,
  render.continuous = my.render.cont
)
