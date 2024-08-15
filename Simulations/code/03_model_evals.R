library(table1)

simulation_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/Simulations"
code_dir <- file.path(simulation_dir, "code")
results_data_dir <- file.path(simulation_dir, "results", "data")


# ----- Load model fitting results -----
load(file = file.path(results_data_dir, "model_metrics_240807.RData"))

table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = sim1_model_perf$pixel)
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = sim1_model_perf$freq)
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = sim2_model_perf$pixel)
table1(~ min_acc + min_auc + `1se_acc` + `1se_auc`, data = sim2_model_perf$freq)
