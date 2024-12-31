sys_info <- Sys.info()
if (!is.null(sys_info) && sys_info["sysname"] == "Linux") {
  proj_dir <- "/projects/sren@xsede.org/ESFGSP_Paper"
  data_dir <- "/scratch/alpine/sren@xsede.org/esfgsp"
} else {
  proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"
  data_dir <- file.path(proj_dir, "results")
}

# Source functions for simulating data
source(file.path(proj_dir, "R", "simulate_data.R"))

# Global parameters for simulations
n_iter <- 500
n_samples <- 1000
n_row <- 16
n_col <- 16
seed <- 42
sparse_level <- 0.1
beta_effect <- 0.1
b_effect <- 0.2


# Simulation 1a. Exp decay adjacency matrix, in pixel space.
sim1a_exp_pixel_data <- list()

cat("Running simulation 1 in pixel space ...\n")
tic()
for (i in seq_len(n_iter)) {
  sim1a_exp_pixel_data[[i]] <- run_pixel_to_freq_simulation(
    n_samples = n_samples,
    n_row = n_row,
    n_col = n_col,
    effect_size = beta_effect,
    c_adj = generate_exp_corr_matrix(n_row, n_col),
    seed = seed
  )
}
cat("Simulation completed. \n")
toc()

cat("Saving simulated data as a .RData file ...\n")
tic()
save(
  sim1a_exp_pixel_data,
  file = file.path(
    data_dir,
    paste0("sim1a_pixel_data_", format(Sys.Date(), "%y%m%d"), ".RData")
  )
)
toc()


# # Simulation 1b. Exp decay adjacency matrix, on freq space.
# sim1b_exp_freq_data <- list()
#
# for (i in seq_len(n_iter)) {
#   sim1b_exp_freq_data[[i]] <- run_simulation_1b(
#     n_samples = n_samples,
#     n_row = n_row,
#     n_col = n_col,
#     sparse_level = sparse_level,
#     effect_size = b_effect,
#     seed = seed
#   )
# }
#
# save(
#   sim1b_exp_freq_data,
#   file = file.path(
#     proj_dir, "results",
#     paste0("sim1b_freq_data_", format(Sys.Date(), "%y%m%d"), ".RData")
#   )
# )
#
#
# # Simulation 2a. 2-neighbor adjacency matrix, on pixel space.
# sim2a_2n_pixel_data <- list()
#
# for (i in seq_len(n_iter)) {
#   sim2a_2n_pixel_data[[i]] <- run_pixel_to_freq_simulation(
#     n_samples = n_samples,
#     n_row = n_row,
#     n_col = n_col,
#     effect_size = beta_effect,
#     c_adj = generate_n_neighbor_matrix(n_row, n_col, d_max = 2),
#     seed = seed
#   )
# }
