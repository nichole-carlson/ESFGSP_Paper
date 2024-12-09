library(tictoc)

proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"
p_data_dir <- "/Volumns/alzheimersdisease/ESFGSPproject/DataLibrary"

# Source the functions for simulating data
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


# ---------- Run Simulation 1 ----------
# Use exponential decay adjacency matrix for eigendecomposition
c_adj_1 <- generate_exp_corr_matrix(n_row, n_col)

cat("Sumulating Data for Simulation 1 in Pixel Space ...\n")
tic()
sim1a_data <- run_pixel_to_freq_simulation(
  n_iter = n_iter,
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  effect_size = beta_effect,
  c_adj = c_adj_1,
  seed = seed
)
toc()

cat("Simulating Data for Simulation 1 in Frequency Space ...\n")
tic()
sim1b_data <- run_simulation_1b(
  n_iter = n_iter,
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  sparse_level = sparse_level,
  effect_size = b_effect,
  seed = seed
)
toc()


# ---------- Run Simulation 2 ----------
# Use 2-neighbor adjacency matrix for eigendecomposition
c_adj_2 <- generate_n_neighbor_matrix(n_row, n_col, d_max = 2)

cat("Simulating Data for Simulation 2 in Pixel Space ...\n")
tic()
sim2a_data <- run_pixel_to_freq_simulation(
  n_iter = n_iter,
  n_samples = n_samples,
  n_row = n_row,
  n_col = n_col,
  effect_size = beta_effect,
  c_adj = c_adj_2,
  seed = seed
)
toc()
