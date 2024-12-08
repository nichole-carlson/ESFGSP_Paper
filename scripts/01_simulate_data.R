proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"

# Source the functions for simulating data
source(file.path(proj_dir, "R", "simulate_data.R"))


# ---------- Choose the optimal effect size ----------
beta_effect_sizes <- c(1, 0.2, 0.1, 0.05, 0.01)
b_effect_sizes <- seq(0.1, 0.3, by = 0.05)




# Choose beta effect size for sim1
png(
  file = file.path(fig_dir, "sim1_p_dist.png"),
  width = 1200, height = 800, res = 150
)

beta_effects <- c(1, 0.2, 0.1, 0.05, 0.01)
par(mfrow = c(ceiling(length(beta_effects) / 2), 2))
for (effect in beta_effects) {
  sim1_1iter <- run_sim1(
    n_iter = 1, n_samples = 1000, img_size = 16,
    beta_effect_size = effect, seed = 42
  )
  x <- sim1_1iter$data$x
  beta <- sim1_1iter$meta_data$beta
  p <- 1 / (1 + exp(-(x %*% beta)))
  hist(
    p,
    breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
    xlab = "Probability p", col = "lightblue", border = "black"
  )
}

dev.off()

# Choose b effect for sim2
png(
  file = file.path(fig_dir, "sim2_p_dist.png"),
  width = 1200, height = 800, res = 150
)

b_effects <- seq(0.3, 0.1, by = -0.05)
par(mfrow = c(ceiling(length(b_effects) / 2), 2))
for (effect in b_effects) {
  sim2_1iter <- run_sim2(
    n_iter = 1, n_samples = 1000, img_size = 16,
    b_effect_size = effect, sparsity = 0.1, seed = 42
  )
  x_freq <- sim2_1iter$data$x_freq
  b <- sim2_1iter$meta_data$b
  p <- 1 / (1 + exp(-(x_freq %*% b)))
  hist(
    p,
    breaks = 30, main = paste("Effect =", effect), xlim = c(0, 1),
    xlab = "Probability p", col = "lightblue", border = "black"
  )
}

dev.off()


# ----- Run Simulations -----

n_sim <- 500
n_samples <- 1000
img_size <- 16
beta_effect <- 0.1
b_effect <- 0.2
b_sparsity <- 0.1
seed <- 42

cat("Simulating Data for", n_sim, "iterations of Simulation 1 ... \n")
tic()
sim1_data <- run_sim1(n_sim, n_samples, img_size, beta_effect, seed)
toc()

cat("Simulating Data for", n_sim, "iterations of Simulation 2 ... \n")
tic()
sim2_data <- run_sim2(n_sim, n_samples, img_size, b_effect, b_sparsity, seed)
toc()

sim1_1iter <- run_sim1(1, n_samples, img_size, beta_effect, seed)
sim2_1iter <- run_sim2(1, n_samples, img_size, b_effect, b_sparsity, seed)

filename <- paste0("simulated_data_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  n_sim, n_samples, img_size, beta_effect, b_effect, b_sparsity, seed,
  sim1_data, sim2_data, sim1_1iter, sim2_1iter,
  file = file.path(results_data_dir, filename)
)

filename2 <- paste0("simulated_data_1iter_", format(Sys.Date(), "%y%m%d"), ".RData")
save(
  n_sim, n_samples, img_size, beta_effect, b_effect, b_sparsity, seed,
  sim1_1iter, sim2_1iter,
  file = file.path(results_data_dir, filename2)
)
