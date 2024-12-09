library(ggplot2)

proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"
figure_path <- file.path(proj_dir, "results", "figures")

# Source the functions for simulating data
source(file.path(proj_dir, "R", "simulate_data.R"))


# ---------- Choose the optimal effect size ----------
beta_effect_sizes <- c(1, 0.2, 0.1, 0.05, 0.01)
b_effect_sizes <- seq(0.1, 0.3, by = 0.05)

n_samples <- 1000
n_row <- 16
n_col <- 16
seed <- 42
c_adj <- generate_exp_corr_matrix(n_row, n_col)
sparse_level <- 0.1

# Choose the optimal effect size for beta (coef in pixel space)
beta_results <- data.frame()

for (effect in beta_effect_sizes) {
  one_iter <- run_pixel_to_freq_simulation(
    n_iter = 1,
    n_samples = n_samples,
    n_row = n_row,
    n_col = n_col,
    effect_size = effect,
    c_adj = c_adj,
    seed = seed
  )
  x <- one_iter$runs[[1]]$x
  beta <- one_iter$beta
  p <- 1 / (1 + exp(-(x %*% beta)))
  beta_results <- rbind(beta_results, data.frame(effect = effect, p = p))
}

beta_effect_histograms <- ggplot(beta_results, aes(x = p)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "lightblue",
    color = "black",
    boundary = 0,
    closed = "right"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(
    ~effect,
    ncol = 2,
    scales = "free_y",
    labeller = labeller(effect = function(x) paste("Effect =", x))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Probability", y = "Frequency")

# Choose optimal effect size for b (coef in frequency space)
b_results <- data.frame()

for (effect in b_effect_sizes) {
  one_iter <- run_simulation_1b(
    n_iter = 1,
    n_samples = n_samples,
    n_row = n_row,
    n_col = n_col,
    sparse_level = sparse_level,
    effect_size = effect,
    seed = seed
  )
  x_freq <- one_iter$runs[[1]]$x_freq
  b <- one_iter$b
  p <- 1 / (1 + exp(-(x_freq %*% b)))
  b_results <- rbind(b_results, data.frame(effect = effect, p = p))
}

b_effect_histograms <- ggplot(b_results, aes(x = p)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "lightblue",
    color = "black",
    boundary = 0,
    closed = "right"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(
    ~effect,
    ncol = 2,
    scales = "free_y",
    labeller = labeller(effect = function(x) paste("Effect =", x))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Probability", y = "Frequency")

ggsave(
  file.path(figure_path, "beta_effect_histograms.pdf"),
  beta_effect_histograms
)

ggsave(
  file.path(figure_path, "b_effect_histograms.pdf"),
  b_effect_histograms
)
