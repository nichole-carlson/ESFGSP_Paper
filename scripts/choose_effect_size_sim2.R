# Choosing coefficient effect size for Simulation 2
#
# Simulate X_freq from a multivariate normal distribution, with covariance
# as a diagonal matrix. Then compute p = logistic (X_freq %*% b), where b
# is a vector with 10% non-zero values. The goal is to find a proper non-zero
# value such that the dist of p is approximately uniform between 0 and 1.

library(ggplot2)

proj_dir <- "/Users/siyangren/Documents/projects/ESFGSP_Paper"
fig_path <- file.path(proj_dir, "results", "figures")

source(file.path(proj_dir, "R", "simulate_data.R")) # simulate_mvn_samples


# Parameters
effect_sizes <- seq(0.01, 0.04, by = 0.005)
sparse_level <- 0.1
seed <- 42

# Covariance matrix D for X_freq, Dii = p + 1 - i
cov_matrix <- matrix(0, nrow = 256, ncol = 256)
diag(cov_matrix) <- 256 + 1 - seq_len(256)

# Simulate X ~ MVN(0, cov_matrix)
x <- simulate_mvn_samples(
  n_sample = 1000,
  cov_matrix = cov_matrix,
  seed = seed
)

# Generate coefficient vector with different effect sizes
b_matrix_list <- lapply(effect_sizes, function(effect) {
  set.seed(seed)
  vec <- rep(0, 256)
  nz_index <- sample(256, size = floor(256 * sparse_level))
  vec[nz_index] <- effect
  matrix(vec, ncol = 1)
})

# Calculate p = logistic(X %*% b)
probs <- lapply(b_matrix_list, function(mat) {
  1 / (1 + exp(-x %*% mat))
})

# Convert list into a data.frame for plotting
probs_df <- do.call(rbind, lapply(seq_along(probs), function(i) {
  data.frame(
    effect = effect_sizes[i],
    prob = probs[[i]]
  )
}))

# Plot histogram for p under each effect size
histograms <- ggplot2::ggplot(probs_df, aes(x = prob)) +
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
  file.path(fig_path, "sim2_p_histograms.pdf"),
  histograms
)
