# Choosing coefficient effect size for Simulation 1
#
# Simulate X once from a multivariate normal distribution.
# Then compute p = logistic(X %*% beta) under different effect sizes.
# The goal is to find an effect size such that the resulting distribution of p
# is approximately uniform between 0 and 1.

library(ggplot2)

proj_dir <- "/Users/siyangren/Documents/projects/ESFGSP_Paper"
fig_path <- file.path(proj_dir, "results", "figures")

source(file.path(proj_dir, "R", "simulate_data.R"))
source(file.path(proj_dir, "R", "utils.R"))


# Parameters
effect_sizes <- c(1, 0.2, 0.1, 0.05, 0.01)
seed <- 42

set.seed(seed)
# Covariance matrix for X
cov_matrix <- gen_exp_corr_mat(n_row = 16, n_col = 16)

# Simulate X ~ MVN(0, cov_matrix)
x <- simulate_mvn_samples(
  n_sample = 1000,
  cov_matrix = cov_matrix,
  seed = seed
)

# Generate 2D coefficient matrices with central regions different
beta_matrix_list <- lapply(effect_sies, function(effect) {
  mat <- matrix(0, nrow = 16, ncol = 16)
  mat[5:12, 5:12] <- effect
  mat
})

# Calculate p = logistic(X %*% beta)
probs <- lapply(beta_matrix_list, function(mat) {
  vec <- as.vector(mat)
  1 / (1 + exp(-x %*% vec))
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
  file.path(fig_path, "sim1_p_histograms.pdf"),
  histograms
)
