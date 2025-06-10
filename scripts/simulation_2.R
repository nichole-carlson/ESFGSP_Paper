# Simulation 2 settings
#
# X_freq and b are defined in the frequency space.
#
# X_freq is a 1000 x 256 matrix simulated from MVN(0, D), where D is a diagonal
# matrix. b has non-zero values in 10% entries.
#
# The binary outcome y is generated from X_freq and b using a logistic model.

packages <- c("optparse", "rprojroot")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
sapply(packages, require, character.only = TRUE)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)

# simulate_mvn_samples, simulate_pixel_outcome
source(file.path(proj_dir, "R", "simulate_data.R"))
# gen_exp_corr_mat, eigen_decomp_cmc
source(file.path(proj_dir, "R", "utils.R"))


# ---------- Capture parameters ----------
option_list <- list(
  optparse::make_option("--sim_id", type = "integer"),
  optparse::make_option("--out_dir", type = "character"),
  optparse::make_option("--effect", type = "double", default = 0.03),
  optparse::make_option("--sparsity", type = "double", default = 0.1),
  optparse::make_option(
    "--n_sample",
    type = "integer",
    default = 1000,
    help = "Samples per iteration"
  ),
  optparse::make_option("--seed", type = "character", default = NULL)
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Convert seed if provided
seed <- if (!is.null(opt$seed)) as.numeric(opt$seed) else NULL

if (!is.null(seed)) set.seed(seed)


# ---------- Simulate Data ----------
# Diagonal covariance matrix (D)
cov_matrix <- matrix(0, nrow = 256, ncol = 256)
diag(cov_matrix) <- 256 + 1 - seq_len(256)

# Simulate X_freq ~ MVN(0, D)
x <- simulate_mvn_samples(n_sample = opt$n_sample, cov_matrix = cov_matrix)

# Coef matrix
b <- {
  vec <- rep(0, 256)
  nz_index <- sample(256, size = floor(256 * opt$sparsity))
  vec[nz_index] <- opt$effect
  matrix(vec, ncol = 1)
}

# covariance matrix in pixel space (Sigma)
cov_sigma <- gen_exp_corr_mat(n_row = 16, n_col = 16)
adj_matrix <- cov_sigma - diag(nrow(cov_sigma))
eigen_mat_model <- eigen_decomp_mcm(adj_matrix)$vectors

# Simulate y
y <- simulate_pixel_outcome(x, b, from_pixel_space = TRUE, e = eigen_mat_model)
