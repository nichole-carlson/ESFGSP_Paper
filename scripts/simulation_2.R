# Simulation 2 settings
#
# X_freq and b are defined in the frequency space.
#
# X_freq is a 1000 x 256 matrix simulated from MVN(0, D), where D is a diagonal
# matrix. b has non-zero values in 10% entries.
#
# The transformation E is calculated from the exponential correlation matrix
# the same as Simulation 1.
#
# The binary outcome y is generated from X_freq and b using a logistic model.

packages <- c("optparse", "rprojroot")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

sapply(packages, require, character.only = TRUE)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)


# ---------- Import functions ----------
# generate_exp_corr, generate_image_data, generate_sparse_coefs,
# generate_outcomes
source(file.path(proj_dir, "R", "data_generation.R"))
# eigen_decomp_mcm, generate_e_matrix, transform_data
source(file.path(proj_dir, "R", "transformations.R"))
# fit_evaluate_lasso
source(file.path(proj_dir, "R", "fit_model.R"))


# ---------- Capture parameters ----------
option_list <- list(
  optparse::make_option("--sim_id", type = "integer"),
  optparse::make_option("--out_dir", type = "character"),
  optparse::make_option("--effect", type = "double"),
  optparse::make_option("--sparsity", type = "double"),
  optparse::make_option("--n_sample", type = "integer"),
  optparse::make_option("--seed", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

set.seed(opt$seed)


# ---------- Simulate Data ----------
# Diagonal covariance matrix (D)
cov_matrix <- diag(256 + 1 - seq_len(256))

# Transformation matrix (E)
exp_corr <- generate_exp_corr(n_row = 16, n_col = 16, rate = 1)
transform_mat <- generate_e_matrix(exp_corr)

# Simulate X_freq ~ MVN(0, D)
x_freq <- generate_image_data(opt$n_sample, cov_matrix)

# Sparse coefficient vector (b)
b <- generate_sparse_coefs(256, opt$sparsity, opt$effect)

# Simulate y
y <- generate_outcomes(x_freq, b)

# Calculate x (pixel)
x <- transform_data(x_freq, transform_mat, to_freq = FALSE)


# ---------- Fit LASSO Model in Pixel and Freq Spaces ----------
fit_results <- list(
  pixel_min = fit_evaluate_lasso(x, y, "lambda.min"),
  pixel_1se = fit_evaluate_lasso(x, y, "lambda.1se"),
  freq_min = fit_evaluate_lasso(x_freq, y, "lambda.min"),
  freq_1se = fit_evaluate_lasso(x_freq, y, "lambda.1se")
)


# ---------- Save as a rds file ----------
data_file <- paste0("data_", sprintf("%03d", opt$sim_id), ".rds")
saveRDS(
  list(x = x, y = y, beta = beta, e = transform_mat, hparams = opt),
  file.path(opt$out_dir, data_file)
)

results_file <- paste0("results_", sprintf("%03d", opt$sim_id), ".rds")
saveRDS(fit_results, file.path(opt$out_dir, results_file))
