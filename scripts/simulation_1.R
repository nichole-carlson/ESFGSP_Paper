# Simulation 1 settings
#
# X and beta are defined in the pixel space.
#
# X is a 1000 x 256 matrix simulated from a multivariate normal N(0, sigma).
# sigma_ij = exp(-d_ij), where d_ij is the Euclidean distance between pixel
# i and j on a 16 x 16 grid.
#
# The coefficent vector beta has non-zero values only within the central 8 x 8
# region.
#
# The binary outcome y is generated from X and beta using a logistic model.
# Two models are then fitted:
#   1. Fit y using X (in pixel space), and transform the estimated coefficients
#      to the frequency space.
#   2. Transform X into frequency space, fit y using X_freq, and convert the
#      estimated coefficients back to pixel space.
# Each model will be fitted on either lambda.min or lambda.1se

packages <- c("optparse", "rprojroot")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

sapply(packages, require, character.only = TRUE)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)


# ---------- Import functions ----------
# generate_exp_corr, generate_image_data, generate_block_coefs,
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
  optparse::make_option("--n_sample", type = "integer"),
  optparse::make_option("--seed", type = "integer")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

set.seed(opt$seed)


# ---------- Simulate Data ----------
# Exponential decay correlation matrix
cov_matrix <- generate_exp_corr(n_row = 16, n_col = 16, rate = 1)

# Transformation matrix (E)
transform_mat <- generate_e_matrix(cov_matrix)

# Simulate X ~ MVN(0, cov_matrix)
x <- generate_image_data(opt$n_sample, cov_matrix)

# Block coefficient vector (beta)
beta <- generate_block_coefs(c(16, 16), c(8, 8), opt$effect)

# Simulate y
y <- generate_outcomes(x, beta)

# Calculate x_freq
x_freq <- transform_data(x, transform_mat, to_freq = TRUE)


# ---------- Fit LASSO Model in Pixel and Freq Spaces ----------
tic()
fit_results <- list(
  pixel_min = fit_evaluate_lasso(x, y, family = "binomial", lambda = "lambda.min"),
  pixel_1se = fit_evaluate_lasso(x, y, family = "binomial", lambda = "lambda.1se"),
  freq_min = fit_evaluate_lasso(x_freq, y, family = "binomial", lambda = "lambda.min"),
  freq_1se = fit_evaluate_lasso(x_freq, y, family = "binomial", lambda = "lambda.1se")
)
toc()


# ---------- Save as a rds file ----------
data_file <- paste0("data_", sprintf("%03d", opt$sim_id), ".rds")
saveRDS(
  list(x = x, y = y, beta = beta, e = transform_mat, hparams = opt),
  file.path(opt$out_dir, data_file)
)

results_file <- paste0("results_", sprintf("%03d", opt$sim_id), ".rds")
saveRDS(fit_results, file.path(opt$out_dir, results_file))
