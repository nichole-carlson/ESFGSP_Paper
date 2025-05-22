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

source(file.path(proj_dir, "R", "simulate_data.R"))
source(file.path(proj_dir, "R", "utils.R"))
source(file.path(proj_dir, "R", "fit_model.R"))


# ---------- Capture parameters ----------
option_list <- list(
  optparse::make_option(
    "--sim_id",
    type = "integer",
    help = "Simulation index"
  ),
  optparse::make_option(
    "--effect",
    type = "double",
    default = 0.1,
    help = "Effect size [default %default]"
  ),
  optparse::make_option(
    "--n_sample",
    type = "integer",
    default = 1000,
    help = "Number of samples per iteration [default %default]"
  ),
  optparse::make_option(
    "--out_dir",
    type = "character",
    help = "Directory to save the output file"
  ),
  optparse::make_option(
    "--seed",
    type = "character",
    default = NULL,
    help = "Random seed (optional)"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Convert seed if provided
seed <- if (!is.null(opt$seed)) as.numeric(opt$seed) else NULL

if (!is.null(seed)) set.seed(seed)


# ---------- Simulate Data ----------
# Exponential decay correlation matrix
cov_matrix <- gen_exp_corr_mat(n_row = 16, n_col = 16)

# Simulate X ~ MVN(0, cov_matrix)
x <- simulate_mvn_samples(n_sample = opt$n_sample, cov_matrix = cov_matrix)

# Generate 2D coeffcient matrix and flatten to 1D
beta_matrix <- matrix(0, nrow = 16, ncol = 16)
beta_matrix[5:12, 5:12] <- opt$effect
beta_vec <- as.vector(beta_matrix)

# Simulate y
y <- simulate_pixel_outcome(x = x, beta = beta_vec, from_pixel_space = TRUE)



# ---------- Calculate Eigenvector for Transformation ----------
# Derive adjacency matrix from correlation matrix
adj_matrix <- cov_matrix - diag(nrow(cov_matrix)) # adj matrix has diag = 0
# Eigendecomposition on MCM, C = adj matrix
eigen_mat_data <- eigen_decomp_mcm(adj_matrix)$vectors

sim_data <- list(x = x, y = y, e = eigen_mat_data)


# ---------- Fit LASSO Model in Pixel and Freq Spaces ----------
spaces <- c("pixel", "freq")
lambdas <- c("lambda.min", "lambda.1se")

lasso_results <- list()

for (space in spaces) {
  lasso_results[[space]] <- list()
  for (lam in lambdas) {
    if (space == "pixel") {
      # models fitted on pixel space
      lasso_results[[space]][[lam]] <- lasso_pixel_or_freq(
        x, y,
        lambda = lam
      )
    } else {
      # models fitted on freq space
      lasso_results[[space]][[lam]] <- lasso_pixel_or_freq(
        x, y,
        in_pixel_space = FALSE, e = eigen_mat_data, lambda = lam
      )
    }
  }
}


# ---------- Save as a rds file ----------
filename <- paste0("sim_", sprintf("%03d", opt$sim_id), ".rds")

saveRDS(
  list(data = sim_data, fit = lasso_results)
  file = file.path(opt$out_dir, filename)
)
