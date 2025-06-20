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
source(file.path(proj_dir, "R", "simulate_data.R"))
# evaluate_lasso_performance, run_lasso_permutations, compute_permutation_pvals
# get_hdi_pvals
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

# Adjacency matrix
adj_matrix <- cov_matrix - diag(nrow(cov_matrix))

# Block coefficient vector (beta)
beta <- generate_block_coefs(c(16, 16), c(8, 8), opt$effect)

# Run simulation
data_sim <- simulate_data(
  opt$n_sample, cov_matrix, beta, adj_matrix,
  on_freq = FALSE
)

# Extract simulated data
y <- data_sim$outcome
b <- data_sim$coef_trans
transform_mat <- data_sim$transform_mat


# ---------- Fit LASSO Model in Pixel and Freq Spaces ----------
space_map <- c(pixel = "orig", freq = "trans")
lambda_choices <- c("lambda.min", "lambda.1se")
p <- ncol(cov_matrix)

# Placeholders
auc_acc_rows <- list()
coef_array <- array(
  NA_real_,
  dim = c(p, 2, length(space_map), length(lambda_choices)),
  dimnames = list(
    feature = seq_len(p),
    coef_type = c("orig", "trans"),
    space = names(space_map),
    lambda = lambda_choices,
  )
)
pval_array <- coef_array

for (space_name in names(space_map)) {
  on_freq <- (space_name == "freq")
  data_mat <- data_sim$x[, space_map[space_name]]

  for (lambda_choice in lambda_choices) {
    # Fit model on the space specified and lambda chosen
    # Return coef and pval on both spaces
    results <- fit_evaluate_lasso(
      data_mat, y, on_freq, transform_mat,
      lambda_choice = lambda_choice
    )

    # Store AUC and accuracy
    auc_acc_rows[[length(auc_acc_rows) + 1]] <- data.frame(
      space = space_name,
      lambda = lambda_choice,
      auc = results$auc,
      accuracy = results$accuracy
    )

    # Process coefs and pvals
    coef_array[, , space_name, lambda_choice] <- results$coef_arr
    pval_array[, , space_name, lambda_choice] <- results$pval_arr
  }
}

auc_acc_df <- do.call(rbind, auc_acc_rows)



# ---------- Save as a rds file ----------
data_path <- paste0("data_", sprintf("%03d", opt$sim_id), ".rds")
list(x = x, y = y, beta = beta, e = transform_mat, hparams = opt) |>
  saveRDS(file = file.path(opt$out_dir, data_path))

results_path <- paste0("results_", sprintf("%03d", opt$sim_id), ".rds")
list(auc_acc = auc_acc_df, coefs = coef_array, pvals = pval_array) |>
  saveRDS(file = file.path(opt$out_dir, results_path))
