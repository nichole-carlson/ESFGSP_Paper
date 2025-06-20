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
source(file.path(proj_dir, "R", "simulate_data.R"))
# fit_evaluate_lasso
source(file.path(proj_dir, "R", "fit_model.R"))


# ---------- Capture parameters ----------
option_list <- list(
  optparse::make_option("--sim_id", type = "integer"),
  optparse::make_option("--out_dir", type = "character"),
  optparse::make_option("--effect", type = "double"),
  optparse::make_option("--sparsity", type = "double"),
  optparse::make_option("--n_sample", type = "integer"),
  optparse::make_option("--seed", type = "integer")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

set.seed(opt$seed)


# ---------- Simulate Data ----------
# Diagonal covariance matrix (D)
cov_matrix <- diag(256 + 1 - seq_len(256))

# Adjacency matrix
exp_corr <- generate_exp_corr(n_row = 16, n_col = 16, rate = 1)
adj_mat <- exp_corr - diag(nrow(exp_corr))

# Sparse coefficient vector (b)
b <- generate_sparse_coefs(256, opt$sparsity, opt$effect)

data_sim <- simulate_data(opt$n_sample, cov_matrix, b, adj_mat, on_freq = TRUE)




# ---------- Fit LASSO Model in Pixel and Freq Spaces ----------
y <- data_sim$outcome
beta <- data_sim$coef_trans
transform_mat <- data_sim$transform_mat
p <- ncol(cov_matrix)

spaces <- list(pixel = data_sim$x_trans, freq = data_sim$x_orig)
lambda_choices <- c("lambda.min", "lambda.1se")

auc_acc_rows <- list()
coef_array <- array(
  NA_real_,
  dim = c(length(spaces), length(lambda_choices), 2, p),
  dimnames = list(
    space = names(spaces),
    lambda = lambda_choices,
    coef_type = c("orig", "trans"),
    feature = seq_len(p)
  )
)
pval_array <- array(
  NA_real_,
  dim = c(length(spaces), length(lambda_choices), 2, p),
  dimnames = list(
    space = names(spaces),
    lambda = lambda_choices,
    pval_type = c("orig", "trans"),
    feature = seq_len(p)
  )
)

# Process data
for (space_name in names(spaces)) {
  to_freq <- (space_name == "pixel")
  data_mat <- spaces[[space_name]]

  for (lambda_choice in lambda_choices) {
    results <- fit_evaluate_lasso(data_mat, y, lambda_choice = lambda_choice)

    # Store AUC and accuracy
    auc_acc_rows[[length(auc_acc_rows) + 1]] <- data.frame(
      space = space_name,
      lambda = lambda_choice,
      auc = results$auc,
      accuracy = results$accuracy
    )

    # Process coefs
    coefs <- results$coefs
    coef_trans <- transform_coef(coefs, transform_mat, to_freq)
    coef_array[space_name, lambda_choice, "orig", ] <- coefs
    coef_array[space_name, lambda_choice, "trans", ] <- coef_trans

    # Process permutation-based pvals
    perm_coefs <- results$perm_coefs
    perm_trans <- t(transform_coef(t(perm_coefs), transform_mat, to_freq))
    pvals_orig <- compute_permutation_pvals(perm_coefs, coefs)
    pvals_trans <- compute_permutation_pvals(perm_trans, coef_trans)
    pval_array[space_name, lambda_choice, "orig", ] <- pvals_orig
    pval_array[space_name, lambda_choice, "trans", ] <- pvals_trans
  }
}

auc_acc_df <- do.call(rbind, auc_acc_rows)


# ---------- Save as a rds file ----------
data_path <- paste0("data_", sprintf("%03d", opt$sim_id), ".rds")
list(x_freq = x_freq, y = y, b = b, e = transform_mat, hparams = opt) |>
  saveRDS(file = file.path(opt$out_dir, data_path))

results_file <- paste0("results_", sprintf("%03d", opt$sim_id), ".rds")
list(auc_acc = auc_acc_df, coefs = coef_array, pvals = pval_array) |>
  saveRDS(file = file.path(opt$out_dir, results_path))
