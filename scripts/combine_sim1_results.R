# This script combines simulation results generated by `simulation_1.R`.
#
# Each simulation is saved as a .rds file containing:
#   - data: a list with simulated inputs
#       - x: n × p input matrix
#       - y: outcome vector
#       - e: p × p eigenbasis matrix
#       - beta
#       - argparse values (e.g., sim_id, effect size, seed)
#   - fit: a nested list of LASSO fits, structured as:
#       - space: "pixel" or "freq"
#         - lambda: "lambda.min" or "lambda.1se"
#           - acc: classification accuracy
#           - auc: area under the ROC curve
#           - coefs: named vector of coefficients
#           - pvals: p-values corresponding to coefficients
#
# This script extracts:
#   - x, y, beta, e, hparams, auc/acc, coefs, pvals.

packages <- c("optparse")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
sapply(packages, require, character.only = TRUE)


# ---------- Capture parameters ----------
option_list <- list(
  optparse::make_option(
    "--result_dir",
    type = "character",
    help = "Directory that saves .rds files for each iter"
  ),
  optparse::make_option(
    "--out_dir",
    type = "character",
    help = "Directory to save the combined .rds file"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# ---------- Extract Results ----------
# List all .rds files
files <- list.files(
  opt$result_dir,
  pattern = "^sim_\\d+\\.rds$",
  full.names = TRUE
)
sim_ids <- as.integer(gsub("\\D", "", basename(files))) # reorder
files <- files[order(sim_ids)]

# Get dimension info from the first file
n_iter <- length(files)
file_1 <- readRDS(files[1])
n <- nrow(file_1$data$x)
p <- ncol(file_1$data$x)
spaces <- names(file_1$fit)
lambdas <- names(file_1$fit[[1]])

# Initial empty vector / matrix / arrays
beta_vec <- NULL
hparams <- NULL
e <- NULL

x_arr <- array(NA_real_, dim = c(n, p, n_iter))
y_arr <- matrix(NA_real_, nrow = n_iter, ncol = n)

auc_acc_df <- data.frame(
  sim = numeric(),
  space = character(),
  lambda = character(),
  auc = numeric(),
  acc = numeric()
)

coef_arr <- array(
  NA_real_,
  dim = c(p, length(lambdas), length(spaces), n_iter),
  dimnames = list(
    coef = seq_len(p),
    lambda = lambdas,
    space = spaces,
    sim = seq_len(n_iter)
  )
)
p_arr <- coef_arr

for (i in seq_along(files)) {
  f <- files[i]
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res)) next

  if (i == 1) {
    beta_vec <- res$data$beta
    hparams <- res$data$hparams[c("effect", "n_sample")]
    e <- res$data$e
  }

  x_arr[, , i] <- res$data$x
  y_arr[i, ] <- res$data$y

  for (s in seq_along(spaces)) {
    space <- spaces[s]
    for (l in seq_along(lambdas)) {
      lam <- lambdas[l]
      res_sl <- res$fit[[space]][[lam]]
      auc_acc_df <- rbind(auc_acc_df, data.frame(
        sim = i,
        space = space,
        lambda = lam,
        auc = as.numeric(res_sl$auc),
        acc = res_sl$acc
      ))
      coef_arr[, l, s, i] <- res_sl$coefs
      p_arr[, l, s, i] <- res_sl$pvals
    }
  }
}



# Save results
saveRDS(
  list(
    x = x_arr,
    y = y_arr,
    beta = beta_vec,
    e = e,
    hparams = hparams,
    auc_acc = auc_acc_df,
    coefs = coef_arr,
    pvals = p_arr
  ),
  file = file.path(opt$out_dir, "sim1_combined_results.rds")
)
