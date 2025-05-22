# This script combines results from simulation_1.R.
#
# Each iteration is saved as a separate .rds file, containing four fits:
#   - two in pixel space (lambda.min, lambda.1se)
#   - two in frequency space (lambda.min, lambda.1se)
#
# Each fit is a list with:
#   - acc: accuracy
#   - auc: area under the ROC curve
#   - coefs: estimated coefficients
#   - pvals: p-values corresponding to each coefficient
#
# This script extracts:
# (1) A summary data.frame with columns: sim_id, space, lambda, acc, auc
# (2) A long-format data.frame with columns: sim_id, space, lambda, coef, pval
#
# The combined results are saved as a .rds file as a list(auc_acc, coefs_pvals)

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

# Initialize lists to collect data
auc_acc_list <- list()
coefs_pvals_list <- list()

for (f in files) {
  sim_id <- as.integer(gsub("\\D", "", basename(f)))
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res)) next

  for (space in names(res)) {
    for (lam in names(res[[space]])) {
      res_sl <- res[[space]][[lam]]

      # Add to auc/acc list
      auc_acc_list[[length(auc_acc_list) + 1]] <- data.frame(
        sim_id = sim_id,
        space = space,
        lambda = lam,
        auc = as.numeric(res_sl$auc),
        acc = res_sl$acc
      )

      # Add to coefs/pvals list
      coefs_pvals_list[[length(coefs_pvals_list) + 1]] <- data.frame(
        sim_id = sim_id,
        space = space,
        lambda = lam,
        coef = res_sl$coefs,
        pval = res_sl$pvals
      )
    }
  }
}

# Combine into data.frames
auc_acc_df <- do.call(rbind, auc_acc_list)
coefs_pvals_df <- do.call(rbind, coefs_pvals_list)

# Save results
saveRDS(
  list(auc_acc = auc_acc_df, coefs_pvals = coefs_pvals_df),
  file = file.path(opt$out_dir, "sim1_combined_results.rds")
)
