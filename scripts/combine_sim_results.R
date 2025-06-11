# This script combines simulation results of iterations.
#
# results_001.rds, results_002.rds, ..., contains list of:
#   - pixel_min, pixel_1se, freq_min, freq_1se, each contains a list of:
#     - acc: classification accuracy
#     - auc: area under the ROC curve
#     - coefs: named vector of coefficients
#     - pvals: p-values corresponding to coefficients

library(optparse)


# ---------- Arguments ----------
option_list <- list(
  optparse::make_option("--indir", type = "character"),
  optparse::make_option("--outdir", type = "character"),
  optparse::make_option("--sim_id", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# ---------- Get files ----------
results_files <- list.files(
  opt$indir,
  pattern = "^results_\\d+\\.rds$", full.names = TRUE
)

# Sort by simulation ID
get_sim_id <- function(files) {
  as.integer(gsub("\\D", "", basename(files)))
}
results_files <- results_files[order(get_sim_id(results_files))]

# Get number of iterations
n_iter <- length(results_files)
if (n_iter == 0) stop("No results files found")


# ---------- Simulation 1 hparams ----------
p <- 256
spaces <- c("pixel", "freq")
lambdas <- c("lambda.min", "lambda.1se")



# ---------- Process files ----------
auc_acc_rows <- list()
coef_arrays <- list()
pval_arrays <- list()

# Create mapping from flat names to array indices
fit_mapping <- expand.grid(
  lambda_idx = seq_along(lambdas), space_idx = seq_along(spaces)
)
fit_mapping$fit_name <- paste0(
  spaces[fit_mapping$space_idx], "_",
  sub("lambda\\.", "", lambdas[fit_mapping$lambda_idx])
)

for (i in seq_len(n_iter)) {
  results <- readRDS(results_files[i])

  # auc acc
  metrics_list <- lapply(fit_mapping$fit_name, function(name) {
    parts <- strsplit(name, "_")[[1]]
    data.frame(
      sim = i,
      space = parts[1],
      lambda = paste0("lambda.", parts[2]),
      auc = as.numeric(results[[name]]$auc),
      acc = results[[name]]$acc
    )
  })
  auc_acc_rows <- c(auc_acc_rows, metrics_list)

  # coefs and pvalues
  coef_matrix <- array(NA_real_, dim = c(p, length(lambdas), length(spaces)))
  pval_matrix <- array(NA_real_, dim = c(p, length(lambdas), length(spaces)))

  for (j in seq_len(nrow(fit_mapping))) {
    fit_name <- fit_mapping$fit_name[j]
    lambda_idx <- fit_mapping$lambda_idx[j]
    space_idx <- fit_mapping$space_idx[j]

    coef_matrix[, lambda_idx, space_idx] <- results[[fit_name]]$coefs
    pval_matrix[, lambda_idx, space_idx] <- results[[fit_name]]$pvals
  }

  coef_arrays[[i]] <- coef_matrix
  pval_arrays[[i]] <- pval_matrix
}


# ---------- Combine results ----------
auc_acc_df <- do.call(rbind, auc_acc_rows)

# Convert to 4D arrays
coef_arr <- array(
  unlist(coef_arrays),
  dim = c(p, length(lambdas), length(spaces), n_iter),
  dimnames = list(
    seq_len(p),
    lambdas,
    spaces,
    seq_len(n_iter)
  )
)


# ---------- Save results ----------
combined_results <- list(
  auc_acc = auc_acc_df,
  coefs = coef_arr,
  pvals = pval_arr
)
output_file <- file.path(
  opt$outdir, paste0(opt$sim_id, "_combined_results.rds")
)
saveRDS(combined_results, output_file)

cat("Saved", n_iter, "simulations to", output_file, "\n")
