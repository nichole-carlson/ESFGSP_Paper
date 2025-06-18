# This script combines simulation results of iterations.
#
# results_001.rds, results_002.rds, ..., contains list of:
#   - pixel_min, pixel_1se, freq_min, freq_1se, each contains a list of:
#     - acc: classification accuracy
#     - auc: area under the ROC curve
#     - coefs: named vector of coefficients
#     - pval_pixel: vector of p-values for pixel-space coefs
#     - pval_freq: for freq-space coefs

packages <- c("optparse", "rprojroot")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

sapply(packages, require, character.only = TRUE)

proj_dir <- rprojroot::find_root(rprojroot::is_git_root)

# process_scalar_field, process_vector_field
source(file.path(proj_dir, "R", "summarize_results.R"))



# ---------- Arguments ----------
option_list <- list(
  optparse::make_option("--indir", type = "character"),
  optparse::make_option("--outdir", type = "character"),
  optparse::make_option("--sim_id", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# ---------- Get files ----------
file_paths <- list.files(
  opt$indir,
  pattern = "^results_\\d+\\.rds$", full.names = TRUE
)

if (length(file_paths) == 0) stop("No file found")

# Sort by simulation ID
get_sim_id <- function(filepath) {
  as.integer(gsub("\\D", "", basename(filepath)))
}
file_paths <- file_paths[order(get_sim_id(file_paths))]



# ---------- Process files ----------
auc_acc_rows <- list()
coef_array <- list()
pval_array <- list()

for (i in seq_along(file_paths)) {
  path <- file_paths[i]
  results <- readRDS(path)

  auc_acc_df <- results$auc_acc
  auc_acc_df$run <- i
  auc_acc_rows[[i]] <- auc_acc_df

  coef_array[[i]] <- results$coefs
  pval_array[[i]] <- results$pvals
}

auc_acc_df <- do.call(rbind, auc_acc_rows)
auc_acc_df <- auc_acc_df[, c("run", setdiff(names(auc_acc_df), "run"))]
coef_array <- simplify2array(coef_array)
pval_array <- simplify2array(pval_array)


# ---------- Save results ----------
combined_results <- list(
  auc_acc = auc_acc_df,
  coefs = coef_array,
  pvals = pval_array
)

out_path <- file.path(opt$outdir, paste0(opt$sim_id, "_combined_results.rds"))
saveRDS(combined_results, out_path)

cat("Saved", length(file_paths), "simulations to", out_path, "\n")
