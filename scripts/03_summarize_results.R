proj_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper"

source(file.path(proj_dir, "R", "summarize_results.R"))


# Load simulation 1a fitted results
load(file = file.path(proj_dir, "results", "sim1a_lasso_results_241229.RData"))
load(file = file.path(proj_dir, "results", "sim1a_pixel_data_241229.RData"))

# Group mean difference in pixel space
with(
  sim1a_exp_pixel_data[[1]],
  vector_to_heatmap(colMeans(x[y == 1, ]) - colMeans(x[y == 0, ]))
)

# Order of eigenvalues
x_index <- with(
  sim1a_exp_pixel_data[[1]], seq_along(x) / length(x)
)

# Group mean difference in the freq space
with(
  sim1a_exp_pixel_data[[1]],
  plot_scatter(
    x = x_index,
    y = colMeans(x_freq[y == 1, ]) - colMeans(x_freq[y == 0, ]),
    xlab = "Order of eigenvalues (smallest to largest)"
  )
)

# True coefs in pixel space
vector_to_heatmap(sim1a_exp_pixel_data[[1]]$beta)

# True coefs in freq space
with(
  sim1a_exp_pixel_data[[1]],
  plot_scatter(
    x = x_index,
    y = b,
    xlab = "Order of eigenvalues (smallest to largest)"
  )
)

# Average coefs (pixel space, model fit in pixel space)
pixel_space_fitted_beta_coefs <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$pixel, function(x) x$coefs)
)
vector_to_heatmap(colMeans(pixel_space_fitted_beta_coefs))

# Average coefs (freq space, model fit in pixel space)
pixel_space_fitted_b_coefs <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$pixel, function(x) x$coefs_freq)
)
plot_scatter(
  x = x_index,
  y = colMeans(pixel_space_fitted_b_coefs),
  xlab = "Order of eigenvalues (smallest to largest)"
)

# Average coefs (pixel space, model fit in freq space)
freq_space_fitted_beta_coefs <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$freq, function(x) x$coefs_pixel)
)
vector_to_heatmap(colMeans(freq_space_fitted_beta_coefs))

# Average coefs (freq space, model fit in freq space)
freq_space_fitted_b_coefs <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$freq, function(x) x$coefs)
)
plot_scatter(
  x = x_index,
  y = colMeans(freq_space_fitted_b_coefs),
  xlab = "Order of eigenvalues (smallest to largest)"
)

# Percentage of p<0.05 for beta
beta_hdi_pvals <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$pixel, function(x) x$pvals)
)
vector_to_heatmap(colMeans(beta_hdi_pvals < 0.05))

# Percentage of p<0.05 for b
b_hdi_pvals <- do.call(
  "rbind",
  lapply(sim1a_lasso_results$freq, function(x) x$pvals)
)
plot_scatter(
  x = x_index,
  y = colMeans(b_hdi_pvals < 0.05),
  xlab = "Order of eigenvalues (smallest to largest)"
)
