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
#
# The simulation is repeated 500 times.

proj_dir <- "/Users/siyangren/Documents/projects/ESFGSP_Paper"

source(file.path(proj_dir, "R", "simulate_data.R"))
source(file.path(proj_dir, "R", "utils.R"))
source(file.path(proj_dir, "R", "fit_model.R"))
