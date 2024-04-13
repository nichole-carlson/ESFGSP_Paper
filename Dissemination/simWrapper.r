# List of required packages
packages <- c("foreach", "doParallel", "parallel")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load packages
sapply(packages, require, character.only = TRUE)


simWrapper <- function(n_sim, f_sim, TF_parallel = FALSE,
                       n_cores = 2, list_export = NULL,
                       list_package = c(), f_clusterCall = function() {}) {
  # Generate random seeds for each simulation to ensure reproducibility
  v_seed <- sample(2^31 - 1, n_sim)

  if (TF_parallel) {
    # Start a cluster with the specified number of cores
    cl <- makeCluster(n_cores)
    # Register the parallel backend to be used with the cluster
    registerDoParallel(cl)

    try({
      # Export specified variables to all workers in the cluster
      clusterExport(cl, list_export)
      # Load the required packages on each worker node
      clusterCall(cl, function() sapply(list_package, library, character.only = TRUE))
      # Run any user-defined initializations on each worker
      clusterCall(cl, f_clusterCall)

      # Perform the simulations in parallel
      output_full <- foreach(i = 1:n_sim, .combine = rbind) %dopar% {
        set.seed(v_seed[i]) # Set the seed for reproducibility
        return(f_sim(i)) # Run the simulation function
      }
    })
    # Stop the cluster once done to free up system resources
    stopCluster(cl)
  } else {
    # Perform the simulations serially
    output_full <- foreach(i = 1:n_sim, .combine = rbind) %do% {
      set.seed(v_seed[i]) # Set the seed for reproducibility
      return(f_sim(i)) # Run the simulation function
    }
  }

  # Return the combined results of all simulations
  return(output_full)
}
