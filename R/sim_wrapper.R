# List of required packages
packages <- c("foreach", "doParallel", "parallel")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load packages
sapply(packages, require, character.only = TRUE)


parallel_wrapper <- function(task_fn, args_list, cores = NULL, pkgs = NULL) {
  # Set the number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  # Set up parallel backend
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # Run tasks in paralle using foreach
  results <- foreach::foreach(
    args = iterators::iter(args_list),
    .packages = pkgs
  ) %dopar% {
    do.call(task_fn, args)
  }

  # Stop the cluster
  parallel::stopCluster(cl)

  return(results)
}
