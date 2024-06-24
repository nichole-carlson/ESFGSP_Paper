library(MASS)
library(glmnet)
library(pROC)

parent_dir <- "/Users/siyangren/Documents/ra-cida/ESFGSP_Paper/"
fig_dir <- file.path(parent_dir, "Figures")

# Define a base class for common functionality
sim_base <- setRefClass(
  "sim_base",
  fields = list(
    img_size = "numeric",
    num_samples = "numeric",
    x = "matrix",
    beta = "numeric",
    y = "numeric",
    w = "matrix",
    v = "matrix",
    seed = "numeric",
    coef_min = "matrix",
    coef_1se = "matrix"
  ),
  methods = list(
    initialize = function(img_size = 16, num_samples = 1000, seed = 42, beta_effect) {
      .self$img_size <- img_size
      .self$num_samples <- num_samples
      .self$seed <- seed
      set.seed(.self$seed)
      .self$gen_cov_matrix()
      .self$define_beta(beta_effect)
    },
    # Assume an exponential correlation structure
    gen_cov_matrix = function() {
      coords <- expand.grid(1:img_size, 1:img_size)
      dist_matrix <- as.matrix(dist(coords))
      .self$w <- exp(-dist_matrix)
    },
    gen_x = function() {
      stop("This method should be implemented by the subclass")
    },
    define_beta = function(beta_effect) {
      p <- .self$img_size^2
      .self$beta <- rep(0, p)
      center_idx <- as.vector(matrix(1:p, nrow = img_size)[5:12, 5:12])
      .self$beta[center_idx] <- beta_effect
    },
    check_p = function() {
      eta <- x %*% beta
      p <- 1 / (1 + exp(-eta))
      return(p)
    },
    gen_response = function() {
      .self$y <- rbinom(num_samples, 1, p)
    },
    fit_lasso = function() {
      train_idx <- sample(1:num_samples, size = floor(0.8 * num_samples))
      test_idx <- setdiff(1:num_samples, train_idx)

      x_train <- x[train_idx, , drop = FALSE]
      x_test <- x[-train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      y_test <- y[-train_idx]

      cv_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
      lambda_min <- cv_model$lambda.min
      lambda_1se <- cv_model$lambda.1se

      eval_perf <- function(l) {
        model <- glmnet(x_train, y_train, alpha = 1, lambda = l, family = "binomial")
        preds <- predict(model, newx = x_test, type = "response")[, 1]
        preds_bin <- ifelse(preds > 0.5, 1, 0)

        acc <- mean(preds_bin == y_test)
        auc <- auc(roc(y_test, preds))

        list(coefs = coef(model), acc = acc, auc = auc)
      }

      perf_min <- eval_perf(lambda_min)
      perf_1se <- eval_perf(lambda_1se)

      # Print performance metrics
      cat("Performance at lambda.min:\n")
      cat("Accuracy:", perf_min["acc"], "\n")
      cat("AUC:", perf_min["auc"], "\n")

      cat("Performance at lambda.1se:\n")
      cat("Accuracy:", perf_1se["acc"], "\n")
      cat("AUC:", perf_1se["auc"], "\n")

      .self$coef_min <- perf_min$coefs
      .self$coef_1se <- perf_1se$coefs
    },
    visualize_coef = function(lambda = "min") {
      coefs <- if (lambda == "min") coef_min else coef_1se
      beta_matrix <- matrix(beta, nrow = img_size)
      lasso_matrix <- matrix(coefs[-1], nrow = img_size)

      par(mfrow = c(1, 2))
      image(
        1:img_size, 1:img_size, beta_matrix,
        main = "True Coefficients",
        xlab = "", ylab = "",
        col = gray.colors(256, start = 0, end = 1, rev = TRUE)
      )
      image(
        1:img_size, 1:img_size, lasso_matrix,
        main = paste("LASSO Coefficients (", lambda, ")", sep = ""),
        xlab = "", ylab = "",
        col = gray.colors(256, start = 0, end = 1, rev = TRUE)
      )
    }
  )
)

# Define a class for Simulation 1
sim1 <- setRefClass(
  "sim1",
  contains = "sim_base",
  methods = list(
    initialize = function(img_size = 16, num_samples = 1000, seed = 42, beta_effect) {
      callSuper(img_size, num_samples, seed, beta_effect)
      .self$gen_x()
    },
    gen_x = function() {
      p <- img_size^2
      .self$x <- mvrnorm(num_samples, mu = rep(0, p), Sigma = w)
    }
  )
)


compare_beta_effects <- function(beta_effects = seq(0.1, 1.0, by = 0.1)) {
  par(mfrow = c(ceiling(length(beta_effects) / 2), 2))
  for (beta_effect in beta_effects) {
    sim1_obj <- sim1$new(beta_effect = beta_effect)
    p <- sim1_obj$check_p()
    title <- paste("beta_effect =", beta_effect)
    hist(p,
      breaks = 30, main = title, , xlim = c(0, 1),
      xlab = "Probability p", col = "lightblue", border = "black"
    )
  }
}

# Run the function to compare different beta effects
png(
  file = file.path(fig_dir, "sim1_p_dist.png"),
  width = 1600, height = 1200, res = 150
)
compare_beta_effects(beta_effects = c(1, 0.2, 0.1, 0.05, 0.01))
dev.off()
