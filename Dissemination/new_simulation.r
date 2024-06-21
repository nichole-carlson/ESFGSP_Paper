library(MASS)
library(glmnet)
library(pROC)

# Define a base class for common functionality
SimBase <- setRefClass(
  "SimBase",
  fields = list(
    img_size = "numeric",
    num_samples = "numeric",
    X = "matrix",
    beta = "numeric",
    y = "numeric",
    W = "matrix",
    V = "matrix",
    seed = "numeric",
    coef_min = "matrix",
    coef_1se = "matrix"
  ),
  methods = list(
    initialize = function(img_size = 16, num_samples = 1000, seed = 42) {
      .self$img_size <- img_size
      .self$num_samples <- num_samples
      .self$seed <- seed
      set.seed(.self$seed)
      .self$gen_cov_matrix()
    },
    gen_cov_matrix = function() {
      coords <- expand.grid(1:img_size, 1:img_size)
      dist_matrix <- as.matrix(dist(coords))
      .self$W <- exp(-dist_matrix)
    },
    gen_X = function() {
      stop("This method should be implemented by the subclass")
    },
    define_beta = function() {
      stop("This method should be implemented by the subclass")
    },
    gen_response = function() {
      eta <- X %*% beta
      p <- 1 / (1 + exp(-eta))
      .self$y <- rbinom(num_samples, 1, p)
    },
    fit_lasso = function() {
      train_idx <- sample(1:num_samples, size = floor(0.8 * num_samples))
      test_idx <- setdiff(1:num_samples, train_idx)

      x_train <- X[train_idx, , drop = FALSE]
      x_test <- X[-train_idx, , drop = FALSE]
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
    },
    run_simulation = function() {
      .self$gen_X()
      .self$define_beta()
      .self$gen_response()
      .self$fit_lasso()
    }
  )
)

# Define a class for Simulation 1
Sim1 <- setRefClass(
  "Sim1",
  contains = "SimBase",
  methods = list(
    gen_X = function() {
      p <- img_size^2
      .self$X <- mvrnorm(num_samples, mu = rep(0, p), Sigma = W)
    },
    define_beta = function() {
      p <- img_size^2
      .self$beta <- rep(0, p)
      center_idx <- as.vector(matrix(1:p, nrow = img_size)[5:12, 5:12])
      .self$beta[center_idx] <- 1
    }
  )
)

# Define a placeholder for Simulation 2
Sim2 <- setRefClass(
  "Sim2",
  contains = "SimBase",
  methods = list(
    gen_X = function() {
      p <- img_size^2
      eigen_decomp <- eigen(W)
      .self$V <- eigen_decomp$vectors
      X_temp <- matrix(rnorm(num_samples * p), nrow = num_samples, ncol = p)
      .self$X <- X_temp %*% V
    },
    define_beta = function() {
      p <- img_size^2
      b <- rep(0, p)
      b[sample(1:p, size = 10)] <- 1
      .self$beta <- t(V) %*% b
    }
  )
)

# Run Simulation 1
sim1 <- Sim1$new()
sim1$gen_response()
lasso_coef1 <- sim1$fit_lasso()
sim1$visualize_coef(lasso_coef1)

# Placeholder for running Simulation 2
# sim2 <- Sim2$new()
# sim2$gen_response()
# lasso_coef2 <- sim2$fit_lasso()
# sim2$visualize_coef(lasso_coef2)
