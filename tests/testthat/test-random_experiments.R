# 1
test_that("error control for inputs method and type works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(random_experiments(
    X = X,
    y = y,
    method = "test"
  ))

  expect_error(random_experiments(
    X = X,
    y = y,
    type = "test"
  ))
})

# 2
test_that("error control for input X works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  X_w_NA <- X
  X_w_NA[sample(prod(dim(X)), size = 100)] <- NA

  # Tests
  expect_error(random_experiments(
    X = drop(X[, 1]),
    y = y
  ),
  "'X' must be a matrix.",
  fixed = TRUE
  )

  expect_error(random_experiments(
    X = matrix(as.character(X), ncol = ncol(X)),
    y = y
  ),
  "'X' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(random_experiments(
    X = matrix(as.factor(X), ncol = ncol(X)),
    y = y
  ),
  "'X' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X_w_NA,
      y = y
    ),
    "'X' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )
})

# 3
test_that("error control for input y works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  y_w_NA <- y
  y_w_NA[sample(length(y), size = 10)] <- NA

  # Tests
  expect_error(random_experiments(
    X = X,
    y = cbind(y, y)
  ),
  "'y' must be a vector.",
  fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = as.character(y)
    ),
    "'y' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(random_experiments(
    X = X,
    y = matrix(as.factor(y), ncol = 1)
  ),
  "'y' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y_w_NA
    ),
    "'y' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = rbind(X, X),
      y = y
    ),
    "Number of rows in X does not match length of y.",
    fixed = TRUE
  )
})

# 4
test_that("error control for input K works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      K = 1
    ),
    "The number of random experiments 'K' must be an integer larger or equal to 2.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      K = 20.3
    ),
    "The number of random experiments 'K' must be an integer larger or equal to 2.",
    fixed = TRUE
  )
})

# 5
test_that("error control for input num_dummies works when method = 'trex'", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  num_dummies <- 2 * p

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      num_dummies = num_dummies + 1e-4,
      method = "trex"
    ),
    "'num_dummies' must be an integer larger or equal to 1.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      num_dummies = 0,
      method = "trex"
    ),
    "'num_dummies' must be an integer larger or equal to 1.",
    fixed = TRUE
  )
})

# 6
test_that("error control for input num_dummies works when method = 'trex+GVS'", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  num_dummies <- 2 * p

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      num_dummies = num_dummies + 1e-4,
      method = "trex+GVS"
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      num_dummies = 1.5 * p,
      method = "trex+GVS"
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      num_dummies = 0,
      method = "trex+GVS"
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )
})

# 7
test_that("the input value of 'T_stop' is valid", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  num_dummies <- 2 * p

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      T_stop = 0,
      num_dummies = num_dummies
    ),
    paste0(
      "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
      num_dummies,
      "."
    ),
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      T_stop = num_dummies + 1,
      num_dummies = num_dummies
    ),
    paste0(
      "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
      num_dummies,
      "."
    ),
    fixed = TRUE
  )
})

# 8
test_that("error control for input corr_max works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      method = "trex+GVS",
      corr_max = -0.1
    ),
    "'corr_max' must have a value between zero and one.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      method = "trex+GVS",
      corr_max = 1.1
    ),
    "'corr_max' must have a value between zero and one.",
    fixed = TRUE
  )
})

# 9
test_that("error control for input lambda_2_lars works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      method = "trex+GVS",
      lambda_2_lars = c(1, 5, 100)
    ),
    "'lambda_2_lars' must be a number larger than zero.",
    fixed = TRUE
  )

  expect_error(
    random_experiments(
      X = X,
      y = y,
      method = "trex+GVS",
      lambda_2_lars = -3
    ),
    "'lambda_2_lars' must be a number larger than zero.",
    fixed = TRUE
  )
})

# 10
test_that("error control for input parallel_max_cores works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      parallel_process = TRUE,
      parallel_max_cores = 1
    ),
    "For parallel processing at least two workers have to be registered:
         'parallel_max_cores' must be an integer larger or equal to 2.",
    fixed = TRUE
  )
})

# 11
test_that("reasonable number of workers is registered for parallel processing", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  K <- 20
  parallel_max_cores <- 21

  # Tests
  expect_message(
    random_experiments(
      X = X,
      y = y,
      K = K,
      parallel_process = TRUE,
      parallel_max_cores = parallel_max_cores
    ),
    paste0(
      "For computing ",
      K,
      " random experiments, it is not useful to register ",
      parallel_max_cores,
      " workers. Setting parallel_max_cores = ",
      min(K, max(
        1, parallel::detectCores(logical = FALSE)
      )),
      " (# physical cores) ...\n"
    ),
    fixed = TRUE
  )
})

# 12
test_that("error control for input parallel_process works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_message(
    random_experiments(
      X = X,
      y = y,
      parallel_process = TRUE
    ),
    "Computing random experiments in parallel...
            Note that this is only advantageous if you have at least a few thousand predictors and/or data points in 'X'.
            Otherwise, the overhead will slow down the computations in parallel. Thus, for small data sizes it is better
            to set parallel_process = FALSE.

      Be careful!",
    fixed = TRUE
  )
})

# 13
test_that("error control for input lars_state_list works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Run K = 21 random experiments with T_stop = 1
  rand_ex_T_stop_1 <- random_experiments(
    X = X,
    y = y,
    K = 21,
    T_stop = 1
  )

  # Extract lars_state_list
  lars_state_list <- rand_ex_T_stop_1$lars_state_list

  # Tests
  expect_error(
    random_experiments(
      X = X,
      y = y,
      K = 20,
      T_stop = 2,
      lars_state_list = lars_state_list
    ),
    "Length of 'lars_state_list' must be equal to number of random experiments 'K'.",
    fixed = TRUE
  )
})
