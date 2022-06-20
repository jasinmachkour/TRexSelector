# 1
test_that("error control for input X works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  p <- ncol(X)
  num_dummies <- 2 * p
  X_w_NA <- X
  X_w_NA[sample(prod(dim(X)), size = 100)] <- NA

  # Tests
  expect_error(add_dummies_GVS(
    X = drop(X[, 1]),
    num_dummies = num_dummies
  ),
  "'X' must be a matrix.",
  fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = matrix(as.character(X), ncol = ncol(X)),
      num_dummies = num_dummies
    ),
    "'X' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = matrix(as.factor(X), ncol = ncol(X)),
      num_dummies = num_dummies
    ),
    "'X' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = X_w_NA,
      num_dummies = num_dummies
    ),
    "'X' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )
})

# 2
test_that("error control for input num_dummies works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  p <- ncol(X)
  num_dummies <- 2 * p

  # Tests
  expect_error(
    add_dummies_GVS(
      X = X,
      num_dummies = num_dummies + 1e-4
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = X,
      num_dummies = 1.5 * p
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = X,
      num_dummies = 0
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )
})

# 3
test_that("error control for input corr_max works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  p <- ncol(X)
  num_dummies <- 2 * p

  # Tests
  expect_error(
    add_dummies_GVS(
      X = X,
      num_dummies = num_dummies,
      corr_max = -0.1
    ),
    "'corr_max' must have a value between zero and one.",
    fixed = TRUE
  )

  expect_error(
    add_dummies_GVS(
      X = X,
      num_dummies = num_dummies,
      corr_max = 1.1
    ),
    "'corr_max' must have a value between zero and one.",
    fixed = TRUE
  )
})

# 4
test_that("dimension of output (i.e., predictor matrix containing dummies) is p + num_dummies", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  p <- ncol(X)
  num_dummies <- 2 * p

  # Create dummies for the T-Rex+GVS selector
  X_Dummy <- add_dummies_GVS(
    X = X,
    num_dummies = num_dummies
  )

  # Tests
  expect_true(ncol(X_Dummy) == p + num_dummies)
})
