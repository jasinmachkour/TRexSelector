# 1
test_that("error control for inputs method and type works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_error(lm_dummy(
    X = X,
    y = y,
    method = "test"
  ))

  expect_error(lm_dummy(
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
  expect_error(lm_dummy(
    X = drop(X[, 1]),
    y = y
  ),
  "'X' must be a matrix.",
  fixed = TRUE
  )

  expect_error(lm_dummy(
    X = matrix(as.character(X), ncol = ncol(X)),
    y = y
  ),
  "'X' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(lm_dummy(
    X = matrix(as.factor(X), ncol = ncol(X)),
    y = y
  ),
  "'X' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(
    lm_dummy(
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
  expect_error(lm_dummy(
    X = X,
    y = cbind(y, y)
  ),
  "'y' must be a vector.",
  fixed = TRUE
  )

  expect_error(lm_dummy(
    X = X,
    y = as.character(y)
  ),
  "'y' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(lm_dummy(
    X = X,
    y = matrix(as.factor(y), ncol = 1)
  ),
  "'y' only allows numerical values.",
  fixed = TRUE
  )

  expect_error(
    lm_dummy(
      X = X,
      y = y_w_NA
    ),
    "'y' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )

  expect_error(lm_dummy(
    X = rbind(X, X),
    y = y
  ),
  "Number of rows in X does not match length of y.",
  fixed = TRUE
  )
})

# 4
test_that(
  "T-LARS model is an object of class tlars_cpp and stays an object of the same class after T-LARS steps",
  {
    # Setup and data generation
    data("Gauss_data")
    X <- Gauss_data$X
    y <- drop(Gauss_data$y)

    # Run one random experiment with T_stop = 1
    mod_tlars <- lm_dummy(
      X = X,
      y = y,
      T_stop = 1
    )

    # Tests
    expect_true(methods::is(object = mod_tlars, class2 = tlars::tlars_cpp))

    # Continue random experiment at the previous T-LARS state (warm start) until T_stop = 3 dummies have been selected
    mod_tlars <- lm_dummy(
      X = X,
      y = y,
      model_tlars = mod_tlars,
      T_stop = 3
    )

    # Tests
    expect_true(methods::is(object = mod_tlars, class2 = tlars::tlars_cpp))
  }
)

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
    lm_dummy(
      X = X,
      y = y,
      num_dummies = num_dummies + 1e-4,
      method = "trex"
    ),
    "'num_dummies' must be an integer larger or equal to 1.",
    fixed = TRUE
  )

  expect_error(
    lm_dummy(
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
    lm_dummy(
      X = X,
      y = y,
      num_dummies = num_dummies + 1e-4,
      method = "trex+GVS"
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    lm_dummy(
      X = X,
      y = y,
      num_dummies = 1.5 * p,
      method = "trex+GVS"
    ),
    "'num_dummies' must be a positive integer multiple of the total number of original predictors in X.",
    fixed = TRUE
  )

  expect_error(
    lm_dummy(
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
    lm_dummy(
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
    lm_dummy(
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
    lm_dummy(
      X = X,
      y = y,
      method = "trex+GVS",
      corr_max = -0.1
    ),
    "'corr_max' must have a value between zero and one.",
    fixed = TRUE
  )

  expect_error(
    lm_dummy(
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
    lm_dummy(
      X = X,
      y = y,
      method = "trex+GVS",
      lambda_2_lars = c(1, 5, 100)
    ),
    "'lambda_2_lars' must be a number larger than zero.",
    fixed = TRUE
  )

  expect_error(
    lm_dummy(
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
test_that("the user is informed that the entire solution path is computed if early_stop = FALSE", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)

  # Tests
  expect_message(
    lm_dummy(
      X = X,
      y = y,
      method = "trex",
      early_stop = FALSE
    ),
    "'T_stop' is ignored. Computing the entire solution path...",
    fixed = TRUE
  )

  expect_message(
    lm_dummy(
      X = X,
      y = y,
      method = "trex+GVS",
      early_stop = FALSE
    ),
    "'T_stop' is ignored. Computing the entire solution path...",
    fixed = TRUE
  )
})
