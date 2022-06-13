# 1
test_that("error control for input beta_hat works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- c(Gauss_data$y)
  beta <- Gauss_data$beta

  set.seed(1234)
  res <- tknock(X, y)
  beta_hat <- res$selected_var

  beta_hat_w_NA <- beta_hat
  beta_hat_w_NA[sample(length(beta_hat_w_NA), size = 10)] <- NA

  # Tests
  expect_error(FDP(
    beta_hat = cbind(beta_hat, beta_hat),
    beta = beta
  ),
  "'beta_hat' must be a vector.",
  fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = as.character(beta_hat),
      beta = beta
    ),
    "'beta_hat' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = matrix(as.factor(beta_hat), ncol = 1),
      beta = beta
    ),
    "'beta_hat' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = beta_hat_w_NA,
      beta = beta
    ),
    "'beta_hat' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )
})

# 2
test_that("error control for input beta works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- c(Gauss_data$y)
  beta <- Gauss_data$beta

  set.seed(1234)
  res <- tknock(X, y)
  beta_hat <- res$selected_var

  beta_w_NA <- beta
  beta_w_NA[sample(length(beta_w_NA), size = 10)] <- NA

  # Tests
  expect_error(FDP(
    beta_hat = beta_hat,
    beta = cbind(beta, beta)
  ),
  "'beta' must be a vector.",
  fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = beta_hat,
      beta = as.character(beta)
    ),
    "'beta' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = beta_hat,
      beta = matrix(as.factor(beta), ncol = 1)
    ),
    "'beta' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = beta_hat,
      beta = beta_w_NA
    ),
    "'beta' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )

  expect_error(
    FDP(
      beta_hat = c(beta_hat, beta_hat),
      beta = beta
    ),
    "Length of beta_hat does not match length of beta.",
    fixed = TRUE
  )
})

# 3
test_that("the value of FDP is an element of the interval [0, 1]", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- c(Gauss_data$y)
  beta <- Gauss_data$beta

  set.seed(1234)
  res <- tknock(X, y)
  beta_hat <- res$selected_var

  fdp <- FDP(
    beta_hat = beta_hat,
    beta = beta
  )

  # Tests
  expect_true(fdp >= 0 && fdp <= 1)
})
