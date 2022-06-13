#' Add dummy predictors to the original predictor matrix
#'
#' Sample num_dummies dummy predictors from the univariate standard normal distribution and append them to the predictor matrix X.
#'
#' @param X Real valued predictor matrix.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#'
#' @return Enlarged predictor matrix.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 100
#' X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' add_dummies(X = X, num_dummies = p)
add_dummies <- function(X,
                        num_dummies) {
  # Error control
  if (!is.matrix(X)) {
    stop("'X' must be a matrix.")
  }

  if (!is.numeric(X)) {
    stop("'X' only allows numerical values.")
  }

  if (anyNA(X)) {
    stop("'X' contains NAs. Please remove or impute them before proceeding.")
  }

  if (length(num_dummies) != 1 ||
      num_dummies %% 1 != 0 ||
      num_dummies < 1) {
    stop("'num_dummies' must be an integer larger or equal to 1.")
  }

  # Number of rows of X
  n <- nrow(X)

  # Create matrix of dummy predictors
  dummies <- matrix(
    stats::rnorm(n * num_dummies),
    nrow = n,
    ncol = num_dummies,
    byrow = FALSE
  )

  # Append dummies to X
  X_Dummy <- cbind(X, dummies)

  return(X_Dummy)
}
