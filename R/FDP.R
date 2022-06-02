#' False discovery proportion (FDP)
#'
#' Computes the FDP based on the estimated and the true regression coefficient vectors.
#'
#' @param beta_hat Estimated regression coefficient vector.
#' @param beta True regression coefficient vector.
#' @param eps Numerical zero.
#'
#' @return False discovery proportion (FDP).
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- c(Gauss_data$y)
#' beta <- Gauss_data$beta
#'
#' set.seed(1234)
#' res <- tknock(X, y)
#' beta_hat <- res$beta.selected
#'
#' FDP(beta_hat = beta_hat, beta = beta)
FDP <- function(beta_hat,
                beta,
                eps = .Machine$double.eps) {
  # Remove all dimension attributes of length one
  beta_hat <- drop(beta_hat)
  beta <- drop(beta)

  # Error control
  if (!is.vector(beta_hat)) {
    stop("'beta_hat' must be a vector.")
  }

  if (!is.numeric(beta_hat)) {
    stop("'beta_hat' only allows numerical values.")
  }

  if (anyNA(beta_hat)) {
    stop("'beta_hat' contains NAs. Please remove or impute them before proceeding.")
  }

  if (!is.vector(drop(beta))) {
    stop("'beta' must be a vector.")
  }

  if (!is.numeric(beta)) {
    stop("'beta' only allows numerical values.")
  }

  if (anyNA(beta)) {
    stop("'beta' contains NAs. Please remove or impute them before proceeding.")
  }

  if (length(beta_hat) != length(beta)) {
    stop("Length of beta_hat does not match length of beta.")
  }

  # Compute FDP
  num_selected_var <- sum(abs(beta_hat) > eps)
  num_false_positives <- sum(abs(beta) < eps & abs(beta_hat) > eps)

  if (num_selected_var == 0) {
    fdp <- 0
  } else {
    fdp <- num_false_positives / num_selected_var
  }

  return(fdp)
}
