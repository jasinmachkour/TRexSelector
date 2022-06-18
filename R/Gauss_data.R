#' Toy data generated from a Gaussian linear model
#'
#' A data set containing a predictor matrix X with n = 50 observations
#' and p = 100 variables (predictors), and a sparse parameter vector beta
#' with associated support vector.
#'
#' @format A list containing a matrix X and vectors y, beta, and support:
#' \describe{
#'   \item{X}{Predictor matrix, n = 50, p = 100.}
#'   \item{y}{Response vector.}
#'   \item{beta}{Parameter vector.}
#'   \item{support}{Support vector.}
#' }
#'
#' @importFrom stats rnorm
#'
#' @examples
#' # Generated as follows:
#' set.seed(789)
#' n <- 50
#' p <- 100
#' X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' beta <- c(rep(5, times = 3), rep(0, times = 97))
#' support <- beta > 0
#' y <- X %*% beta + stats::rnorm(n)
#' Gauss_data <- list(
#'   X = X,
#'   y = y,
#'   beta = beta,
#'   support = support
#' )
"Gauss_data"
