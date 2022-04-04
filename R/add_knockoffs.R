#' Add knockoffs to the original predictor matrix
#'
#' Sample L_val knockoff vectors from the univariate standard normal distribution and append them to the predictor matrix X.
#'
#' @param X Real valued Predictor matrix.
#' @param L_val Number of knockoffs.
#' @param cor.structure TRUE/FALSE.
#' @param empirical TRUE/FALSE.
#' @param eps Numerical zero.
#'
#' @return Enlarged predictor matrix, i.e., original predictor matrix and knockoffs.
#'
#' @import stats
#' @importFrom  Matrix qr
#' @importFrom  Matrix nearPD
#' @import MASS
#' @import mvnfast
#'
#' @export
#'
#' @examples
#' n = 50
#' p = 100
#' add_knockoffs(X = matrix(rnorm(n * p), nrow = n, ncol = p), L_val = p)
add_knockoffs = function(X,
                         L_val,
                         cor.structure = FALSE,
                         empirical = FALSE,
                         eps = .Machine$double.eps) {
  g = 1
  n = nrow(X)
  p = ncol(X)
  mu = rep(0, times = p)
  for (i in seq(g)) {
    if (L_val != p) {
      X_surrogate = matrix(
        stats::rnorm(n * L_val),
        nrow = n,
        ncol = L_val,
        byrow = FALSE
      )
    } else{
      if (cor.structure) {
        Sig = stats::cov(X)
        if (Matrix::qr(Sig)$rank != p) {
          Sig = as.matrix(Matrix::nearPD(Sig)$mat)
        }
        if (empirical) {
          X_surrogate = MASS::mvrnorm(n,
                                      mu = mu,
                                      Sigma = Sig,
                                      empirical = empirical)
        } else{
          X_surrogate = mvnfast::rmvn(
            n,
            mu = mu,
            sigma = Sig,
            ncores = 1,
            isChol = FALSE,
            A = NULL
          )
        }
      } else{
        X_surrogate = matrix(rnorm(n * p),
                             nrow = n,
                             ncol = p,
                             byrow = FALSE)
      }
    }
    X_Knock = cbind(X, X_surrogate)
  }
  return(X_Knock)
}
