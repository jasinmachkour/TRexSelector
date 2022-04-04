#' False discovery proportion (FDP)
#'
#' @param beta_hat Estimated support vector.
#' @param beta True support vector.
#' @param eps Numerical zero.
#'
#' @return FDP
#'
#' @export
FDP = function(beta_hat, beta, eps = .Machine$double.eps) {
  if (sum(is.na(beta_hat)) == length(beta_hat)) {
    fdp = NA
  } else if (sum(abs(beta_hat) > eps) == 0) {
    fdp = 0
  } else{
    fdp = 100 * (sum(abs(beta) < eps &
                       abs(beta_hat) > eps) / sum(abs(beta_hat) > eps))
  }
  return(fdp)
}
