#' True positive proportion (TPP)
#'
#' @param beta_hat Estimated support vector.
#' @param beta True support vector.
#' @param eps Numerical zero.
#'
#' @return TPP.
#'
#' @export
TPP = function(beta_hat, beta, eps = .Machine$double.eps) {
  if (sum(is.na(beta_hat)) == length(beta_hat)) {
    tpp = NA
  } else if (sum(abs(beta) > eps) == 0) {
    tpp = 0
  } else{
    tpp = 100 * (sum(abs(beta) > eps &
                       abs(beta_hat) > eps) / sum(abs(beta) > eps))
  }
  return(tpp)
}
