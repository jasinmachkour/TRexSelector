#' Computes the Deflated Relative Occurrences
#'
#' Computes the matrix of deflated relative occurrences for all variables (i.e., j = 1,..., p) and for T = 1, ..., T_stop.
#'
#' @param p Number of candidate variables.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param num_dummies Number of dummies
#' @param phi_T_mat Matrix of relative occurrences for all variables (i.e., j = 1,..., p) and for T = 1, ..., T_stop.
#' @param Phi Vector of relative occurrences for all variables (i.e., j = 1,..., p) at T = T_stop.
#' @param eps Numerical zero.
#'
#' @return Matrix of deflated relative occurrences for all variables (i.e., j = 1,..., p) and for T = 1, ..., T_stop.
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- c(Gauss_data$y)
#' rand_exp <- random_experiments(X, y)
#' Phi_prime_fun(
#'   p = ncol(X),
#'   T_stop = rand_exp$T_stop,
#'   num_dummies = rand_exp$num_dummies,
#'   phi_T_mat = rand_exp$phi_T_mat,
#'   Phi = rand_exp$Phi,
#'   eps = rand_exp$eps
#' )
Phi_prime_fun <- function(p,
                          T_stop,
                          num_dummies,
                          phi_T_mat,
                          Phi,
                          eps = .Machine$double.eps) {
  av_num_var_sel <- colSums(phi_T_mat)
  fifty_phi_T_mat <- phi_T_mat[Phi > 0.5, , drop = FALSE]
  delta_av_num_var_sel <- colSums(fifty_phi_T_mat)

  if (T_stop > 1) {
    delta_av_num_var_sel[2:T_stop] <- delta_av_num_var_sel[2:T_stop] - delta_av_num_var_sel[1:(T_stop - 1)]
    phi_T_mat[, 2:T_stop] <- phi_T_mat[, 2:T_stop] - phi_T_mat[, 1:(T_stop - 1)]
  }

  phi_scale <- rep(NA, times = length(delta_av_num_var_sel))
  for (t in seq_along(delta_av_num_var_sel)) {
    if (delta_av_num_var_sel[t] > eps) {
      phi_scale[t] <- 1 - (((p - av_num_var_sel[t]) / (num_dummies - t + 1)) / delta_av_num_var_sel[t])
    } else {
      phi_scale[t] <- 0
    }
  }

  Phi_prime <- phi_T_mat %*% phi_scale

  return(Phi_prime)
}
