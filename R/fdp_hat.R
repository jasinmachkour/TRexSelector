#' Computes the conservative FDP estimate of the T-Rex selector (\doi{10.48550/arXiv.2110.06048})
#'
#' @param V Voting level grid.
#' @param Phi Vector of relative occurrences.
#' @param Phi_prime Vector of deflated relative occurrences.
#' @param eps Numerical zero.
#'
#' @return Vector of conservative FDP estimates for each value of the voting level grid.
fdp_hat <- function(V,
                    Phi,
                    Phi_prime,
                    eps = .Machine$double.eps) {
  fdp_h <- rep(NA, times = length(V))
  for (i in seq_along(V)) {
    num_sel_var <- sum(Phi > V[i])
    if (num_sel_var < eps) {
      fdp_h[i] <- 0
    } else {
      fdp_h[i] <- min(1, (sum((1 - Phi_prime)[Phi > V[i]]) / (num_sel_var)))
    }
  }
  return(fdp_h)
}
