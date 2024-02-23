#' Compute set of selected variables
#'
#' Computes the set of selected variables and returns the estimated support vector for the T-Rex selector (\doi{10.48550/arXiv.2110.06048}).
#'
#' @param p Number of candidate variables.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param FDP_hat_mat Matrix whose rows are the vectors of conservative FDP estimates for each value of the voting level grid.
#' @param Phi_mat Matrix of relative occurrences as determined by the T-Rex calibration algorithm.
#' @param V Voting level grid.
#'
#' @return Estimated support vector.
select_var_fun <- function(p,
                           tFDR,
                           T_stop,
                           FDP_hat_mat,
                           Phi_mat,
                           V) {
  # Remove last row in FDP_hat_mat and Phi_mat if T_stop > 1
  if (T_stop > 1) {
    FDP_hat_mat <- FDP_hat_mat[-T_stop, , drop = FALSE]
    Phi_mat <- Phi_mat[-T_stop, , drop = FALSE]
  }

  # Generate R_mat
  R_mat <- matrix(NA,
    nrow = nrow(FDP_hat_mat),
    ncol = ncol(FDP_hat_mat)
  )
  for (TT in seq(nrow(FDP_hat_mat))) {
    for (VV in seq_along(V)) {
      R_mat[TT, VV] <- sum(Phi_mat[TT, ] > V[VV])
    }
  }

  # T-Rex: Select variables
  FDP_hat_mat[FDP_hat_mat > tFDR] <- Inf
  val_max <- suppressWarnings(max(R_mat[!is.infinite(FDP_hat_mat)]))
  ind_max <- matrix(which(R_mat == val_max, arr.ind = TRUE), ncol = 2)
  ind_thresh <- ind_max[nrow(ind_max), ]
  thres_T_dummy <- V[ind_thresh[2]]
  selected_var_T_dummy <- which(Phi_mat[ind_thresh[1], ] > thres_T_dummy)
  selected_var <- rep(0, times = p)
  selected_var[selected_var_T_dummy] <- 1

  res_select_var <- list(
    selected_var = selected_var,
    v_thresh = thres_T_dummy,
    R_mat = R_mat
  )
  return(res_select_var)
}
