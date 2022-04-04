#' T-Knock: Selected Variables
#'
#' Computes the set of selected variables and returns the support vector.
#'
#' @param p Number of candidate variables.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param T_stop Number of included knockoffs after which the random experiments (i.e., forward selection processes) are stopped.
#' @param FDP_hat_mat Matrix whose rows are the vectors of conservative FDP estimates for each value of the voting level grid.
#' @param Phi_mat Matrix of relative occurrences as determined by the T-Knock calibration algorithm.
#' @param V Voting level grid.
#'
#' @return Support Vector.
#'
#' @export
select_var_fun = function(p,
                          tFDR,
                          T_stop,
                          FDP_hat_mat,
                          Phi_mat,
                          V) {
  # Remove last row in FDP_hat_mat and Phi_mat if T_stop > 1
  if (T_stop > 1) {
    FDP_hat_mat = FDP_hat_mat[-T_stop, , drop = FALSE]
    Phi_mat = Phi_mat[-T_stop, , drop = FALSE]
  }

  # Generate R_mat
  R_mat = matrix(NA,
                 nrow = nrow(FDP_hat_mat),
                 ncol = ncol(FDP_hat_mat))
  for (TT in seq(nrow(FDP_hat_mat))) {
    for (VV in seq_along(V)) {
      R_mat[TT, VV] = sum(Phi_mat[TT,] > V[VV])
    }
  }

  # T-Knock: Select variables
  FDP_hat_mat[FDP_hat_mat > tFDR] = Inf
  val_max = max(R_mat[!is.infinite(FDP_hat_mat)])
  ind_max = matrix(which(R_mat == val_max, arr.ind = TRUE), ncol = 2)
  ind_thresh = ind_max[nrow(ind_max),]
  thres_T_Knock = V[ind_thresh[2]]
  selected_var_T_Knock = which(Phi_mat[ind_thresh[1],] > thres_T_Knock)
  beta.selected = rep(0, times = p)
  beta.selected[selected_var_T_Knock] = 1

  res_select_var = list(beta.selected = beta.selected,
                        v_thresh = thres_T_Knock,
                        R_mat = R_mat)
  return(res_select_var)
}
