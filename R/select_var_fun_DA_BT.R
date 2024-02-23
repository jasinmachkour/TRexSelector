#' Compute set of selected variables for the T-Rex+DA+BT selector T-Rex+DA+BT selector (\doi{10.48550/arXiv.2401.15796})
#'
#' Computes the set of selected variables and returns the estimated support vector for the T-Rex+DA+BT selector (\doi{10.48550/arXiv.2401.15796}).
#'
#' @param p Number of candidate variables.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param FDP_hat_array_BT Array containing the conservative FDP estimates for all variables (dimension 1),
#' values of the voting level grid (dimension 2), and values of the dendrogram grid (dimension 3).
#' @param Phi_array_BT Array of relative occurrences as determined by the T-Rex calibration algorithm.
#' @param V Voting level grid.
#' @param rho_grid Dendrogram grid.
#'
#' @return List containing the estimated support vector, etc.
select_var_fun_DA_BT <- function(p,
                                 tFDR,
                                 T_stop,
                                 FDP_hat_array_BT,
                                 Phi_array_BT,
                                 V,
                                 rho_grid) {
  # Remove last row in FDP_hat_array_BT and Phi_array_BT if T_stop > 1
  if (T_stop > 1) {
    FDP_hat_array_BT <- array(FDP_hat_array_BT[-T_stop, , ], dim = dim(FDP_hat_array_BT) - c(1, 0, 0))
    Phi_array_BT <- array(Phi_array_BT[-T_stop, , ], dim = dim(Phi_array_BT) - c(1, 0, 0))
  }

  # Generate R_array
  R_array <- array(NA, dim = dim(FDP_hat_array_BT))
  for (TT in seq(dim(FDP_hat_array_BT)[1])) {
    for (VV in seq_along(V)) {
      for (RR in seq(dim(FDP_hat_array_BT)[3])) {
        R_array[TT, VV, RR] <- sum(Phi_array_BT[TT, , RR] > V[VV])
      }
    }
  }

  # T-Rex+DA+BT: Select variables
  FDP_hat_array_BT[FDP_hat_array_BT > tFDR] <- Inf
  val_max <- suppressWarnings(max(R_array[!is.infinite(FDP_hat_array_BT)]))
  ind_max <- matrix(which(R_array == val_max, arr.ind = TRUE), ncol = 3)
  ind_thresh <- matrix(ind_max[which(ind_max[, 2] == max(ind_max[, 2])), ], ncol = 3)
  ind_thresh <- ind_thresh[nrow(ind_thresh), ]
  thresh_T_dummy <- V[ind_thresh[2]]
  rho_thresh <- rho_grid[ind_thresh[3]]
  selected_var_T_dummy <- which(Phi_array_BT[ind_thresh[1], , ind_thresh[3]] > thresh_T_dummy)
  selected_var <- rep(0, times = p)
  selected_var[selected_var_T_dummy] <- 1

  res_select_var <- list(
    selected_var = selected_var,
    v_thresh = thresh_T_dummy,
    rho_thresh = rho_thresh,
    R_array = R_array
  )
  return(res_select_var)
}
