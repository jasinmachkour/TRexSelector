#' Runs the T-Knock filter with the extended calibration algorithm
#'
#' The T-Knock filter performs high-dimensional variable selection for a user-defined target false discovery rate (target FDR).
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param K Number of random experiments.
#' @param max_LL Integer factor determining the maximum number of dummies as a multiple of p (i.e., num_dummies = max_LL * p).
#' @param max_T_stop If TRUE the maximum number of dummies that can be included before stopping is set to ceiling(n / 2).
#' @param method 'tknock' for T-Knock filter and 'tknock+GVS' for T-Knock+GVS filter.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param parallel_process If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing (default: number of physical cores).
#' @param seed Seed for random number generator (only if parallel_process = TRUE).
#' @param eps Numerical zero.
#' @param verbose If TRUE progress in computations is shown.
#'
#' @return A list containing the support vector, and the values of T and L as determined by the extended calibration algorithm.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X = Gauss_data$X
#' y = c(Gauss_data$y)
#' set.seed(1234)
#' res = tknock(X, y)
#' res$beta.selected
tknock = function(X,
                  y,
                  tFDR = 0.2,
                  K = 20,
                  max_LL = 10,
                  max_T_stop = TRUE,
                  method = 'tknock',
                  type = 'lar',
                  corr_max = 0.5,
                  lambda_2_lars = NULL,
                  parallel_process = FALSE,
                  parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                  seed = NULL,
                  eps = .Machine$double.eps,
                  verbose = TRUE) {
  # Scale X and center y
  X = scale(X)
  y = y - mean(y)

  # Voting level grid
  V = seq(0.5, 1 - eps, by = 1 / K)
  V_len = length(V)

  # Initialize L-loop
  n = nrow(X)
  p = ncol(X)

  LL = 1
  T_stop = 1
  FDP_hat = rep(NA, times = V_len)

  # Setup cluster
  if (parallel_process && foreach::getDoParWorkers() == 1) {
    cl = parallel::makeCluster(parallel_max_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  while ((LL <= max_LL &&
          FDP_hat[which(abs(V - 0.75) < eps)] > tFDR) ||
         sum(!is.na(FDP_hat)) == 0) {
    num_dummies = LL * p
    LL = LL + 1

    # K Random experiments
    rand_exp = random_experiments(
      X = X,
      y = y,
      K = K,
      T_stop = T_stop,
      num_dummies = num_dummies,
      method = method,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      earlyStop = TRUE,
      verbose = verbose,
      intercept = FALSE,
      normalize = TRUE,
      cor.structure = FALSE,
      empirical = FALSE,
      parallel_process = parallel_process,
      parallel_max_cores = parallel_max_cores,
      seed = seed,
      eps = eps
    )
    phi_T_mat = rand_exp$phi_T_mat
    rep_osp.mat = rand_exp$rep_osp.mat
    Phi = rand_exp$Phi

    # Phi_prime
    Phi_prime = Phi_prime_fun(
      p = p,
      T_stop = T_stop,
      num_dummies = num_dummies,
      phi_T_mat = phi_T_mat,
      Phi = Phi,
      eps = eps
    )

    FDP_hat = fdp_hat(
      V = V,
      Phi = Phi,
      Phi_prime = Phi_prime,
      T_stop = T_stop,
      num_dummies = num_dummies
    )

    if (verbose) {
      cat(paste('\n L =', num_dummies, '\n'))
    }
  }

  # Initialize T-loop
  FDP_hat_mat = matrix(FDP_hat, nrow = 1)
  Phi_mat = matrix(Phi, nrow = 1)

  if (max_T_stop) {
    max_T = ceiling(n / 2)
  } else{
    max_T = num_dummies
  }

  # Reset seed
  if(!is.null(seed)){
    seed = seed + 12345
  }

  while (FDP_hat[V_len] <= tFDR && T_stop < max_T) {
    T_stop = T_stop + 1

    # K Random experiments
    rand_exp = random_experiments(
      X = X,
      y = y,
      K = K,
      T_stop = T_stop,
      num_dummies = num_dummies,
      method = method,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      earlyStop = TRUE,
      lars_state_list = rand_exp$lars_state_list,
      verbose = verbose,
      intercept = FALSE,
      normalize = TRUE,
      cor.structure = FALSE,
      empirical = FALSE,
      parallel_process = parallel_process,
      parallel_max_cores = parallel_max_cores,
      seed = seed,
      eps = eps
    )
    phi_T_mat = rand_exp$phi_T_mat
    rep_osp.mat = rand_exp$rep_osp.mat
    Phi = rand_exp$Phi
    lars_state_list = rand_exp$lars_state_list

    Phi_mat = rbind(Phi_mat, Phi)

    # Phi_prime
    Phi_prime = Phi_prime_fun(
      p = p,
      T_stop = T_stop,
      num_dummies = num_dummies,
      phi_T_mat = phi_T_mat,
      Phi = Phi,
      eps = eps
    )

    FDP_hat = fdp_hat(
      V = V,
      Phi = Phi,
      Phi_prime = Phi_prime,
      T_stop = T_stop,
      num_dummies = num_dummies
    )
    FDP_hat_mat = rbind(FDP_hat_mat, FDP_hat)

    if (verbose) {
      cat(paste('\n T =' , T_stop, '\n'))
    }
  }

  # T-Knock: Select variables
  res_T_dummy = select_var_fun(
    p = p,
    tFDR = tFDR,
    T_stop = T_stop,
    FDP_hat_mat = FDP_hat_mat,
    Phi_mat = Phi_mat,
    V = V
  )
  beta.selected = res_T_dummy$beta.selected
  v_thresh = res_T_dummy$v_thresh
  R_mat = res_T_dummy$R_mat

  res = list(
    beta.selected = beta.selected,
    tFDR = tFDR,
    T_stop = T_stop,
    num_dummies = num_dummies,
    V = V,
    v_thresh = v_thresh,
    FDP_hat_mat = FDP_hat_mat,
    Phi_mat = Phi_mat,
    R_mat = R_mat,
    phi_T_mat = phi_T_mat,
    Phi_prime = Phi_prime,
    method = method,
    type = type,
    corr_max = corr_max,
    lambda_2_lars = lambda_2_lars
  )

  return(res)
}
