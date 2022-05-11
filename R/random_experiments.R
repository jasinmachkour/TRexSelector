#' Runs K random experiments and computes the matrix of relative occurrences
#'
#' @param X Real valued Predictor matrix.
#' @param y Response vector.
#' @param K Number of random experiments.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param L_val Number of dummies
#' @param method 'tknock' for T-Knock filter and 'tknock+GVS' for T-Knock+GVS filter.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param earlyStop If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise the entire solution path is computed.
#' @param lars_state_list List of variables associated with previous stopping points of the K random experiments (necessary to restart forward selection selection process exactly where it was previously terminated).
#' @param verbose If TRUE progress in computations is shown.
#' @param intercept If TRUE an intercept is included.
#' @param normalize If TRUE the predictors are standardized and the response is centered.
#' @param cor.structure TRUE/FALSE.
#' @param empirical TRUE/FALSE.
#' @param parallel_process If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing (default: number of physical cores).
#' @param seed Seed for random number generator (only if parallel_process = TRUE).
#' @param eps Numerical zero.
#'
#' @return Results of the K random experiments.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import doRNG
#' @import methods
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data('Gauss_data')
#' X = Gauss_data$X
#' y = c(Gauss_data$y)
#' res = random_experiments(X, y)
#' res$phi_T_mat
random_experiments = function(X,
                              y,
                              K = 20,
                              T_stop = 1,
                              L_val = ncol(X),
                              method = 'tknock',
                              type = 'lar',
                              corr_max = 0.5,
                              lambda_2_lars = NULL,
                              earlyStop = TRUE,
                              lars_state_list,
                              verbose = TRUE,
                              intercept = FALSE,
                              normalize = TRUE,
                              cor.structure = FALSE,
                              empirical = TRUE,
                              parallel_process = FALSE,
                              parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                              seed = NULL,
                              eps = .Machine$double.eps) {
  p = ncol(X)

  if (missing(lars_state_list) || is.null(lars_state_list)) {
    lars_state_list = vector(mode = 'list', length = K)
  }

  # Combines Output Lists of Parallel 'foreach' Loop
  comb_fun <- function(x, ...) {
    lapply(seq_along(x),
           function(i)
             c(x[[i]], lapply(list(...), function(y)
               y[[i]])))
  }

  # Setup cluster
  if (parallel_process && foreach::getDoParWorkers() == 1) {
    cl = parallel::makeCluster(parallel_max_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  `%par_exe%` = ifelse(parallel_process, doRNG::`%dorng%`, foreach::`%do%`)
  h = NULL
  res = foreach::foreach (
    h = seq(K),
    .combine = comb_fun,
    .multicombine = TRUE,
    .init = list(list(), list(), list()),
    .options.RNG = seed
  ) %par_exe% {

    # Recreate tlarsCpp object
    if (parallel_process && !is.null(lars_state_list[[h]]) && !methods::is(object = lars_state_list[[h]], class2 = tlars::tlars_cpp)) {
      lars_state = tlars::tlars_model(lars_state = lars_state_list[[h]])
    }else{
      lars_state = lars_state_list[[h]]
    }

    lars_state = lm_dummy(
      X = X,
      y = y,
      T_stop = T_stop,
      num_dummies = L_val,
      method = method,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      early_stop = earlyStop,
      lars_state = lars_state,
      verbose = verbose,
      intercept = intercept,
      normalize = normalize,
      cor.structure = cor.structure,
      empirical = empirical
    )

    rep_osp_dummy.selected = do.call(cbind, lars_state$get_beta_path())

    # Extract content of object lars_state if performing parallel computations
    if (parallel_process) {
      lars_state = lars_state$get_all()
    }

    dummy_num_path = colSums(matrix(
      abs(rep_osp_dummy.selected[(p + 1):(p + L_val), ]) > eps,
      nrow = L_val,
      ncol = ncol(rep_osp_dummy.selected)
    ))
    var_num_path = colSums(matrix(
      abs(rep_osp_dummy.selected[1:p, ]) > eps,
      nrow = p,
      ncol = ncol(rep_osp_dummy.selected)
    ))

    phi_T_mat = matrix(0, nrow = p, ncol = T_stop)
    for (c in seq(T_stop)) {
      if (!any(dummy_num_path == c)) {
        ind_sol_path = length(dummy_num_path)
        warning(
          paste(
            'For T_stop = ',
            c,
            ' LARS is running until k = min(n, p) and stops there before selecting ',
            c,
            ' dummy/dummies.',
            sep = ''
          )
        )
      } else{
        ind_sol_path = which(as.numeric(dummy_num_path) == c)[1]
      }
      phi_T_mat[, c] = (1 / K) * (abs(rep_osp_dummy.selected[1:p, ind_sol_path]) > eps)
    }

    rep_osp.mat = rep_osp_dummy.selected[1:p, ncol(rep_osp_dummy.selected)]

    list(phi_T_mat,
         rep_osp.mat,
         lars_state)
  }

  lars_state_list = res[[3]]
  names(lars_state_list) = paste('lars_state (K = ', seq(K) , ')', sep = '')

  phi_T_mat = Reduce('+', res[[1]])
  rep_osp.mat = unname(Reduce(rbind, res[[2]]))
  Phi = apply(abs(rep_osp.mat) > eps, 2, sum) / K

  rand_exp_res = list(
    phi_T_mat = phi_T_mat,
    rep_osp.mat = rep_osp.mat,
    Phi = Phi,
    lars_state_list = lars_state_list,
    K = K,
    T_stop = T_stop,
    L_val = L_val,
    method = method,
    type = type,
    corr_max = corr_max,
    lambda_2_lars = lambda_2_lars,
    seed = seed,
    eps = eps
  )

  return(rand_exp_res)
}
