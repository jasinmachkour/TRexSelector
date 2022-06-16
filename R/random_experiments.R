#' Run K random experiments
#'
#' Run K random experiments and compute the matrix of relative occurrences for all variables
#' and all numbers of included variables before stopping.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param K Number of random experiments.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#' @param method 'tknock' for the T-Knock filter and 'tknock+GVS' for the T-Knock+GVS filter.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise
#' the entire solution path is computed.
#' @param lars_state_list If parallel_process = TRUE: List of state variables of the previous T-LARS steps of the K random experiments
#' (necessary for warm-starts, i.e., restarting the forward selection process exactly where it was previously terminated).
#' If parallel_process = FALSE: List of objects of the class tlars_cpp associated with the K random experiments
#' (necessary for warm-starts, i.e., restarting the forward selection process exactly where it was previously terminated).
#' @param verbose Logical. If TRUE progress in computations is shown.
#' @param intercept Logical. If TRUE an intercept is included.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#' @param parallel_process Logical. If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing
#' (default: minimum{Number of random experiments K, number of physical cores}).
#' @param seed Seed for random number generator (ignored if parallel_process = FALSE).
#' @param eps Numerical zero.
#'
#' @return List containing the results of the K random experiments.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParWorkers registerDoSEQ `%do%` `%dopar%` foreach
#' @importFrom doRNG `%dorng%`
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- c(Gauss_data$y)
#' res <- random_experiments(X = X, y = y)
#' relative_occurrences_matrix <- res$phi_T_mat
#' relative_occurrences_matrix
random_experiments <- function(X,
                               y,
                               K = 20,
                               T_stop = 1,
                               num_dummies = ncol(X),
                               method = "tknock",
                               type = "lar",
                               corr_max = 0.5,
                               lambda_2_lars = NULL,
                               early_stop = TRUE,
                               lars_state_list,
                               verbose = TRUE,
                               intercept = FALSE,
                               standardize = TRUE,
                               parallel_process = FALSE,
                               parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                               seed = NULL,
                               eps = .Machine$double.eps) {
  # Error control
  method <- match.arg(method, c("tknock", "tknock+GVS"))

  type <- match.arg(type, c("lar", "lasso"))

  if (!is.matrix(X)) {
    stop("'X' must be a matrix.")
  }

  if (!is.numeric(X)) {
    stop("'X' only allows numerical values.")
  }

  if (anyNA(X)) {
    stop("'X' contains NAs. Please remove or impute them before proceeding.")
  }

  if (!is.vector(drop(y))) {
    stop("'y' must be a vector.")
  }

  if (!is.numeric(y)) {
    stop("'y' only allows numerical values.")
  }

  if (anyNA(y)) {
    stop("'y' contains NAs. Please remove or impute them before proceeding.")
  }

  if (nrow(X) != length(drop(y))) {
    stop("Number of rows in X does not match length of y.")
  }

  if (length(K) != 1 ||
    K < 2 ||
    K %% 1 != 0) {
    stop("The number of random experiments 'K' must be an integer larger or equal to 2.")
  }

  if (method == "tknock") {
    if (length(num_dummies) != 1 ||
      num_dummies %% 1 != 0 ||
      num_dummies < 1) {
      stop("'num_dummies' must be an integer larger or equal to 1.")
    }
  }

  # Number of variables in X
  p <- ncol(X)

  # Continue error control
  if (method == "tknock+GVS") {
    if (length(num_dummies) != 1 ||
      num_dummies %% p != 0 ||
      num_dummies < 1) {
      stop(
        "'num_dummies' must be a positive integer multiple of the total number of original predictors in X."
      )
    }
  }

  if (length(T_stop) != 1 ||
    !(T_stop %in% seq(1, num_dummies))) {
    stop(
      paste0(
        "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
        num_dummies,
        "."
      )
    )
  }

  if (method == "tknock+GVS") {
    if (length(corr_max) != 1 ||
      corr_max < 0 ||
      corr_max > 1) {
      stop("'corr_max' must have a value between zero and one.")
    }

    if (!is.null(lambda_2_lars)) {
      if (length(lambda_2_lars) != 1 ||
        lambda_2_lars < eps) {
        stop("'lambda_2_lars' must be a number larger than zero.")
      }
    }
  }

  if (parallel_process &&
    (
      length(parallel_max_cores) != 1 ||
        parallel_max_cores %% 1 != 0 ||
        parallel_max_cores < 2
    )) {
    stop(
      "For parallel processing at least two workers have to be registered:
         'parallel_max_cores' must be an integer larger or equal to 2."
    )
  }

  if (parallel_process &&
    parallel_max_cores > min(K, max(1, parallel::detectCores(logical = FALSE)))) {
    parallel_max_cores_modified <-
      min(K, max(1, parallel::detectCores(logical = FALSE)))
    message(
      paste0(
        "For computing ",
        K,
        " random experiments, it is not useful to register ",
        parallel_max_cores,
        " workers. Setting parallel_max_cores = ",
        min(K, max(
          1, parallel::detectCores(logical = FALSE)
        )),
        " (# physical cores) ...\n"
      )
    )
    parallel_max_cores <-
      min(K, max(1, parallel::detectCores(logical = FALSE)))
  }

  if (parallel_process && T_stop == 1 && num_dummies <= p) {
    message(
      "Computing random experiments in parallel...
            Note that this is only advantageous if you have at least a few thousand predictors and/or data points in 'X'.
            Otherwise, the overhead will slow down the computations in parallel. Thus, for small data sizes it is better
            to set parallel_process = FALSE.

      Be careful!"
    )
  }

  if (!(missing(lars_state_list) || is.null(lars_state_list))) {
    if (length(lars_state_list) != K) {
      stop("Length of 'lars_state_list' must be equal to number of random experiments 'K'.")
    }
  }

  # Create empty lars_state_list if missing or NULL
  if (missing(lars_state_list) || is.null(lars_state_list)) {
    lars_state_list <- vector(mode = "list", length = K)
  }

  # Combines Output Lists of Parallel 'foreach' Loop
  comb_fun <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) {
        c(x[[i]], lapply(list(...), function(y) {
          y[[i]]
        }))
      }
    )
  }

  # Setup cluster
  if (parallel_process && foreach::getDoParWorkers() == 1) {
    cl <- parallel::makeCluster(parallel_max_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  `%par_exe%` <-
    ifelse(parallel_process, doRNG::`%dorng%`, foreach::`%do%`)
  h <- NULL
  res <- foreach::foreach(
    h = seq(K),
    .combine = comb_fun,
    .multicombine = TRUE,
    .init = list(list(), list(), list()),
    .options.RNG = seed
  ) %par_exe% {
    # Recreate tlarsCpp object if necessary
    if (parallel_process &&
      !is.null(lars_state_list[[h]]) &&
      !methods::is(object = lars_state_list[[h]], class2 = tlars::tlars_cpp)) {
      lars_state <- tlars::tlars_model(lars_state = lars_state_list[[h]])
    } else {
      lars_state <- lars_state_list[[h]]
    }

    # Run random experiment
    lars_state <- lm_dummy(
      X = X,
      y = y,
      model_tlars = lars_state,
      T_stop = T_stop,
      num_dummies = num_dummies,
      method = method,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      early_stop = early_stop,
      verbose = verbose,
      intercept = intercept,
      standardize = standardize
    )

    # Extract T-LARS path
    lars_path <- do.call(cbind, lars_state$get_beta_path())

    # Extract content of object lars_state if performing parallel computations
    if (parallel_process) {
      lars_state <- lars_state$get_all()
    }

    # Number of included dummies along solution path
    dummy_num_path <- colSums(matrix(
      abs(lars_path[(p + 1):(p + num_dummies), ]) > eps,
      nrow = num_dummies,
      ncol = ncol(lars_path)
    ))

    # Number of included original variables along solution path
    var_num_path <- colSums(matrix(
      abs(lars_path[1:p, ]) > eps,
      nrow = p,
      ncol = ncol(lars_path)
    ))

    # Relative occurrences
    phi_T_mat <- matrix(0, nrow = p, ncol = T_stop)
    for (c in seq(T_stop)) {
      if (!any(dummy_num_path == c)) {
        ind_sol_path <- length(dummy_num_path)
        warning(
          paste(
            "For T_stop = ",
            c,
            " LARS is running until k = min(n, p) and stops there before selecting ",
            c,
            " dummies.",
            sep = ""
          )
        )
      } else {
        ind_sol_path <- which(as.numeric(dummy_num_path) == c)[1]
      }
      phi_T_mat[, c] <-
        (1 / K) * (abs(lars_path[1:p, ind_sol_path]) > eps)
    }

    rand_exp_last_betas_mat <- lars_path[1:p, ncol(lars_path)]

    list(
      phi_T_mat,
      rand_exp_last_betas_mat,
      lars_state
    )
  }

  # Merging results of all random experiments
  lars_state_list <- res[[3]]
  names(lars_state_list) <-
    paste("lars_state (K = ", seq(K), ")", sep = "")

  phi_T_mat <- Reduce("+", res[[1]])
  rand_exp_last_betas_mat <- unname(Reduce(rbind, res[[2]]))
  Phi <- apply(abs(rand_exp_last_betas_mat) > eps, 2, sum) / K

  # List of results
  rand_exp_res <- list(
    phi_T_mat = phi_T_mat,
    rand_exp_last_betas_mat = rand_exp_last_betas_mat,
    Phi = Phi,
    lars_state_list = lars_state_list,
    K = K,
    T_stop = T_stop,
    num_dummies = num_dummies,
    method = method,
    type = type,
    corr_max = corr_max,
    lambda_2_lars = lambda_2_lars,
    seed = seed,
    eps = eps
  )

  return(rand_exp_res)
}
