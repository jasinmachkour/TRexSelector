#' Perform one random experiment
#'
#' Run one random experiment of the T-Rex selector (\doi{10.48550/arXiv.2110.06048}), i.e., generates dummies, appends them to the predictor matrix, and runs
#' the forward selection algorithm until it is terminated after T_stop dummies have been selected.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param model_tlars Object of the class tlars_cpp. It contains all state variables of the previous T-LARS step (necessary for warm-starts, i.e., restarting
#' the forward selection process exactly where it was previously terminated).
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#' @param method 'trex' for the T-Rex selector (\doi{10.48550/arXiv.2110.06048}),
#' 'trex+GVS' for the T-Rex+GVS selector (\doi{10.23919/EUSIPCO55093.2022.9909883}),
#' 'trex+DA+AR1' for the T-Rex+DA+AR1 selector,
#' 'trex+DA+equi' for the T-Rex+DA+equi selector,
#' 'trex+DA+BT' for the T-Rex+DA+BT selector (\doi{10.48550/arXiv.2401.15796}),
#' 'trex+DA+NN' for the T-Rex+DA+NN selector (\doi{10.48550/arXiv.2401.15139}).
#' @param GVS_type 'IEN' for the Informed Elastic Net (\doi{10.1109/CAMSAP58249.2023.10403489}),
#' 'EN' for the ordinary Elastic Net (\doi{10.1111/j.1467-9868.2005.00503.x}).
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise
#' the entire solution path is computed.
#' @param verbose Logical. If TRUE progress in computations is shown when performing T-LARS steps on the created model.
#' @param intercept Logical. If TRUE an intercept is included.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#'
#' @return Object of the class tlars_cpp.
#'
#' @importFrom tlars tlars_model tlars tlars_cpp
#' @importFrom glmnet cv.glmnet
#' @importFrom stats rnorm
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' eps <- .Machine$double.eps
#' n <- 75
#' p <- 100
#' X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' beta <- c(rep(3, times = 3), rep(0, times = 97))
#' y <- X %*% beta + rnorm(n)
#' res <- lm_dummy(X = X, y = y, T_stop = 1, num_dummies = 5 * p)
#' beta_hat <- res$get_beta()[seq(p)]
#' support <- abs(beta_hat) > eps
#' support
lm_dummy <- function(X,
                     y,
                     model_tlars,
                     T_stop = 1,
                     num_dummies = ncol(X),
                     method = "trex",
                     GVS_type = "IEN",
                     type = "lar",
                     corr_max = 0.5,
                     lambda_2_lars = NULL,
                     early_stop = TRUE,
                     verbose = TRUE,
                     intercept = FALSE,
                     standardize = TRUE) {
  # Numerical zero
  eps <- .Machine$double.eps

  # Error control
  method <- match.arg(method, c("trex", "trex+GVS", "trex+DA+AR1", "trex+DA+equi", "trex+DA+BT", "trex+DA+NN"))

  type <- match.arg(type, c("lar", "lasso"))

  GVS_type <- match.arg(GVS_type, c("IEN", "EN"))

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

  if (!(missing(model_tlars) || is.null(model_tlars))) {
    if (!methods::is(object = model_tlars, class2 = tlars::tlars_cpp)) {
      stop("'model_tlars' must be an object of class tlars_cpp.")
    }
  }

  if (method == "trex" || method == "trex+DA+AR1" || method == "trex+DA+equi" || method == "trex+DA+BT" || method == "trex+DA+NN") {
    if (length(num_dummies) != 1 ||
      num_dummies %% 1 != 0 ||
      num_dummies < 1) {
      stop("'num_dummies' must be an integer larger or equal to 1.")
    }
  }

  if (method == "trex+GVS") {
    if (length(num_dummies) != 1 ||
      num_dummies %% ncol(X) != 0 ||
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

  if (method == "trex+GVS") {
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

  # Create T-LARS model if it is not supplied (from a previous T-LARS step)
  if (T_stop == 1 ||
    missing(model_tlars) ||
    is.null(model_tlars)) {
    if (method == "trex" || method == "trex+DA+AR1" || method == "trex+DA+equi" || method == "trex+DA+BT" || method == "trex+DA+NN") {
      X_Dummy <- add_dummies(
        X = X,
        num_dummies = num_dummies
      )
    } else {
      GVS_dummies <- add_dummies_GVS(
        X = X,
        num_dummies = num_dummies,
        corr_max = corr_max
      )
      X_Dummy <- GVS_dummies$X_Dummy

      # Ridge regression to determine lambda_2 for elastic net
      if (is.null(lambda_2_lars)) {
        n <- ncol(X)
        alpha <- 0
        cvfit <-
          glmnet::cv.glmnet(
            x = X,
            y = y,
            intercept = intercept,
            standardize = standardize,
            alpha = alpha,
            type.measure = "mse",
            family = "gaussian",
            nfolds = 10
          )
        lambda_2_glmnet <- cvfit$lambda.1se
        lambda_2_lars <- lambda_2_glmnet * n * (1 - alpha) / 2
      }

      # Data modification for Elastic Net (EN)
      if (GVS_type == "EN") {
        p_dummy <- ncol(X_Dummy)
        X_Dummy <-
          (1 / sqrt(1 + lambda_2_lars)) * rbind(X_Dummy, diag(rep(sqrt(
            lambda_2_lars
          ), times = p_dummy)))
        y <- append(y, rep(0, times = p_dummy))
      }

      # Data modification for Informed Elastic Net (IEN)
      if (GVS_type == "IEN") {
        p <- ncol(X)
        p_dummy <- ncol(X_Dummy)
        max_clusters <- GVS_dummies$max_clusters
        cluster_sizes <- GVS_dummies$cluster_sizes
        IEN_cl_id_vectors <- GVS_dummies$IEN_cl_id_vectors
        X_Dummy <-
          sqrt(lambda_2_lars) * rbind((1 / sqrt(lambda_2_lars)) * X_Dummy,
                                      (1 / sqrt(cluster_sizes)) *
                                        matrix(rep(IEN_cl_id_vectors, times = p_dummy / p),
                                               ncol = ncol(IEN_cl_id_vectors) * p_dummy / p))
        y <- append(y, rep(0, times = max_clusters))
      }

      # Scale data again
      X_Dummy <- scale(X_Dummy)
      y <- y - mean(y)
    }

    # Create new LARS object
    model_tlars <- tlars::tlars_model(
      X = X_Dummy,
      y = y,
      num_dummies = num_dummies,
      verbose = FALSE,
      info = FALSE,
      intercept = intercept,
      standardize = standardize,
      type = type
    )
  }

  # Execute LARS step
  tlars::tlars(
    model = model_tlars,
    T_stop = T_stop,
    early_stop = early_stop,
    info = FALSE
  )

  return(model_tlars)
}
