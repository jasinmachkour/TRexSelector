#' Perform one random experiment
#'
#' Run one random experiment of the T-Knock filter, i.e., generates dummies, appends them to the predictor matrix, and runs the forward selection algorithm until it is terminated after T_stop dummies have been selected.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param num_dummies Number of dummies
#' @param method 'tknock' for T-Knock filter and 'tknock+GVS' for T-Knock+GVS filter.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param early_stop If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise the entire solution path is computed.
#' @param lars_state LARS object (i.e., object of the class tlarsCpp). Contains variables associated with previous stopping point (necessary to restart forward selection selection process exactly where it was previously terminated).
#' @param verbose If TRUE progress in computations is shown.
#' @param intercept If TRUE an intercept is included.
#' @param normalize If TRUE the predictors are standardized and the response is centered.
#' @param cor.structure TRUE/FALSE.
#' @param empirical TRUE/FALSE.
#'
#' @return LARS object (i.e., object of the class tlars_cpp).
#'
#' @import tlars
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
#' res <- lm_dummy(X, y, T_stop = 1, num_dummies = 5 * p)
#' support <- abs(res$get_beta()[seq(p)]) > eps
#' support
lm_dummy <- function(X,
                     y,
                     T_stop = 1,
                     num_dummies = ncol(X),
                     method = "tknock",
                     type = "lar",
                     corr_max = 0.5,
                     lambda_2_lars = NULL,
                     early_stop = TRUE,
                     lars_state,
                     verbose = TRUE,
                     intercept = FALSE,
                     normalize = TRUE,
                     cor.structure = FALSE,
                     empirical = FALSE) {
  method <- match.arg(method, c("tknock", "tknock+GVS"))

  if (T_stop == 1 || missing(lars_state) || is.null(lars_state)) {
    if (method == "tknock") {
      X_dummy <- add_dummies(
        X,
        num_dummies = num_dummies,
        cor.structure = cor.structure,
        empirical = empirical
      )
    } else {
      X_dummy <- add_dummies_GVS(X,
        num_dummies = num_dummies,
        corr_max = corr_max
      )
      # Data modification for Elastic Net
      p_dummy <- ncol(X_dummy)
      X_dummy <- (1 / sqrt(1 + lambda_2_lars)) * rbind(X_dummy, diag(rep(sqrt(
        lambda_2_lars
      ), times = p_dummy)))
      y <- append(y, rep(0, times = p_dummy))
      # Scale data again
      X_dummy <- scale(X_dummy)
      y <- y - mean(y)
    }
    # Create new LARS object
    lars_state <- tlars::tlars_model(
      X = X_dummy,
      y = y,
      num_dummies = num_dummies,
      verbose = FALSE,
      info = FALSE,
      intercept = intercept,
      standardize = normalize,
      type = type
    )
  }

  # Execute LARS step
  tlars::tlars(
    model = lars_state,
    T_stop = T_stop,
    early_stop = early_stop,
    info = FALSE
  )

  return(lars_state)
}
