#' Perform one random experiment
#'
#' Run one random experiment of the T-Knock filter, i.e., generates knockoffs, appends them to the predictor matrix, and runs the forward selection algorithm until it is terminated after T_stop knockoffs have been selected.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param T_stop Number of included knockoffs after which the random experiments (i.e., forward selection processes) are stopped.
#' @param L_val Number of knockoffs.
#' @param method 'tknock' for T-Knock filter and 'tknock+GVS' for T-Knock+GVS filter.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param earlyStop If TRUE, then the forward selection process is stopped after T_stop knockoffs have been included. Otherwise the entire solution path is computed.
#' @param lars_state LARS object (i.e., object of the class tlarsCpp). Contains variables associated with previous stopping point (necessary to restart forward selection selection process exactly where it was previously terminated).
#' @param verbose If TRUE progress in computations is shown.
#' @param intercept If TRUE an intercept is included.
#' @param normalize If TRUE the predictors are standardized and the response is centered.
#' @param cor.structure TRUE/FALSE.
#' @param empirical TRUE/FALSE.
#'
#' @return LARS object (i.e., object of the class tlarsCpp).
#'
#' @import tlars
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' eps = .Machine$double.eps
#' n = 75
#' p = 100
#' X = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' beta = c(rep(3, times = 3), rep(0, times = 97))
#' y = X %*% beta + rnorm(n)
#' res = lm_knockoff(X, y, T_stop = 1, L_val = 5 * p)
#' support = abs(res$getLastBeta()[seq(p)]) > eps
#' support
lm_knockoff = function(X,
                       y,
                       T_stop = 1,
                       L_val = ncol(X),
                       method = 'tknock',
                       type = 'lar',
                       corr_max = 0.5,
                       lambda_2_lars = NULL,
                       earlyStop = TRUE,
                       lars_state,
                       verbose = TRUE,
                       intercept = FALSE,
                       normalize = TRUE,
                       cor.structure = FALSE,
                       empirical = FALSE) {
  method = match.arg(method, c('tknock', 'tknock+GVS'))

  if (T_stop == 1 || missing(lars_state) || is.null(lars_state)) {
    if (method == 'tknock') {
      X_Knock = add_knockoffs(
        X,
        L_val = L_val,
        cor.structure = cor.structure,
        empirical = empirical
      )
    } else{
      X_Knock = add_knockoffs_GVS(X,
                                  L_val = L_val,
                                  corr_max = corr_max)
      # Data modification for Elastic Net
      p_Knock = ncol(X_Knock)
      X_Knock = (1 / sqrt(1 + lambda_2_lars)) * rbind(X_Knock, diag(rep(sqrt(
        lambda_2_lars
      ), times = p_Knock)))
      y = append(y, rep(0, times = p_Knock))
      # Scale data again
      X_Knock = scale(X_Knock)
      y = y - mean(y)
    }
    # Create new LARS object
    lars_state = tlars::tlars_model(
      X = X_Knock,
      y = y,
      num_dummies = L_val,
      verbose = FALSE,
      info = FALSE,
      intercept = intercept,
      standardize = normalize,
      type = type
    )
  }

  # Execute LARS step
  tlars::tlars(model = lars_state,
               T_stop = T_stop,
               early_stop = earlyStop,
               info = FALSE)

  return(lars_state)
}
