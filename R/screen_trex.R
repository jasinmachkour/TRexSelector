#' Run the Screen-T-Rex selector (\doi{10.1109/SSP53291.2023.10207957})
#'
#' The Screen-T-Rex selector (\doi{10.1109/SSP53291.2023.10207957}) performs very fast variable selection in high-dimensional settings while
#' informing the user about the automatically selected false discovery rate (FDR).
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param K Number of random experiments.
#' @param R Number of bootstrap resamples.
#' @param method 'trex' for the T-Rex selector (\doi{10.48550/arXiv.2110.06048}),
#' 'trex+GVS' for the T-Rex+GVS selector (\doi{10.23919/EUSIPCO55093.2022.9909883}),
#' 'trex+DA+AR1' for the T-Rex+DA+AR1 selector,
#' 'trex+DA+equi' for the T-Rex+DA+equi selector.
#' @param bootstrap Logical. If TRUE Screen-T-Rex is carried out with bootstrapping.
#' @param conf_level_grid Confidence level grid for the bootstrap confidence intervals.
#' @param cor_coef AR(1) autocorrelation coefficient for the T-Rex+DA+AR1 selector or equicorrelation coefficient for the T-Rex+DA+equi selector.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param rho_thr_DA Correlation threshold for the T-Rex+DA+AR1 selector and the T-Rex+DA+equi selector (i.e., method = 'trex+DA+AR1' or 'trex+DA+equi').
#' @param parallel_process Logical. If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing.
#' @param seed Seed for random number generator (ignored if parallel_process = FALSE).
#' @param eps Numerical zero.
#' @param verbose Logical. If TRUE progress in computations is shown.
#'
#' @return A list containing the estimated support vector, the automatically selected false discovery rate (FDR) and additional information.
#'
#' @importFrom parallel detectCores
#' @importFrom stats coef arima
#' @importFrom boot boot boot.ci
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- c(Gauss_data$y)
#' set.seed(123)
#' res <- screen_trex(X = X, y = y)
#' selected_var <- res$selected_var
#' selected_var
screen_trex <- function(X,
                        y,
                        K = 20,
                        R = 1000,
                        method = "trex",
                        bootstrap = FALSE,
                        conf_level_grid = seq(0, 1, by = 0.001),
                        cor_coef = NA,
                        type = "lar",
                        corr_max = 0.5,
                        lambda_2_lars = NULL,
                        rho_thr_DA = 0.02,
                        parallel_process = FALSE,
                        parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                        seed = NULL,
                        eps = .Machine$double.eps,
                        verbose = TRUE) {

  # Error control
  method <- match.arg(method, c("trex", "trex+GVS", "trex+DA+AR1", "trex+DA+equi"))

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

  if (bootstrap) {
    if(min(conf_level_grid) < 0 ||
       max(conf_level_grid) > 1) {
      stop("conf_level_grid must be a vector containing numbers between 0 and 1.")
    }
  }

  # Scale X and center y
  X <- scale(X)
  y <- y - mean(y)

  # Number of candidate variables
  p <- ncol(X)

  # Set T_stop = 1
  T_stop = 1

  # Run random experiments
  rand_exp <- random_experiments(X = X,
                                 y = y,
                                 K = K,
                                 T_stop  = T_stop,
                                 num_dummies = ncol(X),
                                 method = method,
                                 type = type,
                                 corr_max = corr_max,
                                 lambda_2_lars = lambda_2_lars,
                                 early_stop = TRUE,
                                 verbose = verbose,
                                 intercept = FALSE,
                                 standardize = TRUE,
                                 dummy_coef = TRUE,
                                 parallel_process = parallel_process,
                                 parallel_max_cores = parallel_max_cores,
                                 seed = seed,
                                 eps = eps)

  # Extract matrix containing the terminal candidate coefficient vectors of all K random experiments as rows
  beta_mat <- rand_exp$rand_exp_last_betas_mat

  # Extract matrix containing the terminal dummy coefficient vectors of all K random experiments as rows
  dummy_beta_mat <- rand_exp$dummy_rand_exp_last_betas_mat

  # Extract relative occurrences matrix and vector of relative occurrences at T_stop = 1
  phi_T_mat <- rand_exp$phi_T_mat
  Phi <- rand_exp$Phi

  # Compute the autocorrelation coefficient to be used by the T-Rex+DA selector
  if (method == "trex+DA+AR1" && is.na(cor_coef)) {
    cor_coef <- abs(mean(apply(X, 1, function(smpl) {
      stats::coef(stats::arima(smpl, order = c(1, 0, 0), include.mean = FALSE, method = "ML"))
    })))
  }
  if (method == "trex+DA+equi" && is.na(cor_coef)) {
    cor_mat <- stats::cor(X)
    cor_coef <- mean(cor_mat[lower.tri(cor_mat, diag = FALSE)])
    rm(cor_mat)
  }

  # Dependency aware relative occurrences for T-Rex+DA+AR1 selector
  if (method == "trex+DA+AR1") {
    kap <- ceiling(log(rho_thr_DA) / log(cor_coef))
    DA_delta_mat <- matrix(NA, nrow = p, ncol = T_stop)
    for (t in seq(T_stop)) {
      for (j in seq(1, p)) {
        sliding_window <- c(
          seq(max(1, j - kap), max(1, j - 1)),
          seq(min(p, j + 1), min(p, j + kap))
        )
        if (j %in% c(1, p)) {
          sliding_window <- sliding_window[-which(sliding_window == j)]
        }
        DA_delta_mat[j, t] <- 2 - min(abs(phi_T_mat[j, t] - phi_T_mat[sliding_window, t]))
      }
    }
    phi_T_mat <- phi_T_mat / DA_delta_mat
    Phi <- Phi / DA_delta_mat[, T_stop]
  }

  # Dependency aware relative occurrences for T-Rex+DA+equi selector
  if (method == "trex+DA+equi") {
    if (abs(cor_coef) > rho_thr_DA) {
      DA_delta_mat <- matrix(NA, nrow = p, ncol = T_stop)
      for (t in seq(T_stop)) {
        for (j in seq(1, p)) {
          sliding_window <- seq(1, p)[-j]
          DA_delta_mat[j, t] <- 2 - min(abs(phi_T_mat[j, t] - phi_T_mat[sliding_window, t]))
        }
      }
      phi_T_mat <- phi_T_mat / DA_delta_mat
      Phi <- Phi / DA_delta_mat[, T_stop]
    }
  }


  if (bootstrap) {
    # Screen-T-Rex with bootstrap
    # Determine the confidence interval grid based on the dummy coefficients using the non-parametric bootstrap
    boot_data <- c(dummy_beta_mat[abs(dummy_beta_mat) > eps ])
    boot_fun <- function(boot_data, idx) {
      mean(boot_data[idx])
    }
    bootstrap_dummy_coef <- boot::boot(data = boot_data,
                                       statistic = boot_fun,
                                       R = R)
    bootstrap_conf_int <- boot::boot.ci(boot.out = bootstrap_dummy_coef,
                                        conf = conf_level_grid,
                                        type = "norm")$normal[ , c(2,3)]

    # Remove all variables whose average coefficient value lies within the bootstrap CI
    # or that have a relative occurrence not larger than 0.5 from the selected active set
    maj_vote_set_zero <- Phi <= 0.5
    beta_means <- colMeans(beta_mat)
    R_no_boot <- sum(!maj_vote_set_zero)
    R_with_boot <- sapply( seq(nrow(bootstrap_conf_int)), function(b) sum(!((beta_means >= bootstrap_conf_int[b, 1]) & (beta_means <= bootstrap_conf_int[b, 2]))) )
    conf_level_index <- which.max(R_with_boot <= R_no_boot)
    conf_level <- conf_level_grid[conf_level_index]
    boot_set_zero <- (beta_means >= bootstrap_conf_int[conf_level_index, 1]) & (beta_means <= bootstrap_conf_int[conf_level_index, 2])
    set_zero <- boot_set_zero # boot_set_zero | maj_vote_set_zero
  } else {
    # Screen-T-Rex without bootstrap
    maj_vote_set_zero <- Phi <= 0.5
    set_zero <- maj_vote_set_zero
  }

  # Selected variables
  selected_var <- rep(0, times = p)
  selected_var[!set_zero] <- 1

  # FDR estimate
  if (sum(!set_zero) == 0) {
    FDR_estimate <- 0
  } else {
    FDR_estimate <- T_stop / sum(!set_zero)
  }

  # List of results
  res <- list(
    selected_var = selected_var,
    FDR_estimate = FDR_estimate,
    dummy_beta_mat = dummy_beta_mat,
    beta_mat = beta_mat,
    Phi = Phi,
    K = K,
    R = R,
    method = method,
    bootstrap = bootstrap,
    conf_level = if (exists("conf_level")) {conf_level} else {NA},
    cor_coef,
    type = type,
    corr_max = corr_max,
    lambda_2_lars = lambda_2_lars,
    rho_thr_DA = rho_thr_DA
  )

  return(res)
}



