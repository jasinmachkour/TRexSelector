#' Run the T-Rex selector (\doi{10.48550/arXiv.2110.06048})
#'
#' The T-Rex selector (\doi{10.48550/arXiv.2110.06048}) performs fast variable selection in high-dimensional settings while
#' controlling the false discovery rate (FDR) at a user-defined target level.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param K Number of random experiments.
#' @param max_num_dummies Integer factor determining the maximum number of dummies as a multiple of the number of original variables p
#' (i.e., num_dummies = max_num_dummies * p).
#' @param max_T_stop If TRUE the maximum number of dummies that can be included before stopping is set to ceiling(n / 2),
#' where n is the number of data points/observations.
#' @param method 'trex' for the T-Rex selector (\doi{10.48550/arXiv.2110.06048}),
#' 'trex+GVS' for the T-Rex+GVS selector (\doi{10.23919/EUSIPCO55093.2022.9909883}),
#' 'trex+DA+AR1' for the T-Rex+DA+AR1 selector,
#' 'trex+DA+equi' for the T-Rex+DA+equi selector,
#' 'trex+DA+BT' for the T-Rex+DA+BT selector (\doi{10.48550/arXiv.2401.15796}),
#' 'trex+DA+NN' for the T-Rex+DA+NN selector (\doi{10.48550/arXiv.2401.15139}).
#' @param GVS_type 'IEN' for the Informed Elastic Net (\doi{10.1109/CAMSAP58249.2023.10403489}),
#' 'EN' for the ordinary Elastic Net (\doi{10.1111/j.1467-9868.2005.00503.x}).
#' @param cor_coef AR(1) autocorrelation coefficient for the T-Rex+DA+AR1 selector or equicorrelation coefficient for the T-Rex+DA+equi selector.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters (for method = 'trex+GVS').
#' @param lambda_2_lars lambda_2-value for LARS-based Elastic Net.
#' @param rho_thr_DA Correlation threshold for the T-Rex+DA+AR1 selector and the T-Rex+DA+equi selector (i.e., method = 'trex+DA+AR1' or 'trex+DA+equi').
#' @param hc_dist Distance measure of the hierarchical clustering/dendrogram (only for trex+DA+BT):
#' 'single' for single-linkage, "complete" for complete linkage, "average" for average linkage (see [hclust] for more options).
#' @param hc_grid_length Length of the height-cutoff-grid for the dendrogram (integer between 1 and the number of original variables p).
#' @param parallel_process Logical. If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing.
#' @param seed Seed for random number generator (ignored if parallel_process = FALSE).
#' @param eps Numerical zero.
#' @param verbose Logical. If TRUE progress in computations is shown.
#'
#' @return A list containing the estimated support vector and additional information, including the number of used dummies and the number of included dummies before stopping.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParWorkers registerDoSEQ
#' @importFrom stats coef arima cor as.dist hclust
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- c(Gauss_data$y)
#' set.seed(1234)
#' res <- trex(X = X, y = y)
#' selected_var <- res$selected_var
#' selected_var
trex <- function(X,
                 y,
                 tFDR = 0.2,
                 K = 20,
                 max_num_dummies = 10,
                 max_T_stop = TRUE,
                 method = "trex",
                 GVS_type = "IEN",
                 cor_coef = NA,
                 type = "lar",
                 corr_max = 0.5,
                 lambda_2_lars = NULL,
                 rho_thr_DA = 0.02,
                 hc_dist = "single",
                 hc_grid_length = min(20, ncol(X)),
                 parallel_process = FALSE,
                 parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                 seed = NULL,
                 eps = .Machine$double.eps,
                 verbose = TRUE) {
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

  if (length(tFDR) != 1 ||
      tFDR < 0 ||
      tFDR > 1) {
    stop("'tFDR' must be a number between 0 and 1 (including 0 and 1).")
  }

  if (length(K) != 1 ||
      K < 2 ||
      K %% 1 != 0) {
    stop("The number of random experiments 'K' must be an integer larger or equal to 2.")
  }

  if (length(max_num_dummies) != 1 ||
      max_num_dummies < 1 ||
      max_num_dummies %% 1 != 0) {
    stop("'max_num_dummies' must be an integer larger or equal to 1.")
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

  if (parallel_process &&
      (length(parallel_max_cores) != 1 ||
       parallel_max_cores %% 1 != 0 ||
       parallel_max_cores < 2)) {
    stop(
      "For parallel processing at least two workers have to be registered:
         'parallel_max_cores' must be an integer larger or equal to 2."
    )
  }

  if (parallel_process &&
      parallel_max_cores > min(K, max(
        1, parallel::detectCores(logical = FALSE)
      ))) {
    parallel_max_cores_modified <-
      min(K, max(1, parallel::detectCores(logical = FALSE)))
    message(
      paste0(
        "For computing ",
        K,
        " random experiments, it is not useful/possible to register ",
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

  # Scale X and center y
  X <- scale(X)
  y <- y - mean(y)

  # Number of rows n and columns p of X
  n <- nrow(X)
  p <- ncol(X)

  # T-Rex+DA (binary tree)
  if (method == "trex+DA+BT") {
    # Dendrogram for trex+DA+BT
    if (system.file(package = "WGCNA") == "") {
      message("To speed up computations, please install the package 'WGCNA'.")
      cor_mat <- stats::cor(X)
    } else {
      cor_mat <- WGCNA::cor(X)
    }
    cor_mat_distance <- stats::as.dist(1 - abs(cor_mat))
    if (system.file(package = "fastcluster") == "") {
      message("To speed up computations, please install the package 'fastcluster'.")
      dendrogram <- stats::hclust(cor_mat_distance, method = hc_dist)
    } else {
      dendrogram <- fastcluster::hclust(cor_mat_distance, method = hc_dist)
    }
    rho_grid_subsample <- round(seq(1, p, length.out = hc_grid_length))
    rho_grid_len <- hc_grid_length
    rho_grid <- c(1 - rev(dendrogram$height), 1)[rho_grid_subsample]
    clusters <- stats::cutree(dendrogram, h = 1 - rho_grid)

    gr_j_list <-
      lapply(seq(1, p), FUN = function (j) {
        lapply(seq(1, rho_grid_len), FUN = function(x) {
          gr_num_j <- clusters[j, x]
          gr_j <- which(clusters[ , x] == gr_num_j)
          gr_j <- gr_j[-which(gr_j == j)]
        })
      })

    # Closest correlation point to reference point (only for trex+DA+BT) for determining number of dummies
    opt_point_BT <- round(0.75 * rho_grid_len)
  }

  # T-Rex+DA (nearest neighbors)
  if (method == "trex+DA+NN") {
    # Nearest neighbors (NN) groups for trex+DA+NN
    if (system.file(package = "WGCNA") == "") {
      message("To speed up computations, please install the package 'WGCNA'.")
      cor_mat <- stats::cor(X)
    } else {
      cor_mat <- WGCNA::cor(X)
    }
    rho_grid_len <- hc_grid_length
    rho_grid <- seq(0, 1, length.out = rho_grid_len)
    gr_j_list <-
      lapply(seq(1, p), FUN = function (j) {
        lapply(seq(1, rho_grid_len), FUN = function(x) {
          gr_j <- which(abs(cor_mat[ , j]) >= rho_grid[x])
          gr_j <- gr_j[-which(gr_j == j)]
        })
      })
    opt_point_BT <- round(0.75 * rho_grid_len)
  }

  # Voting level grid
  V <- seq(0.5, 1 - eps, by = 1 / K)
  V_len <- length(V)

  # Initialize L-loop
  LL <- 1
  T_stop <- 1

  if (method == "trex+DA+AR1" && is.na(cor_coef)) {
    cor_coef <- abs(mean(apply(X, 1, function(smpl) {
      stats::coef(stats::arima(smpl, order = c(1, 0, 0), include.mean = FALSE, method = "ML"))
    })))
  }
  if (method == "trex+DA+equi" && is.na(cor_coef)) {
    if (system.file(package = "WGCNA") == "") {
      message("To speed up computations, please install the package 'WGCNA'.")
      cor_mat <- stats::cor(X)
    } else {
      cor_mat <- WGCNA::cor(X)
    }
    cor_coef <- mean(cor_mat[lower.tri(cor_mat, diag = FALSE)])
    rm(cor_mat)
  }
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    FDP_hat <- matrix(NA, nrow = V_len, ncol = rho_grid_len)
  } else {
    FDP_hat <- rep(NA, times = V_len)
  }

  # 75% voting reference point for determining number of dummies
  opt_point <- which(abs(V - 0.75) < eps)
  if (length(opt_point) == 0) {
    # If 75% optimization point does not exist, choose closest optimization point lower than 75%
    opt_point <- length(V[V < 0.75])
  }

  # Setup cluster
  if (parallel_process && foreach::getDoParWorkers() == 1) {
    cl <- parallel::makeCluster(parallel_max_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  # FDP larger than tFDR
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    fdp_larger_tFDR <- FDP_hat[opt_point, opt_point_BT] > tFDR
  } else {
    fdp_larger_tFDR <- FDP_hat[opt_point] > tFDR
  }

  while ((LL <= max_num_dummies && fdp_larger_tFDR) ||
         sum(!is.na(FDP_hat)) == 0) {

    num_dummies <- LL * p
    LL <- LL + 1

    # K Random experiments
    suppressWarnings(
      rand_exp <- random_experiments(
        X = X,
        y = y,
        K = K,
        T_stop = T_stop,
        num_dummies = num_dummies,
        method = method,
        GVS_type = GVS_type,
        type = type,
        corr_max = corr_max,
        lambda_2_lars = lambda_2_lars,
        early_stop = TRUE,
        verbose = verbose,
        intercept = FALSE,
        standardize = TRUE,
        parallel_process = parallel_process,
        parallel_max_cores = parallel_max_cores,
        seed = seed,
        eps = eps
      )
    )
    phi_T_mat <- rand_exp$phi_T_mat
    Phi <- rand_exp$Phi

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

    # Dependency aware relative occurrences for T-Rex+DA+BT or T-Rex+DA+NN selector
    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      DA_delta_mat_BT <- matrix(NA, nrow = p, ncol = rho_grid_len)
      for (j in seq(1, p)) {
        DA_delta_mat_BT[j, ] <-
          sapply(gr_j_list[[j]], FUN = function(x) {
            if (length(x) == 0) {
              2
            } else {
              2 - min(abs(phi_T_mat[j, T_stop] - phi_T_mat[x, T_stop]))
            }
          })
      }
      phi_T_array_BT <- array(apply(DA_delta_mat_BT, 2, function(x) phi_T_mat / x), dim = c(p, T_stop, rho_grid_len))
      Phi_BT <- Phi / DA_delta_mat_BT
    }

    # Phi_prime and FDP_hat
    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      # Phi_prime
      Phi_prime <- matrix(sapply(seq(rho_grid_len), FUN = function(x) {
        Phi_prime_fun(
          p = p,
          T_stop = T_stop,
          num_dummies = num_dummies,
          phi_T_mat = matrix(phi_T_array_BT[ , , x], nrow = p, ncol = T_stop),
          Phi = Phi_BT[ , x],
          eps = eps
        )
      }), nrow = p, ncol = rho_grid_len)

      # FDP_hat
      FDP_hat <- matrix(sapply(seq(rho_grid_len), FUN = function(x) {
        fdp_hat(
          V = V,
          Phi = Phi_BT[ , x],
          Phi_prime = Phi_prime[ , x]
        )
      }), nrow = V_len, ncol = rho_grid_len)
    } else {
      # Phi_prime
      Phi_prime <- Phi_prime_fun(
        p = p,
        T_stop = T_stop,
        num_dummies = num_dummies,
        phi_T_mat = phi_T_mat,
        Phi = Phi,
        eps = eps
      )

      # FDP_hat
      FDP_hat <- fdp_hat(
        V = V,
        Phi = Phi,
        Phi_prime = Phi_prime
      )
    }

    # FDP larger than tFDR
    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      fdp_larger_tFDR <- FDP_hat[opt_point, opt_point_BT] > tFDR
    } else {
      fdp_larger_tFDR <- FDP_hat[opt_point] > tFDR
    }

    # Print number of appended dummies by the extended calibration algorithm of the T-Rex selector
    if (verbose) {
      cat(paste("\n Appended dummies:", num_dummies, "\n"))
    }
  }

  # Initialize T-loop
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    FDP_hat_array_BT <- array(FDP_hat, dim = c(dim(FDP_hat), 1))
    Phi_array_BT <- array(Phi_BT, dim = c(dim(Phi_BT), 1))
  } else {
    FDP_hat_mat <- matrix(FDP_hat, nrow = 1)
    Phi_mat <- matrix(Phi, nrow = 1)
  }

  if (max_T_stop) {
    max_T <- min(num_dummies, ceiling(n / 2))
  } else {
    max_T <- num_dummies
  }

  # Reset seed
  if (!is.null(seed)) {
    seed <- seed + 12345
  }

  # FDP lower than target FDR?
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    fdp_lower_tFDR <- FDP_hat[V_len, opt_point_BT] <= tFDR
  } else {
    fdp_lower_tFDR <- FDP_hat[V_len] <= tFDR
  }

  while (fdp_lower_tFDR && (T_stop < max_T)) {
    T_stop <- T_stop + 1

    # K Random experiments
    suppressWarnings(
      rand_exp <- random_experiments(
        X = X,
        y = y,
        K = K,
        T_stop = T_stop,
        num_dummies = num_dummies,
        method = method,
        GVS_type = GVS_type,
        type = type,
        corr_max = corr_max,
        lambda_2_lars = lambda_2_lars,
        early_stop = TRUE,
        lars_state_list = rand_exp$lars_state_list,
        verbose = verbose,
        intercept = FALSE,
        standardize = TRUE,
        parallel_process = parallel_process,
        parallel_max_cores = parallel_max_cores,
        seed = seed,
        eps = eps
      )
    )
    phi_T_mat <- rand_exp$phi_T_mat
    Phi <- rand_exp$Phi

    # Dependency aware relative occurrences for the T-Rex+DA selector
    if (method == "trex+DA+AR1") {
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

    # Dependency aware relative occurrences for T-Rex+DA+BT or T-Rex+DA+NN selector
    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      DA_delta_mat_BT <- matrix(NA, nrow = p, ncol = rho_grid_len)
      for (j in seq(1, p)) {
        DA_delta_mat_BT[j, ] <-
          sapply(gr_j_list[[j]], FUN = function(x) {
            if (length(x) == 0) {
              2
            } else {
              2 - min(abs(phi_T_mat[j, T_stop] - phi_T_mat[x, T_stop]))
            }
          })
      }
      phi_T_array_BT <- array(apply(DA_delta_mat_BT, 2, function(x) phi_T_mat / x), dim = c(p, T_stop, rho_grid_len))
      Phi_BT <- Phi / DA_delta_mat_BT
    }

    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      Phi_array_BT <- array(c(Phi_array_BT, Phi_BT), dim = dim(Phi_array_BT) + c(0, 0, 1))

      # Phi_prime
      Phi_prime <- matrix(sapply(seq(rho_grid_len), FUN = function(x) {
        Phi_prime_fun(
          p = p,
          T_stop = T_stop,
          num_dummies = num_dummies,
          phi_T_mat = matrix(phi_T_array_BT[ , , x], nrow = p, ncol = T_stop),
          Phi = Phi_BT[ , x],
          eps = eps
        )
      }), nrow = p, ncol = rho_grid_len)

      # FDP_hat
      FDP_hat <- matrix(sapply(seq(rho_grid_len), FUN = function(x) {
        fdp_hat(
          V = V,
          Phi = Phi_BT[ , x],
          Phi_prime = Phi_prime[ , x]
        )
      }), nrow = V_len, ncol = rho_grid_len)

      FDP_hat_array_BT <- array(c(FDP_hat_array_BT, FDP_hat), dim = dim(FDP_hat_array_BT) + c(0, 0, 1))
    } else {
      Phi_mat <- rbind(Phi_mat, Phi)

      # Phi_prime
      Phi_prime <- Phi_prime_fun(
        p = p,
        T_stop = T_stop,
        num_dummies = num_dummies,
        phi_T_mat = phi_T_mat,
        Phi = Phi,
        eps = eps
      )

      # FDP_hat
      FDP_hat <- fdp_hat(
        V = V,
        Phi = Phi,
        Phi_prime = Phi_prime
      )

      FDP_hat_mat <- rbind(FDP_hat_mat, FDP_hat)
    }

    # FDP lower than target FDR?
    if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
      fdp_lower_tFDR <- FDP_hat[V_len, opt_point_BT] <= tFDR
    } else {
      fdp_lower_tFDR <- FDP_hat[V_len] <= tFDR
    }

    # Print the number of by the extended calibration algorithm of the T-Rex selector included dummies along the selection process before stopping
    if (verbose) {
      cat(paste("\n Included dummies before stopping:", T_stop, "\n"))
    }
  }

  # "Transpose" Phi_array_BT and FDP_hat_array_BT (only for "trex+DA+BT" and "trex+DA+NN")
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    Phi_array_BT <- aperm(Phi_array_BT, perm = c(3, 1, 2))
    FDP_hat_array_BT <- aperm(FDP_hat_array_BT, perm = c(3, 1, 2))
  }

  # T-Rex: Select variables
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    res_T_dummy <- select_var_fun_DA_BT(
      p = p,
      tFDR = tFDR,
      T_stop = T_stop,
      FDP_hat_array_BT = FDP_hat_array_BT,
      Phi_array_BT = Phi_array_BT,
      V = V,
      rho_grid = rho_grid
    )
    selected_var <- res_T_dummy$selected_var
    v_thresh <- res_T_dummy$v_thresh
    rho_thresh <- res_T_dummy$rho_thresh
    R_array <- res_T_dummy$R_array
  } else {
    res_T_dummy <- select_var_fun(
      p = p,
      tFDR = tFDR,
      T_stop = T_stop,
      FDP_hat_mat = FDP_hat_mat,
      Phi_mat = Phi_mat,
      V = V
    )
    selected_var <- res_T_dummy$selected_var
    v_thresh <- res_T_dummy$v_thresh
    rho_thresh <- NA
    R_mat <- res_T_dummy$R_mat
  }

  # List of results
  if (method %in% c("trex+DA+BT", "trex+DA+NN")) {
    res <- list(
      selected_var = selected_var,
      tFDR = tFDR,
      T_stop = T_stop,
      num_dummies = num_dummies,
      V = V,
      rho_grid = rho_grid,
      v_thresh = v_thresh,
      rho_thresh = rho_thresh,
      #
      FDP_hat_array_BT = FDP_hat_array_BT,
      Phi_array_BT = Phi_array_BT,
      R_array = R_array,
      phi_T_array_BT = phi_T_array_BT,
      #
      Phi_prime = Phi_prime,
      method = method,
      GVS_type = GVS_type,
      cor_coef = cor_coef,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      rho_thr_DA = rho_thr_DA,
      hc_dist = hc_dist
    )
  } else {
    res <- list(
      selected_var = selected_var,
      tFDR = tFDR,
      T_stop = T_stop,
      num_dummies = num_dummies,
      V = V,
      v_thresh = v_thresh,
      rho_thresh = rho_thresh,
      #
      FDP_hat_mat = FDP_hat_mat,
      Phi_mat = Phi_mat,
      R_mat = R_mat,
      phi_T_mat = phi_T_mat,
      #
      Phi_prime = Phi_prime,
      method = method,
      GVS_type = GVS_type,
      cor_coef = cor_coef,
      type = type,
      corr_max = corr_max,
      lambda_2_lars = lambda_2_lars,
      rho_thr_DA = rho_thr_DA,
      hc_dist = hc_dist
    )
  }

  return(res)
}

