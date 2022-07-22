#' Add dummy predictors to the original predictor matrix, as required by the T-Rex+GVS selector
#'
#' Generate num_dummies dummy predictors as required for the T-Rex+GVS selector and append them to the predictor matrix X.
#'
#' @param X Real valued predictor matrix.
#' @param num_dummies Number of dummies that are appended to the predictor matrix. Has to be a multiple of the number of original variables.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#'
#' @return Enlarged predictor matrix for the T-Rex+GVS selector.
#'
#' @importFrom stats cov as.dist hclust cutree aggregate
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 100
#' X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' add_dummies_GVS(X = X, num_dummies = p)
add_dummies_GVS <- function(X,
                            num_dummies,
                            corr_max = 0.5) {
  # Error control
  if (!is.matrix(X)) {
    stop("'X' must be a matrix.")
  }

  if (!is.numeric(X)) {
    stop("'X' only allows numerical values.")
  }

  if (anyNA(X)) {
    stop("'X' contains NAs. Please remove or impute them before proceeding.")
  }

  # Dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # Continue error control
  if (length(num_dummies) != 1 ||
    num_dummies %% p != 0 ||
    num_dummies < 1) {
    stop(
      "'num_dummies' must be a positive integer multiple of the total number of original predictors in X."
    )
  }

  if (length(corr_max) != 1 ||
    corr_max < 0 ||
    corr_max > 1) {
    stop("'corr_max' must have a value between zero and one.")
  }

  # Add variable names
  if (is.null(names(X))) {
    colnames(X) <- paste0("V", seq(p))
  }

  # Single linkage hierarchical clustering using the sample correlation as the similarity measure
  sigma_X <- stats::cov(X)
  sigma_X_dist <- stats::as.dist(1 - abs(stats::cov2cor(sigma_X)))
  fit <- stats::hclust(sigma_X_dist, method = "single")
  clusters <- stats::cutree(fit, h = 1 - corr_max)
  max_clusters <- max(clusters)
  clusters <- data.frame(
    "Var" = names(clusters),
    "Cluster_Nr." = unname(clusters)
  )
  clusters <-
    stats::aggregate(clusters$"Var" ~ clusters$"Cluster_Nr.",
      FUN = "c",
      simplify = FALSE
    )
  cluster_sizes <- vector("numeric", length = max_clusters)
  for (j in seq(max_clusters)) {
    cluster_sizes[j] <- length(clusters$`clusters$Var`[[j]])
  }

  # Generate dummy predictors and append them to the original predictor matrix X
  w_max <- num_dummies / p
  X_p_sub_dummy <- matrix(NA, nrow = n, ncol = p)
  X_Dummy <- matrix(NA, nrow = n, ncol = p + num_dummies)
  X_Dummy[, seq(p)] <- X

  for (w in seq(w_max)) {
    for (z in seq(max_clusters)) {
      sub_X <- X[, clusters$`clusters$Var`[[z]], drop = FALSE]
      sigma_sub_X <- stats::cov(sub_X)
      idx <- cumsum(cluster_sizes)
      mu <- rep(0, times = cluster_sizes[z])
      if (z == 1) {
        X_p_sub_dummy[, seq(idx[z])] <- MASS::mvrnorm(n,
          mu = mu,
          Sigma = sigma_sub_X,
          empirical = FALSE
        )
      } else {
        X_p_sub_dummy[, seq(idx[z - 1] + 1, idx[z])] <- MASS::mvrnorm(n,
          mu = mu,
          Sigma = sigma_sub_X,
          empirical = FALSE
        )
      }
    }
    X_Dummy[, seq(w * p + 1, (w + 1) * p)] <- X_p_sub_dummy
  }

  return(X_Dummy)
}
