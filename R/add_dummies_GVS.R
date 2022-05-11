#' Add GVS-dummies to the original predictor matrix
#'
#' Sample num_dummies dummy vectors for the T-Knock+GVS filter from the univariate standard normal distribution and append them to the predictor matrix X.
#'
#' @param X Real valued predictor matrix.
#' @param num_dummies Number of dummies.
#' @param corr_max Maximum allowed correlation between any two predictors from different clusters.
#'
#' @return Enlarged predictor matrix for the T-Knock+GVS filter, i.e., original predictor matrix and dummies
#'
#' @import stats
#' @import MASS
#'
#' @export
#'
#' @examples
#' n = 50
#' p = 100
#' add_dummies_GVS(X = matrix(rnorm(n * p), nrow = n, ncol = p), L_val = p)
add_dummies_GVS = function(X,
                           num_dummies,
                           corr_max = 0.5) {
  n = nrow(X)
  p = ncol(X)
  if (is.null(names(X))) {
    colnames(X) = paste0('V', seq(p))
  }

  # Single linkage hierarchical clustering using the sample correlation as the similarity measure
  sigma_X = stats::cov(X)
  sigma_X_dist = stats::as.dist(1 - abs(stats::cov2cor(sigma_X)))
  fit = stats::hclust(sigma_X_dist, method = 'single')
  clusters = stats::cutree(fit, h = 1 - corr_max)
  max_clusters = max(clusters)
  clusters = data.frame('Var' = names(clusters),
                        'Cluster_Nr.' = unname(clusters))
  clusters = stats::aggregate(clusters$'Var' ~ clusters$'Cluster_Nr.', FUN = 'c')
  cluster_sizes = vector('numeric', length = max_clusters)
  for (j in seq(max_clusters)) {
    cluster_sizes[j] = length(clusters$`clusters$Var`[[j]])
  }
  w_max = num_dummies / p

  X_p_surrogate = matrix(NA, nrow = n, ncol = p)
  X_Dummy = matrix(NA, nrow = n, ncol = p + num_dummies)
  X_Dummy[, seq(p)] = X

  for (w in seq(w_max)) {
    for (z in seq(max_clusters)) {
      sub_X = X[, clusters$`clusters$Var`[[z]], drop = FALSE]
      sigma_sub_X = stats::cov(sub_X)
      idx = cumsum(cluster_sizes)
      mu = rep(0, times = cluster_sizes[z])
      if (z == 1) {
        X_p_surrogate[, seq(idx[z])] = MASS::mvrnorm(n,
                                                     mu = mu,
                                                     Sigma = sigma_sub_X,
                                                     empirical = FALSE)
      } else{
        X_p_surrogate[, seq(idx[z - 1] + 1, idx[z])] = MASS::mvrnorm(n,
                                                                     mu = mu,
                                                                     Sigma = sigma_sub_X,
                                                                     empirical = FALSE)
      }
    }
    X_Dummy[, seq(w * p + 1, (w + 1) * p)] = X_p_surrogate
  }
  return(X_Dummy)
}
