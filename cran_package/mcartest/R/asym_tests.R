#' Asymptotic L2 test for MCAR
#'
#' Tests \eqn{H_0:\ \mu_A=\mu_B} via \eqn{T_{\mu,L2}=n\\lVert \hat\mu_A-\hat\mu_B\\rVert^2_{L2}}.
#' p-values are obtained from a KL-mixture; the number of components is chosen by FVE.
#'
#' @inheritParams mcar_common-params
#' @param fve Fraction of variance explained (0-1) to choose \eqn{q}.
#' @param n_sim Number of Monte Carlo draws for the KL-mixture.
#' @param seed RNG seed.
#' @return `htest` (with extras).
#' @examples
#' set.seed(1)
#' m <- 50
#' n <- 40
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) {
#'   d <- diff(g)[1]
#'   c(0, cumsum(rnorm(length(g) - 1, sd = sqrt(d))))
#' }
#' X <- t(replicate(n, bm(grid)))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X[O == 0L] <- NA_real_
#'
#' # Asymptotic L2 test
#' res_L2 <- asym_mean_L2_test(
#'   X = X,
#'   n_sim = 2000,
#'   seed = 123
#' )
#' res_L2$p.value
#' @export
asym_mean_L2_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, n_sim = 10000,
                              min_frac = 0.10, seed = NULL) {
  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)

  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid <- grid[idx]
  X_sub <- X[, idx, drop = FALSE]
  O_sub <- O[, idx, drop = FALSE]

  est <- .group_mean_estimators(X_sub, O_sub, group_A)
  mean_A <- est$mean_A
  mean_B <- est$mean_B
  mean_diff <- mean_A - mean_B

  mean_A_tfd <- tf::tfd(matrix(mean_A, nrow = 1), arg = subgrid)
  mean_B_tfd <- tf::tfd(matrix(mean_B, nrow = 1), arg = subgrid)
  mean_diff_tfd <- tf::tfd(matrix(mean_diff, nrow = 1), arg = subgrid)

  T_L2 <- n * tf::tf_integrate(mean_diff_tfd^2, arg = subgrid)

  cov_matrix <- .covariance_estimator(X_sub, O_sub, group_A, mean_A, mean_B, est$pA, est$pB)
  KL_decomp <- .kl_decomposition(cov_matrix, subgrid)
  eigenvalues <- KL_decomp$eigenvalues

  if (!is.null(seed)) set.seed(seed)
  cum_var <- cumsum(eigenvalues) / sum(eigenvalues)
  q <- which(cum_var >= fve)[1]
  eigenvalues_fve <- eigenvalues[seq_len(q)]

  Z <- matrix(rnorm(q * n_sim), nrow = q)
  sim_stats <- colSums((Z^2) * eigenvalues_fve)

  p_value <- (sum(sim_stats >= T_L2) + 1) / (length(sim_stats) + 1)

  data_name <- if (!is.null(fd)) "fd" else "X"

  output <- .create_output(stat_name = paste0("T_\u03bc,L\u00B2"), 
                           stat_value = T_L2, 
                           p_value = p_value,
                           method = "L2 test",
                           data_name = data_name,
                           estimate = list(
                             mean_A = mean_A_tfd,
                             mean_B = mean_B_tfd,
                             mean_diff = mean_diff_tfd
                             ),
                           parameter = c(q = q, n_sim = n_sim)
                           )
  output
}


#' Asymptotic supremum test for MCAR (optional simultaneous bands)
#'
#' Tests \eqn{H_0:\ \mu_A=\mu_B} via \eqn{T_{\mu,D}=\sqrt{n}\,\lVert \hat\mu_A-\hat\mu_B\rVert_\infty}.
#' KL-based GP approximation; with `compute_bands=TRUE` also returns simultaneous bands.
#'
#' @inheritParams mcar_common-params
#' @param fve Fraction of variance explained (0-1) to select KL components.
#' @param n_sim Number of Monte Carlo draws for the GP approx.
#' @param seed RNG seed.
#' @param alpha Significance level for bands.
#' @param compute_bands Logical: also compute simultaneous bands?
#' @param bands_only Logical: return only band info?
#' @return `htest` (with extras) or a light band list.
#' @examples
#' set.seed(1)
#' m <- 50
#' n <- 40
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) {
#'   d <- diff(g)[1]
#'   c(0, cumsum(rnorm(length(g) - 1, sd = sqrt(d))))
#' }
#' X <- t(replicate(n, bm(grid)))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X[O == 0L] <- NA_real_
#'
#' # Asymptotic Supremum test with simultaneous bands
#' res_sup <- asym_mean_sup_test(
#'   X = X,
#'   n_sim = 2000,
#'   compute_bands = TRUE,
#'   seed = 123
#' )
#' res_sup$p.value
#'
#' # Quick plot: mean difference + 95% simultaneous band
#' diff_hat <- tf::tf_evaluate(res_sup$estimate$mean_diff, arg = res_sup$bands$grid)[[1]]
#' plot(res_sup$bands$grid, diff_hat,
#'   type = "l",
#'   ylim = range(c(res_sup$bands$lower, res_sup$bands$upper)),
#'   xlab = "t", ylab = "mean difference"
#' )
#' abline(h = 0, lty = 3)
#' lines(res_sup$bands$grid, res_sup$bands$lower, lty = 2)
#' lines(res_sup$bands$grid, res_sup$bands$upper, lty = 2)
#' @export
asym_mean_sup_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                               fve = 0.99, n_sim = 10000,
                               min_frac = 0.10, seed = NULL, alpha = 0.05,
                               compute_bands = TRUE, bands_only = FALSE) {
  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)

  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid <- grid[idx]
  X_sub <- X[, idx, drop = FALSE]
  O_sub <- O[, idx, drop = FALSE]

  est <- .group_mean_estimators(X_sub, O_sub, group_A)
  mean_A <- est$mean_A
  mean_B <- est$mean_B
  mean_diff <- mean_A - mean_B

  mean_A_tfd <- tf::tfd(matrix(mean_A, nrow = 1), arg = subgrid)
  mean_B_tfd <- tf::tfd(matrix(mean_B, nrow = 1), arg = subgrid)
  mean_diff_tfd <- tf::tfd(matrix(mean_diff, nrow = 1), arg = subgrid)

  T_D <- sqrt(n) * max(abs(mean_diff))

  cov_matrix <- .covariance_estimator(X_sub, O_sub, group_A, mean_A, mean_B, est$pA, est$pB)
  KL_decomp <- .kl_decomposition(cov_matrix, subgrid)
  eigenvalues <- KL_decomp$eigenvalues
  eigenfunctions <- KL_decomp$eigenfunctions

  if (!is.null(seed)) set.seed(seed)
  cum_var <- cumsum(eigenvalues) / sum(eigenvalues)
  q <- which(cum_var >= fve)[1]
  eigenvalues_fve <- eigenvalues[seq_len(q)]
  eigenfunctions_fve <- eigenfunctions[, seq_len(q), drop = FALSE]
  A <- sweep(eigenfunctions_fve, 2, sqrt(eigenvalues_fve), "*")

  Z <- matrix(rnorm(q * n_sim), nrow = q)
  gp_vals <- A %*% Z
  sim_stats <- apply(abs(gp_vals), 2, max)

  p_value <- (sum(sim_stats >= T_D) + 1) / (length(sim_stats) + 1)

  if (isTRUE(compute_bands)) {
    bands <- .confidence_bands("D", mean_diff, sim_stats, n, alpha, subgrid)
  } else {
    bands <- list(lower = NULL, upper = NULL, band = NULL)
  }

  data_name <- if (!is.null(fd)) "fd" else "X"

  if (isTRUE(bands_only)) {
    return(list(
      estimate = list(
        mean_A = mean_A_tfd,
        mean_B = mean_B_tfd,
        mean_diff = mean_diff_tfd
      ),
      parameter = c(q = q, n_sim = n_sim),
      bands = bands
    ))
  }

  output <- .create_output(stat_name = paste0("T_\u03bc,D"), 
                           stat_value = T_D,
                           p_value = p_value,
                           method = "Supremum test",
                           data_name = data_name,
                           estimate = list(
                             mean_A = mean_A_tfd,
                             mean_B = mean_B_tfd,
                             mean_diff = mean_diff_tfd
                             ),
                           parameter = c(q = q, n_sim = n_sim),
                           bands = bands
                           )
  output
}
