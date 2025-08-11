#' Tidyfun extensions: Algorithms 1, 2, and 3 (auto inputs)
#'
#' Implements:
#'  - Algorithm 1: L2 test of equality of mean functions with available-data estimators
#'  - Algorithm 2: Supremum test (T_{mu,D}) with KL-based Gaussian process approximation
#'  - Algorithm 3: Simultaneous confidence bands for mu_A - mu_B
#'
#' Changes vs. original: O_mat, group_A, and grid are created **internally** from X_obs.
#'  - O_mat := 1 if X_obs is observed (non-NA) and 0 otherwise.
#'  - group_A (Example 1 in Ofner et al., 2025): row is in A iff it is **complete** (no NA),
#'    otherwise it is in B (incomplete).
#'  - grid := seq(0, 1, length.out = ncol(X_obs)).
#'
#' The exported functions allow optional overrides for O_mat, group_A, and grid, but when NULL
#' the above automatic construction is used. This keeps backward-compatibility and supports
#' one-argument calls with just X_obs.
#'
#' The functions operate on matrices of observed values (with NA for missing) and the common grid.
#' Integration uses tf::tf_integrate for numerical stability and API consistency.
#'
#' @name tidyfun_ext_algo1to3
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------- internal utilities (no export) -------------------------------------

#' @keywords internal
#' @noRd
.tfu_assert_inputs <- function(X_obs, O_mat, group_A, grid) {
  stopifnot(is.matrix(X_obs) || is.data.frame(X_obs))
  if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
  stopifnot(is.matrix(O_mat), is.numeric(grid))
  stopifnot(nrow(X_obs) == nrow(O_mat), ncol(X_obs) == ncol(O_mat))
  stopifnot(length(group_A) == nrow(X_obs))
  stopifnot(all(is.finite(grid)), length(unique(grid)) == length(grid))
}

#' @keywords internal
#' @noRd
.tfu_group_from_delta <- function(O_mat, delta) {
  stopifnot(is.numeric(delta), length(delta) == 1L, delta > 0, delta <= 1)
  rowMeans(O_mat != 0) >= delta
}

#' Build/validate O_mat, group_A, grid from X_obs when missing
#' @keywords internal
#' @noRd
.tfu_prepare_inputs <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL, delta_A = 1) {
  # ensure numeric matrix
  if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
  storage.mode(X_obs) <- "numeric"
  n <- nrow(X_obs); m <- ncol(X_obs)
  
  # O_mat: 1 if observed, 0 if NA
  if (is.null(O_mat)) {
    O_mat <- 1L * (!is.na(X_obs))
  } else {
    stopifnot(is.matrix(O_mat), nrow(O_mat) == n, ncol(O_mat) == m)
    O_mat <- 1L * (O_mat != 0)
  }
  
  # group_A: Example 1 (delta_A = 1) oder Example 2 (delta_A < 1)
  if (is.null(group_A)) {
    group_A <- .tfu_group_from_delta(O_mat, delta_A)
  } else {
    stopifnot(length(group_A) == n)
    group_A <- as.logical(group_A)
  }
  
  # spezifische Fehlermeldungen, falls eine Gruppe leer ist
  nA <- sum(group_A)
  nB <- sum(!group_A)
  if (nA == 0L || nB == 0L) {
    obs_frac <- rowMeans(O_mat != 0)
    rng <- range(obs_frac)
    if (nA == 0L) {
      stop(sprintf(
        "Automatische Gruppierung ergab eine leere Gruppe A (keine Zeilen mit Beobachtungsanteil >= delta_A).\nTipp: 'delta_A' verringern oder 'group_A' explizit übergeben.\nDiag: delta_A = %.3f, range(obs_frac) = [%.3f, %.3f], nA = %d, nB = %d.",
        delta_A, rng[1], rng[2], nA, nB
      ))
    } else {
      stop(sprintf(
        "Automatische Gruppierung ergab eine leere Gruppe B (alle Zeilen haben Beobachtungsanteil >= delta_A).\nTipp: 'delta_A' erhöhen oder 'group_A' explizit übergeben.\nDiag: delta_A = %.3f, range(obs_frac) = [%.3f, %.3f], nA = %d, nB = %d.",
        delta_A, rng[1], rng[2], nA, nB
      ))
    }
  }
  
  # grid: equidistant on [0,1]
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = m)
  } else {
    stopifnot(is.numeric(grid), length(grid) == m)
  }
  
  .tfu_assert_inputs(X_obs, O_mat, group_A, grid)
  list(X_obs = X_obs, O_mat = O_mat, group_A = group_A, grid = grid)
}

#' Weighted (trapezoid) quadrature weights on possibly non-uniform grid
#' @keywords internal
#' @noRd
.tfu_trap_weights <- function(x) {
  m <- length(x)
  if (m == 1L) return(1)
  w <- numeric(m)
  w[1] <- (x[2] - x[1]) / 2
  w[m] <- (x[m] - x[m-1]) / 2
  if (m > 2L) w[2:(m-1)] <- (x[3:m] - x[1:(m-2)]) / 2
  w
}

#' Available-data means and observation probabilities p-hat_A, p-hat_B
#' (denominator is total n, not n_A / n_B)
#' @keywords internal
#' @noRd
.tfu_available_means <- function(X_obs, O_mat, group_A, eps = 1e-8) {
  n  <- nrow(X_obs)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  pA_hat <- colSums(O_mat * IA) / n
  pB_hat <- colSums(O_mat * IB) / n
  pA_hat <- pmax(pA_hat, eps)
  pB_hat <- pmax(pB_hat, eps)
  
  muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_hat)
  
  list(muA = muA_hat, muB = muB_hat, pA = pA_hat, pB = pB_hat)
}

#' Corrected covariance estimator K(s,t) per paper, vectorized
#' @keywords internal
#' @noRd
.tfu_corrected_cov <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X_obs)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  Xc <- X_obs
  A_idx <- which(IA == 1)
  B_idx <- which(IB == 1)
  if (length(A_idx)) Xc[A_idx, ] <- sweep(X_obs[A_idx, , drop = FALSE], 2, muA_hat, `-`)
  if (length(B_idx)) Xc[B_idx, ] <- sweep(X_obs[B_idx, , drop = FALSE], 2, muB_hat, `-`)
  
  # zero-out unobserved & NA
  Xc[O_mat == 0] <- 0
  Xc[is.na(Xc)]  <- 0
  
  # group splits
  XA <- Xc[A_idx, , drop = FALSE]
  XB <- Xc[B_idx, , drop = FALSE]
  
  # column-wise scaling by p-hat
  if (nrow(XA)) XA <- sweep(XA, 2, pA_hat, "/")
  if (nrow(XB)) XB <- sweep(XB, 2, pB_hat, "/")
  
  Sum_A <- if (nrow(XA)) crossprod(XA) else matrix(0, ncol(X_obs), ncol(X_obs))
  Sum_B <- if (nrow(XB)) crossprod(XB) else matrix(0, ncol(X_obs), ncol(X_obs))
  
  K <- (Sum_A + Sum_B) / n
  # symmetrize (numeric hygiene)
  K <- (K + t(K)) / 2
  K
}

#' KL basis from covariance with trapezoid weights (W^{1/2} K W^{1/2})
#' returns eigenvalues lam (>=0) and eigenfunctions phi (L2_w-orthonormal)
#' @keywords internal
#' @noRd
.tfu_kl_from_cov <- function(K, grid) {
  w <- .tfu_trap_weights(grid)
  sw <- sqrt(w)
  
  # S = W^{1/2} K W^{1/2}
  S <- (sw * t(sw * K))
  S <- (S + t(S)) / 2
  
  ev <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  U   <- ev$vectors
  
  # phi = W^{-1/2} U, then normalize in L2_w
  phi <- sweep(U, 1, sw, "/")
  norms <- sqrt(colSums(phi^2 * w))
  phi   <- sweep(phi, 2, norms, "/")
  
  list(lam = lam, phi = phi, w = w)
}

#' Optional subdomain selection: keep t with at least min_frac in both groups
#' @keywords internal
#' @noRd
.tfu_subdomain_idx <- function(O_mat, group_A, min_frac = 0.10) {
  if (is.null(min_frac)) return(seq_len(ncol(O_mat)))
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  which(colSums(O_mat * IA) >= min_frac * n &
          colSums(O_mat * IB) >= min_frac * n)
}

# --------- Algorithm 1: L2 test ---------------------------------------------

#' Algorithm 1: L2 test for equality of mean functions
#'
#' @param X_obs numeric matrix (n x m) of observed values; use NA for missing.
#' @param O_mat optional 0/1 observation matrix (n x m). If NULL (default), it is
#'   constructed internally as 1 where X_obs is observed and 0 where it is NA.
#' @param group_A optional logical/numeric length n; if NULL (default), rows with no NA are
#'   assigned to group A (complete) and rows with any NA to group B (incomplete).
#' @param grid optional numeric vector of length m (strictly increasing); defaults to
#'   seq(0, 1, length.out = ncol(X_obs)).
#' @param fve numeric in (0,1]; fraction of variance explained to truncate KL (default 0.99).
#' @param B integer; Monte-Carlo draws for the null mixture (default 5000).
#' @param eps numeric; lower bound for \eqn{\hat p} (default 1e-8).
#' @param min_frac optional numeric in (0,1]; if given, restricts to t with
#'   at least \code{min_frac} observed in both groups (available-data subdomain).
#' @param seed optional integer for reproducibility; if not NULL, sets seed for simulation.
#' @param delta_A numeric in (0,1]; threshold for grouping rows into A vs B when
#'   \code{group_A} is NULL. Default 1 (Example 1: only fully observed rows in A).
#'   For Example 2, set e.g. 0.7 to put rows with ≥70% observed entries into A.
#' @return list with elements: \code{stat}, \code{p_value}, \code{grid}, \code{muA}, \code{muB},
#'   \code{lam} (used eigenvalues), and \code{idx} (subdomain indices).
#' @export
tfu_algo1_L2_test <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                              fve = 0.99, B = 5000, eps = 1e-8,
                              min_frac = 0.10, seed = NULL, delta_A = 1) {
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  idx <- .tfu_subdomain_idx(O_mat, group_A, min_frac)
  if (length(idx) < 2L) stop("Subdomain too small; relax 'min_frac' or check inputs.")
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  g  <- grid[idx]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  # Observed statistic via tidyfun integration
  diff_tfd <- tf::tfd(matrix(muA - muB, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  # Covariance & KL
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam
  if (!is.null(seed)) set.seed(seed)
  # choose q by FVE
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
  }
  # mixture of chi^2: sum(lam_j Z_j^2)
  Z <- matrix(rnorm(length(lam_q) * B), nrow = length(lam_q))
  W <- colSums((Z^2) * lam_q)
  p  <- (sum(W >= T_L2) + 1) / (B + 1)
  
  list(stat = as.numeric(T_L2), p_value = p,
       grid = g, muA = muA, muB = muB, lam = lam_q, idx = idx)
}

# --------- Algorithm 2: Supremum test ---------------------------------------

#' Algorithm 2: Supremum test (T_{mu,D})
#'
#' @inheritParams tfu_algo1_L2_test
#' @return list with elements: \code{stat}, \code{p_value}, \code{grid}, \code{muA}, \code{muB},
#'   \code{lam}, \code{phi} (first q eigenfunctions), and \code{idx}.
#' @export
tfu_algo2_sup_test <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                               fve = 0.95, B = 5000, eps = 1e-8,
                               min_frac = 0.10, seed = NULL, delta_A = 1) {
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  idx <- .tfu_subdomain_idx(O_mat, group_A, min_frac)
  if (length(idx) < 2L) stop("Subdomain too small; relax 'min_frac' or check inputs.")
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  g  <- grid[idx]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  # Observed statistic
  T_D <- sqrt(n) * max(abs(muA - muB))
  
  # Covariance & KL
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam; phi <- KL$phi
  
  # choose q by FVE
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0; phi_q <- matrix(0, nrow = length(g), ncol = 1L)
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
    phi_q <- phi[, seq_len(q), drop = FALSE]
  }
  
  # Simulate sup | sum_j sqrt(lam_j) Z_j phi_j(t) |
  if (!is.null(seed)) set.seed(seed)
  Z   <- matrix(rnorm(q * B), nrow = q)              # q x B
  lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")       # m x q
  gp_vals <- lam_phi %*% Z                           # m x B
  W <- apply(abs(gp_vals), 2, max)
  p <- (sum(W >= T_D) + 1) / (B + 1)
  
  list(stat = as.numeric(T_D), p_value = p,
       grid = g, muA = muA, muB = muB, lam = lam_q, phi = phi_q, idx = idx)
}

# --------- Algorithm 3: Simultaneous confidence bands ------------------------

#' Algorithm 3: Simultaneous confidence bands for mu_A - mu_B
#'
#' @inheritParams tfu_algo1_L2_test
#' @param alpha numeric in (0,1), confidence level (default 0.05).
#' @return list with elements: \code{grid}, \code{muA}, \code{muB}, \code{diff},
#'   \code{lower}, \code{upper}, \code{q_alpha}, \code{lam}, \code{phi}, \code{idx}.
#' @export
tfu_algo3_conf_bands <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                                 alpha = 0.05, fve = 0.95, B = 5000, eps = 1e-8,
                                 min_frac = 0.10, seed = NULL, delta_A = 1) {
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  idx <- .tfu_subdomain_idx(O_mat, group_A, min_frac)
  if (length(idx) < 2L) stop("Subdomain too small; relax 'min_frac' or check inputs.")
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  g  <- grid[idx]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  # Covariance & KL
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam; phi <- KL$phi
  
  # choose q by FVE
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0; phi_q <- matrix(0, nrow = length(g), ncol = 1L)
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
    phi_q <- phi[, seq_len(q), drop = FALSE]
  }
  
  # simulate sup for quantile
  if (!is.null(seed)) set.seed(seed)
  Z   <- matrix(rnorm(q * B), nrow = q)
  lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")
  gp_vals <- lam_phi %*% Z
  W <- apply(abs(gp_vals), 2, max)
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
  
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth
  upper <- diff + halfwidth
  
  list(grid = g, muA = muA, muB = muB, diff = diff,
       lower = lower, upper = upper, q_alpha = q_alpha,
       lam = lam_q, phi = phi_q, idx = idx)
}
