#' Tidyfun extension: Algorithm 5 (bootstrap for L2 and Supremum) — auto inputs
#'
#' Implements Algorithm 5 with automatic construction of O_mat, group_A, and grid.
#'  - O_mat := 1 if X_obs is observed (non-NA) and 0 otherwise.
#'  - group_A (Example 1 / Example 2): by default rows with observed fraction ≥ delta_A
#'    are assigned to A, the rest to B. With delta_A = 1 this is the "complete vs. incomplete"
#'    split (Example 1). For Example 2 choose delta_A in (0,1), e.g. 0.7.
#'  - grid := seq(0, 1, length.out = ncol(X_obs)).
#'
#' You may still override any of these by passing the argument explicitly.
#'
#' @name tidyfun_ext_algo5
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# ---------------- internal utilities ----------------

# accept matrix or data.frame for X_obs
.tfu_assert_inputs <- function(X_obs, O_mat, group_A, grid) {
  stopifnot(is.matrix(X_obs) || is.data.frame(X_obs))
  if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
  stopifnot(is.matrix(O_mat), is.numeric(grid))
  stopifnot(nrow(X_obs) == nrow(O_mat), ncol(X_obs) == ncol(O_mat))
  stopifnot(length(group_A) == nrow(X_obs))
  stopifnot(all(is.finite(grid)), length(unique(grid)) == length(grid))
}

# grouping helper for Example 1/2
.tfu_group_from_delta <- function(O_mat, delta) {
  stopifnot(is.numeric(delta), length(delta) == 1L, delta > 0, delta <= 1)
  rowMeans(O_mat != 0) >= delta
}

# build/validate O_mat, group_A, grid from X_obs when missing
.tfu_prepare_inputs <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL, delta_A = 1) {
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

.tfu_subdomain_idx <- function(O_mat, group_A, min_frac = 0.10) {
  if (is.null(min_frac)) return(seq_len(ncol(O_mat)))
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  which(colSums(O_mat * IA) >= min_frac * n &
          colSums(O_mat * IB) >= min_frac * n)
}

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

.tfu_boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
  p_hat <- colSums(O_grp) / n_total
  p_hat <- pmax(p_hat, eps)
  numer <- colSums(replace(X_grp, is.na(X_grp), 0))
  numer / (n_total * p_hat)
}

# --------------- Algorithm 5 (bootstrap) ---------------

#' Algorithm 5: Bootstrap p-values for L2 and Supremum tests (auto inputs)
#'
#' @param X_obs numeric matrix (n x m) of observed values; use NA for missing.
#' @param O_mat optional 0/1 observation matrix (n x m). If NULL, built from X_obs.
#' @param group_A optional logical/numeric length n. If NULL, rows with observed fraction
#'   ≥ \code{delta_A} are assigned to group A (Example 1 corresponds to \code{delta_A = 1}).
#' @param grid optional numeric length m; if NULL, seq(0, 1, length.out = ncol(X_obs)).
#' @param B integer; number of bootstrap draws (default 5000).
#' @param eps numeric; lower bound for \\eqn{\\hat p} (default 1e-8).
#' @param min_frac optional numeric in (0,1]; if given, restricts to t with
#'   at least \\code{min_frac} observed in both groups (available-data subdomain).
#' @param return_boot logical; if TRUE, also return bootstrap distributions.
#' @param seed optional integer for reproducibility; if not NULL, sets seed.
#' @param delta_A numeric in (0,1]; threshold for grouping rows into A vs B when
#'   \code{group_A} is NULL. Default 1 (Example 1). For Example 2, set e.g. 0.7.
#' @return list with elements:
#'   \\itemize{
#'     \\item \\code{grid}, \\code{muA}, \\code{muB}
#'     \\item \\code{T_L2}, \\code{T_D} (observed statistics)
#'     \\item \\code{p_L2}, \\code{p_D} (bootstrap p-values)
#'     \\item \\code{boot_L2}, \\code{boot_D} (if \\code{return_boot=TRUE})
#'     \\item \\code{idx} (subdomain indices)
#'   }
#' @export
tfu_algo5_bootstrap <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                                B = 5000, eps = 1e-8, min_frac = 0.10,
                                return_boot = FALSE, seed = NULL, delta_A = 1) {
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
  
  # Observed statistics
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2 <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  T_D  <- sqrt(n) * max(abs(diff))
  
  # Center residuals by group
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  # Bootstrap
  if (!is.null(seed)) set.seed(seed)
  nA <- length(A_idx); nB <- length(B_idx)
  boot_L2 <- numeric(B); boot_D <- numeric(B)
  
  for (b in seq_len(B)) {
    sampA <- if (nA) sample(A_idx, nA, replace = TRUE) else integer(0)
    sampB <- if (nB) sample(B_idx, nB, replace = TRUE) else integer(0)
    
    muA_b <- if (nA) .tfu_boot_group_mean(X_cent[sampA, , drop = FALSE], O[sampA, , drop = FALSE], n, eps) else rep(0, ncol(X))
    muB_b <- if (nB) .tfu_boot_group_mean(X_cent[sampB, , drop = FALSE], O[sampB, , drop = FALSE], n, eps) else rep(0, ncol(X))
    
    d_b <- muA_b - muB_b
    d_tfd <- tf::tfd(matrix(d_b, nrow = 1), arg = g)
    boot_L2[b] <- n * tf::tf_integrate(d_tfd^2, arg = g)
    boot_D[b]  <- sqrt(n) * max(abs(d_b))
  }
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (B + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (B + 1)
  
  out <- list(grid = g, muA = muA, muB = muB,
              T_L2 = as.numeric(T_L2), T_D = as.numeric(T_D),
              p_L2 = p_L2, p_D = p_D, idx = idx)
  if (isTRUE(return_boot)) {
    out$boot_L2 <- boot_L2
    out$boot_D  <- boot_D
  }
  out
}
