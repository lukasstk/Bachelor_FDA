# -----------------------------------------------------------------------------
# Tidyfun extension: Algorithm 5 + Bootstrap Confidence Bands (Algorithm 7)
# -----------------------------------------------------------------------------
#
# This file extends the previous `tfu_algo5_bootstrap()` by optionally producing
# bootstrap simultaneous confidence bands for mu_A - mu_B (Algorithm 7 in the
# paper). In addition, a convenience function `tfu_algo7_conf_bands_bootstrap()`
# is provided to compute only the bands.
#
# Assumptions:
# - The helper utilities `.tfu_prepare_inputs()`, `.tfu_subdomain_idx_paper()`,
#   `.tfu_available_means()`, and `.tfu_boot_group_mean()` from your existing
#   code are available in the same package/file set.
#   If `.tfu_boot_group_mean()` is not defined yet, an implementation is
#   included below (guarded to avoid re-definition).
#
# -----------------------------------------------------------------------------

# define only if not already defined in the session/package
if (!exists(".tfu_boot_group_mean", mode = "function")) {
  # group-wise available-data mean inside bootstrap (X_grp may contain NAs; O_grp same shape)
  .tfu_boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
    p_hat <- colSums(O_grp) / n_total
    p_hat <- pmax(p_hat, eps)
    numer <- colSums(replace(X_grp, is.na(X_grp), 0))
    numer / (n_total * p_hat)
  }
}

# -----------------------------------------------------------------------------
# Extended Algorithm 5: add bootstrap confidence bands (Algorithm 7 style)
# -----------------------------------------------------------------------------

#' Algorithm 5 (extended): Bootstrap p-values for L2 and Supremum tests
#'                        + simultaneous confidence bands (Algorithm 7)
#'
#' @param X_obs numeric matrix (n x m) of observed values; use NA for missing.
#' @param O_mat optional 0/1 observation matrix (n x m). If NULL, built from X_obs.
#' @param group_A optional logical/numeric length n. If NULL, rows with observed
#'        fraction ≥ delta_A are assigned to group A (Example 1 corresponds to 1).
#' @param grid optional numeric length m; if NULL, seq(0, 1, length.out = m).
#' @param B integer; number of bootstrap draws (default 5000).
#' @param eps numeric; lower bound for p-hat (default 1e-8).
#' @param min_frac numeric in (0,1]; Paper subdomain (default 0.10).
#' @param return_boot logical; if TRUE, also return bootstrap distributions.
#' @param seed optional integer for reproducibility.
#' @param delta_A numeric in (0,1]; threshold for grouping A vs. B when group_A is NULL.
#' @param compute_bands logical; if TRUE, compute bootstrap simultaneous bands
#'        for mu_A - mu_B following Algorithm 7 (default TRUE).
#' @param alpha numeric; 1 - confidence level for the bands (default 0.05).
#' @return list with elements:
#'   - grid, muA, muB, diff
#'   - T_L2, T_D, p_L2, p_D
#'   - (if compute_bands) lower, upper, q_alpha
#'   - (if return_boot) boot_L2, boot_D
#'   - idx, min_frac_used, fallback, delta_A
#' @export
tfu_algo5_bootstrap <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                                B = 5000, eps = 1e-8, min_frac = 0.10,
                                return_boot = FALSE, seed = NULL, delta_A = 1,
                                compute_bands = TRUE, alpha = 0.05) {
  # --- prepare inputs & subdomain ---
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  # --- available-data means ---
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  # --- observed statistics ---
  # L2 via tidyfun integrate (trapz) and D via sup
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2 <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  T_D  <- sqrt(n) * max(abs(diff))
  
  # --- build centered residuals by group (keep NAs; handled by O in means) ---
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  # --- bootstrap ---
  if (!is.null(seed)) set.seed(seed)
  nA <- length(A_idx); nB <- length(B_idx)
  boot_L2 <- numeric(B); boot_D <- numeric(B)
  
  for (b in seq_len(B)) {
    sampA <- if (nA) sample(A_idx, nA, replace = TRUE) else integer(0)
    sampB <- if (nB) sample(B_idx, nB, replace = TRUE) else integer(0)
    
    muA_b <- if (nA) .tfu_boot_group_mean(X_cent[sampA, , drop = FALSE],
                                          O[sampA, , drop = FALSE], n, eps) else rep(0, ncol(X))
    muB_b <- if (nB) .tfu_boot_group_mean(X_cent[sampB, , drop = FALSE],
                                          O[sampB, , drop = FALSE], n, eps) else rep(0, ncol(X))
    
    d_b <- muA_b - muB_b
    d_tfd <- tf::tfd(matrix(d_b, nrow = 1), arg = g)
    boot_L2[b] <- n * tf::tf_integrate(d_tfd^2, arg = g)
    boot_D[b]  <- sqrt(n) * max(abs(d_b))
  }
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (B + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (B + 1)
  
  out <- list(grid = g, muA = muA, muB = muB, diff = diff,
              T_L2 = as.numeric(T_L2), T_D = as.numeric(T_D),
              p_L2 = p_L2, p_D = p_D, idx = idx,
              min_frac_used = sub$min_frac_used, fallback = sub$fallback,
              delta_A = delta_A)
  
  # --- Algorithm 7 bands ---
  if (isTRUE(compute_bands)) {
    q_alpha <- as.numeric(stats::quantile(boot_D, probs = 1 - alpha, names = FALSE))
    halfwidth <- q_alpha / sqrt(n)
    out$lower <- diff - halfwidth
    out$upper <- diff + halfwidth
    out$q_alpha <- q_alpha
    out$alpha <- alpha
  }
  
  if (isTRUE(return_boot)) {
    out$boot_L2 <- boot_L2
    out$boot_D  <- boot_D
  }
  out
}

# -----------------------------------------------------------------------------
# Dedicated convenience wrapper: ONLY bootstrap confidence bands (Algorithm 7)
# -----------------------------------------------------------------------------

#' Algorithm 7: Bootstrap simultaneous confidence bands for mu_A - mu_B
#'
#' Thin wrapper that reuses `tfu_algo5_bootstrap()` to compute only the bands.
#'
#' @inheritParams tfu_algo5_bootstrap
#' @return list with elements: grid, muA, muB, diff, lower, upper, q_alpha,
#'         idx, min_frac_used, fallback, delta_A, alpha. Optionally also p-values
#'         and bootstrap distributions if requested via `with_tests` / `return_boot`.
#' @param with_tests logical; if TRUE, also return L2 / Sup-test stats & p-values
#'        (they are computed anyway). Default FALSE keeps the output lean.
#' @export
tfu_algo7_conf_bands_bootstrap <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                                           B = 10000, eps = 1e-8, min_frac = 0.10,
                                           seed = NULL, delta_A = 1, alpha = 0.05,
                                           return_boot = FALSE, with_tests = FALSE) {
  res <- tfu_algo5_bootstrap(X_obs = X_obs, O_mat = O_mat, group_A = group_A, grid = grid,
                             B = B, eps = eps, min_frac = min_frac, seed = seed,
                             delta_A = delta_A, return_boot = return_boot,
                             compute_bands = TRUE, alpha = alpha)
  if (!isTRUE(with_tests)) {
    # keep only the band-related essentials
    keep <- c("grid","muA","muB","diff","lower","upper","q_alpha",
              "idx","min_frac_used","fallback","delta_A","alpha")
    return(res[keep])
  }
  res
}