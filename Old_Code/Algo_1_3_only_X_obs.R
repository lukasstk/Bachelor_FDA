#' Tidyfun extensions: Algorithms 1, 2, and 3 (auto inputs, X_obs-only)
#'
#' Implementiert:
#'  - Algorithmus 1: L2-Test (T_{mu,L2}) mit available-data Schätzern
#'  - Algorithmus 2: Supremum-Test (T_{mu,D}) via KL-basierte GP-Approx.
#'  - Algorithmus 3: Simultane Konfidenzbänder für mu_A - mu_B
#'
#' Änderungen:
#'  - Öffentliche Signaturen erwarten NUR X_obs (numeric Matrix n x m mit NA).
#'  - O_mat, group_A (delta_A=1, automatische Anpassung), grid und Subdomain werden intern gebaut.
#'  - Nur Minimal-Check für X_obs.
#'
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------------------------- interne Utilities (schlank) -----------------------

.tfu_group_from_delta <- function(O_mat, delta) {
  rowMeans(O_mat != 0) >= delta
}

# Subdomain I nach Paper (10%) mit einmaligem Fallback (Warnings)
.tfu_subdomain_idx_paper <- function(O_mat, group_A, min_frac = 0.10, m_min = 2L) {
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  cA <- colSums(O_mat * IA)
  cB <- colSums(O_mat * IB)
  
  idx <- which(cA >= min_frac * n & cB >= min_frac * n)
  if (length(idx) >= max(2L, m_min)) {
    return(list(idx = idx, min_frac_used = min_frac, fallback = NULL))
  }
  
  idx2 <- which(cA > 0 & cB > 0)
  if (length(idx2) >= max(2L, m_min)) {
    warning("Subdomain (10%): keine Punkte mit ≥ min_frac in beiden Gruppen -> verwende gemeinsame Überlappung (>0).")
    return(list(idx = idx2, min_frac_used = NA_real_, fallback = "overlap"))
  }
  
  idx3 <- which(colSums(O_mat) > 0)
  if (length(idx3) >= max(2L, m_min)) {
    warning("Subdomain: keine gemeinsame Überlappung -> verwende alle beobachteten Zeitpunkte (evtl. geringe Aussagekraft).")
    return(list(idx = idx3, min_frac_used = NA_real_, fallback = "any-observed"))
  }
  
  warning("Subdomain: zu wenige Beobachtungen insgesamt. Ergebnisse nicht interpretierbar (trivialer Fallback).")
  list(idx = seq_len(min(2L, ncol(O_mat))), min_frac_used = NA_real_, fallback = "trivial")
}

# Baut O_mat, group_A, grid aus X_obs (nur X_obs-Check)
.tfu_prepare_inputs <- function(X_obs) {
  if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
  storage.mode(X_obs) <- "numeric"
  n <- nrow(X_obs); m <- ncol(X_obs)
  if (m < 2L) stop("X_obs muss >= 2 Spalten (Zeitpunkte) haben.")
  
  O_mat <- 1L * (!is.na(X_obs))
  
  # automatische Gruppierung mit Start delta_A=1 und einmaliger Lockerung falls nötig
  group_A <- .tfu_group_from_delta(O_mat, delta = 1)
  nA <- sum(group_A); nB <- sum(!group_A)
  if (nA == 0L || nB == 0L) {
    obs_frac <- rowMeans(O_mat != 0)
    thr_seq <- sort(unique(c(seq(1, 0.5, by = -0.05), 0.49, 0.4, 0.3, 0.2, 0.1)), decreasing = TRUE)
    for (thr in thr_seq) {
      cand <- .tfu_group_from_delta(O_mat, thr)
      if (sum(cand) > 0L && sum(!cand) > 0L) { group_A <- cand; break }
    }
    if (sum(group_A) == 0L || sum(!group_A) == 0L) {
      medf <- median(obs_frac)
      group_A <- obs_frac >= medf
      warning("Automatische Gruppierung: Median-Split, da zuvor eine leere Gruppe vorlag.")
    } else {
      warning("Automatische Gruppierung: delta_A intern reduziert, um beide Gruppen zu füllen.")
    }
  }
  
  grid <- seq(0, 1, length.out = m)
  list(X_obs = X_obs, O_mat = O_mat, group_A = group_A, grid = grid)
}

# Trapezgewichte für evtl. nicht-uniformes Grid
.tfu_trap_weights <- function(x) {
  m <- length(x)
  if (m == 1L) return(1)
  w <- numeric(m)
  w[1] <- (x[2] - x[1]) / 2
  w[m] <- (x[m] - x[m-1]) / 2
  if (m > 2L) w[2:(m-1)] <- (x[3:m] - x[1:(m-2)]) / 2
  w
}

# Available-Data Mittelwerte & Beobachtungswahrsch. (denom: n)
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

# Korrekturschätzer der Kovarianz
.tfu_corrected_cov <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X_obs)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  Xc <- X_obs
  A_idx <- which(IA == 1)
  B_idx <- which(IB == 1)
  if (length(A_idx)) Xc[A_idx, ] <- sweep(X_obs[A_idx, , drop = FALSE], 2, muA_hat, `-`)
  if (length(B_idx)) Xc[B_idx, ] <- sweep(X_obs[B_idx, , drop = FALSE], 2, muB_hat, `-`)
  
  Xc[O_mat == 0] <- 0
  Xc[is.na(Xc)]  <- 0
  
  XA <- Xc[A_idx, , drop = FALSE]
  XB <- Xc[B_idx, , drop = FALSE]
  
  if (nrow(XA)) XA <- sweep(XA, 2, pA_hat, "/")
  if (nrow(XB)) XB <- sweep(XB, 2, pB_hat, "/")
  
  Sum_A <- if (nrow(XA)) crossprod(XA) else matrix(0, ncol(X_obs), ncol(X_obs))
  Sum_B <- if (nrow(XB)) crossprod(XB) else matrix(0, ncol(X_obs), ncol(X_obs))
  
  K <- (Sum_A + Sum_B) / n
  (K + t(K)) / 2
}

# KL-Basis aus Kovarianz mit Trapez-Gewichten
.tfu_kl_from_cov <- function(K, grid) {
  w <- .tfu_trap_weights(grid)
  sw <- sqrt(w)
  S <- (sw * t(sw * K))          # entspricht diag(sw) %*% K %*% diag(sw)
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  U   <- ev$vectors
  phi <- sweep(U, 1, sw, "/")    # W^{-1/2} U
  norms <- sqrt(colSums(phi^2 * w))
  phi   <- sweep(phi, 2, norms, "/")  # L2_w-orthonormal
  list(lam = lam, phi = phi, w = w)
}

# ----------------------------- Algorithmus 1: L2-Test --------------------------

#' L2-Test (nur X_obs)
#' @return Liste: stat, p_value, grid, muA, muB, lam, idx, min_frac_used, fallback, delta_A
tfu_algo1_L2_test <- function(X_obs) {
  # feste Defaults (bei Bedarf im Code anpassen)
  fve <- 0.99; B <- 5000L; eps <- 1e-8; min_frac <- 0.10
  delta_A <- 1
  
  prep <- .tfu_prepare_inputs(X_obs)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g  <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  diff_tfd <- tf::tfd(matrix(muA - muB, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam
  
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
  }
  Z <- matrix(rnorm(length(lam_q) * B), nrow = length(lam_q))
  W <- colSums((Z^2) * lam_q)
  p  <- (sum(W >= T_L2) + 1) / (B + 1)
  
  list(stat = as.numeric(T_L2), p_value = p,
       grid = g, muA = muA, muB = muB, lam = lam_q, idx = idx,
       min_frac_used = sub$min_frac_used, fallback = sub$fallback, delta_A = delta_A)
}

# ----------------------------- Algorithmus 2: Supremum-Test --------------------

#' Supremum-Test (nur X_obs)
#' @return Liste: stat, p_value, grid, muA, muB, lam, phi, idx, min_frac_used, fallback, delta_A
tfu_algo2_sup_test <- function(X_obs) {
  # feste Defaults
  fve <- 0.95; B <- 5000L; eps <- 1e-8; min_frac <- 0.10
  delta_A <- 1
  
  prep <- .tfu_prepare_inputs(X_obs)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g  <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  T_D <- sqrt(n) * max(abs(muA - muB))
  
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam; phi <- KL$phi
  
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0; phi_q <- matrix(0, nrow = length(g), ncol = 1L)
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
    phi_q <- phi[, seq_len(q), drop = FALSE]
  }
  
  Z   <- matrix(rnorm(q * B), nrow = q)           # q x B
  lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")    # m x q
  gp_vals <- lam_phi %*% Z                        # m x B
  W <- apply(abs(gp_vals), 2, max)
  p <- (sum(W >= T_D) + 1) / (B + 1)
  
  list(stat = as.numeric(T_D), p_value = p,
       grid = g, muA = muA, muB = muB, lam = lam_q, phi = phi_q, idx = idx,
       min_frac_used = sub$min_frac_used, fallback = sub$fallback, delta_A = delta_A)
}

# ----------------------------- Algorithmus 3: Konfidenzbänder ------------------

#' Simultane Konfidenzbänder (nur X_obs) – nutzt die KL/Eigenbasis aus tfu_algo2_sup_test
#' @return Liste: grid, muA, muB, diff, lower, upper, q_alpha, lam, phi, idx, min_frac_used, fallback, delta_A
tfu_algo3_conf_bands <- function(X_obs) {
  alpha <- 0.05; B <- 5000L
  X_mat <- if (!is.matrix(X_obs)) as.matrix(X_obs) else X_obs
  n <- nrow(X_mat)
  
  a2 <- tfu_algo2_sup_test(X_obs)
  g <- a2$grid
  muA <- a2$muA
  muB <- a2$muB
  lam_q <- a2$lam
  phi_q <- a2$phi
  
  diff <- muA - muB
  q <- length(lam_q)
  
  Z <- matrix(rnorm(q * B), nrow = q)
  lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")
  gp_vals <- lam_phi %*% Z
  W <- apply(abs(gp_vals), 2, max)
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
  
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth
  upper <- diff + halfwidth
  
  list(
    grid = g, muA = muA, muB = muB, diff = diff,
    lower = lower, upper = upper, q_alpha = q_alpha,
    lam = lam_q, phi = phi_q, idx = a2$idx,
    min_frac_used = a2$min_frac_used, fallback = a2$fallback, delta_A = a2$delta_A
  )
}
