#' Tidyfun extensions: Algorithms 1, 2, and 3 (auto inputs)
#'
#' Implements:
#'  - Algorithm 1: L2 test of equality of mean functions with available-data estimators
#'  - Algorithm 2: Supremum test (T_{mu,D}) with KL-based Gaussian process approximation
#'  - Algorithm 3: Simultaneous confidence bands for mu_A - mu_B
#'
#' Änderungen vs. Original: O_mat, group_A und grid werden **intern** aus X_obs erzeugt.
#'  - O_mat := 1, wenn X_obs beobachtet (nicht-NA), sonst 0.
#'  - group_A (Bsp. 1/2 in Ofner et al., 2025): Zeile in A, wenn beobachteter Anteil ≥ delta_A;
#'    Default delta_A = 1 (nur vollständig beobachtete in A). Für Bsp. 2 z.B. delta_A = 0.7.
#'  - grid := seq(0, 1, length.out = ncol(X_obs)).
#'  - Subdomain **Paper-Regel**: wählt direkt die Punkte mit mind. 10% Beobachtung in **beiden**
#'    Gruppen. Falls leer, einmaliger Fallback auf gemeinsame Überlappung (>0), sonst alle
#'    beobachteten Punkte, sonst minimal 2 Punkte. In diesen Fällen gibt es **Warnings**.
#'
#' Hinweis: Integration nutzt tidyfun::tf_integrate für numerische Stabilität.
#'
#' @name tidyfun_ext_algo1to3_auto
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------------------------- interne Utilities ---------------------------------

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

#' Subdomain I nach Paper (10%) mit **einmaligem** Fallback (Warnings)
#' @keywords internal
#' @noRd
.tfu_subdomain_idx_paper <- function(O_mat, group_A, min_frac = 0.10, m_min = 2L) {
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  cA <- colSums(O_mat * IA)
  cB <- colSums(O_mat * IB)
  
  # Paper-Regel: ≥ min_frac * n in beiden Gruppen
  idx <- which(cA >= min_frac * n & cB >= min_frac * n)
  if (length(idx) >= max(2L, m_min)) {
    return(list(idx = idx, min_frac_used = min_frac, fallback = NULL))
  }
  
  # Fallback 1: gemeinsame Überlappung (>0 in beiden Gruppen)
  idx2 <- which(cA > 0 & cB > 0)
  if (length(idx2) >= max(2L, m_min)) {
    warning("Subdomain (10%): keine Punkte mit ≥ min_frac in beiden Gruppen -> verwende gemeinsame Überlappung (>0).")
    return(list(idx = idx2, min_frac_used = NA_real_, fallback = "overlap"))
  }
  
  # Fallback 2: alle beobachteten Zeitpunkte (geringe Aussagekraft)
  idx3 <- which(colSums(O_mat) > 0)
  if (length(idx3) >= max(2L, m_min)) {
    warning("Subdomain: keine gemeinsame Überlappung -> verwende alle beobachteten Zeitpunkte (evtl. geringe Aussagekraft).")
    return(list(idx = idx3, min_frac_used = NA_real_, fallback = "any-observed"))
  }
  
  # Härtester Fallback: nimm die ersten 2 Spalten
  warning("Subdomain: zu wenige Beobachtungen insgesamt. Ergebnisse nicht interpretierbar (trivialer Fallback).")
  list(idx = seq_len(min(2L, ncol(O_mat))), min_frac_used = NA_real_, fallback = "trivial")
}

#' Build/validate O_mat, group_A, grid from X_obs when missing
#' @keywords internal
#' @noRd
.tfu_prepare_inputs <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL, delta_A = 1) {
  # ensure numeric matrix
  if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
  storage.mode(X_obs) <- "numeric"
  n <- nrow(X_obs); m <- ncol(X_obs)
  if (m < 2L) stop("X_obs muss mindestens 2 Spalten (Zeitpunkte) besitzen.")
  
  # O_mat: 1 wenn beobachtet, 0 wenn NA
  if (is.null(O_mat)) {
    O_mat <- 1L * (!is.na(X_obs))
  } else {
    stopifnot(is.matrix(O_mat), nrow(O_mat) == n, ncol(O_mat) == m)
    O_mat <- 1L * (O_mat != 0)
  }
  
  # group_A: Example 1/2 via delta_A
  if (is.null(group_A)) {
    group_A <- .tfu_group_from_delta(O_mat, delta_A)
  } else {
    stopifnot(length(group_A) == n)
    group_A <- as.logical(group_A)
  }
  
  # Warnen, falls eine Gruppe leer ist -> automatisch Delta lockern/verschärfen (einmalig)
  nA <- sum(group_A); nB <- sum(!group_A)
  if (nA == 0L || nB == 0L) {
    obs_frac <- rowMeans(O_mat != 0)
    if (nA == 0L) {
      # Delta automatisch reduzieren, bis beide Gruppen belegt sind
      thr_seq <- sort(unique(c(seq(delta_A, 0.5, by = -0.05), 0.49, 0.4, 0.3, 0.2, 0.1)), decreasing = TRUE)
      for (thr in thr_seq) {
        cand <- .tfu_group_from_delta(O_mat, thr)
        if (sum(cand) > 0L && sum(!cand) > 0L) { group_A <- cand; break }
      }
      if (sum(group_A) == 0L || sum(!group_A) == 0L) {
        # Fallback: Median-Split
        medf <- median(obs_frac)
        group_A <- obs_frac >= medf
        warning("Automatische Gruppierung: delta_A angepasst (Median-Split), da zuvor eine leere Gruppe vorlag.")
      } else {
        warning("Automatische Gruppierung: delta_A reduziert, um beide Gruppen zu füllen.")
      }
    } else {
      # Analog: delta erhöhen – hier direkt Median-Split
      medf <- median(obs_frac)
      group_A <- obs_frac >= medf
      warning("Automatische Gruppierung: delta_A erhöht (Median-Split), da Gruppe B leer war.")
    }
  }
  
  # grid
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = m)
  } else {
    stopifnot(is.numeric(grid), length(grid) == m)
  }
  
  .tfu_assert_inputs(X_obs, O_mat, group_A, grid)
  list(X_obs = X_obs, O_mat = O_mat, group_A = group_A, grid = grid)
}

#' Trapezgewichte auf ggf. nicht-uniformem Grid
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

#' Available-Data Mittelwerte & Beobachtungswahrsch. (denom: n)
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

#' Korrekturschätzer der Kovarianz K(s,t)
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
  
  # unobserved/NA auf 0
  Xc[O_mat == 0] <- 0
  Xc[is.na(Xc)]  <- 0
  
  # Spaltung nach Gruppen
  XA <- Xc[A_idx, , drop = FALSE]
  XB <- Xc[B_idx, , drop = FALSE]
  
  # spaltenweise Skalierung
  if (nrow(XA)) XA <- sweep(XA, 2, pA_hat, "/")
  if (nrow(XB)) XB <- sweep(XB, 2, pB_hat, "/")
  
  Sum_A <- if (nrow(XA)) crossprod(XA) else matrix(0, ncol(X_obs), ncol(X_obs))
  Sum_B <- if (nrow(XB)) crossprod(XB) else matrix(0, ncol(X_obs), ncol(X_obs))
  
  K <- (Sum_A + Sum_B) / n
  K <- (K + t(K)) / 2  # symmetrisieren
  K
}

#' KL-Basis aus Kovarianz mit Trapez-Gewichten (W^{1/2} K W^{1/2})
#' @keywords internal
#' @noRd
.tfu_kl_from_cov <- function(K, grid) {
  w <- .tfu_trap_weights(grid)
  sw <- sqrt(w)
  S <- (sw * t(sw * K))              # diag(sw) %*% K %*% diag(sw)
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  U   <- ev$vectors
  phi <- sweep(U, 1, sw, "/")       # W^{-1/2} U
  norms <- sqrt(colSums(phi^2 * w))
  phi   <- sweep(phi, 2, norms, "/")  # L2_w-orthonormal
  list(lam = lam, phi = phi, w = w)
}

# ----------------------------- Algorithmus 1: L2-Test --------------------------

#' Algorithmus 1: L2-Test auf Gleichheit der Mittelwerte (verfügbare Daten)
#'
#' @param X_obs numeric Matrix (n x m) mit Beobachtungen; NA für fehlend.
#' @param O_mat optional 0/1-Matrix (n x m). Default: intern aus X_obs.
#' @param group_A optional logischer Vektor Länge n; Default: intern via delta_A.
#' @param grid optional numeric Vektor Länge m; Default: seq(0,1,len=m).
#' @param fve Anteil erklärter Varianz für KL-Trunkierung (Default 0.99).
#' @param B  Monte-Carlo Ziehungen für Mischungsapproximation (Default 5000).
#' @param eps Untergrenze für p̂ (Default 1e-8).
#' @param min_frac numeric in [0,1]; **Paper-Regel** Mindestanteil je Gruppe für Subdomain (Default 0.10).
#' @param seed integer für Reproduzierbarkeit.
#' @param delta_A Schwelle in (0,1]; Gruppierung A vs. B, wenn group_A fehlend (Default 1).
#' @return Liste: stat, p_value, grid, muA, muB, lam, idx, min_frac_used, fallback, delta_A
#' @export
tfu_algo1_L2_test <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                              fve = 0.99, B = 5000, eps = 1e-8,
                              min_frac = 0.10, seed = NULL, delta_A = 1) {
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g  <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  # beobachtete Statistik
  diff_tfd <- tf::tfd(matrix(muA - muB, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  # Kovarianz & KL
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam
  if (!is.null(seed)) set.seed(seed)
  # q via FVE
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

#' Algorithmus 2: Supremum-Test (T_{mu,D})
#'
#' @inheritParams tfu_algo1_L2_test
#' @return Liste: stat, p_value, grid, muA, muB, lam, phi, idx, min_frac_used, fallback, delta_A
#' @export
tfu_algo2_sup_test <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                               fve = 0.95, B = 5000, eps = 1e-8,
                               min_frac = 0.10, seed = NULL, delta_A = 1) {
  prep <- .tfu_prepare_inputs(X_obs, O_mat, group_A, grid, delta_A = delta_A)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g  <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]
  O  <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  
  # beobachtete Statistik
  T_D <- sqrt(n) * max(abs(muA - muB))
  
  # Kovarianz & KL
  K   <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL  <- .tfu_kl_from_cov(K, g)
  lam <- KL$lam; phi <- KL$phi
  
  # q via FVE
  if (sum(lam) <= 0) {
    q <- 1L; lam_q <- 0; phi_q <- matrix(0, nrow = length(g), ncol = 1L)
  } else {
    cum <- cumsum(lam) / sum(lam)
    q <- which(cum >= fve)[1]
    lam_q <- lam[seq_len(q)]
    phi_q <- phi[, seq_len(q), drop = FALSE]
  }
  
  if (!is.null(seed)) set.seed(seed)
  Z   <- matrix(rnorm(q * B), nrow = q)              # q x B
  lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")       # m x q
  gp_vals <- lam_phi %*% Z                           # m x B
  W <- apply(abs(gp_vals), 2, max)
  p <- (sum(W >= T_D) + 1) / (B + 1)
  
  list(stat = as.numeric(T_D), p_value = p,
       grid = g, muA = muA, muB = muB, lam = lam_q, phi = phi_q, idx = idx,
       min_frac_used = sub$min_frac_used, fallback = sub$fallback, delta_A = delta_A)
}

# ----------------------------- Algorithmus 3: Konfidenzbänder ------------------

#' Algorithmus 3: Simultane Konfidenzbänder für mu_A - mu_B
#'  -> nutzt tfu_algo2_sup_test für KL-Eigenwerte/-funktionen und Supremum-Approximation
#'
#' @inheritParams tfu_algo1_L2_test
#' @param alpha Konfidenzniveau (Default 0.05)
#' @return Liste: grid, muA, muB, diff, lower, upper, q_alpha, lam, phi, idx, min_frac_used, fallback, delta_A
#' @export
tfu_algo3_conf_bands <- function(X_obs, O_mat = NULL, group_A = NULL, grid = NULL,
                                 alpha = 0.05, fve = 0.95, B = 5000, eps = 1e-8,
                                 min_frac = 0.10, seed = NULL, delta_A = 1) {
  # Hole alle benötigten Größen direkt aus Algorithmus 2 (KL-Basis & Trunkierung bereits dort bestimmt)
  a2 <- tfu_algo2_sup_test(
    X_obs = X_obs, O_mat = O_mat, group_A = group_A, grid = grid,
    fve = fve, B = B, eps = eps, min_frac = min_frac, seed = seed, delta_A = delta_A
  )
  
  # Entpacken
  g    <- a2$grid
  muA  <- a2$muA
  muB  <- a2$muB
  lam  <- a2$lam      # bereits nach FVE getrimmt
  phi  <- a2$phi      # dito
  idx  <- a2$idx
  diff <- muA - muB
  
  # Monte-Carlo für Supremum der Grenz-GP nur mit (lam, phi) aus Algo 2
  if (!is.null(seed)) set.seed(seed)
  q <- length(lam)
  if (q < 1) {
    lam <- 0
    phi <- matrix(0, nrow = length(g), ncol = 1L)
    q <- 1L
  }
  Z <- matrix(rnorm(q * B), nrow = q)            # q x B
  lam_phi <- sweep(phi, 2, sqrt(lam), "*")       # m x q
  gp_vals <- lam_phi %*% Z                       # m x B
  W <- apply(abs(gp_vals), 2, max)
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
  
  # Bandbreite und Bänder
  n <- nrow(as.matrix(X_obs))
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth
  upper <- diff + halfwidth
  
  list(
    grid = g, muA = muA, muB = muB, diff = diff,
    lower = lower, upper = upper, q_alpha = q_alpha,
    lam = lam, phi = phi, idx = idx,
    min_frac_used = a2$min_frac_used, fallback = a2$fallback, delta_A = a2$delta_A
  )
}
