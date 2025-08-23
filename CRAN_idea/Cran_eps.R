#' Tidyfun extensions: Mean-based MCAR tests & bands (auto inputs; fd-or-matrix)
#'
#' Implements:
#'  - asym_mean_L2_test(): L2 test (asymptotic, KL-mixture, optional parallel MC)
#'  - asym_mean_sup_test(): Supremum test (asymptotic, KL-based GP approx, optional parallel MC)
#'  - asym_conf_bands(): Simultaneous confidence bands (asymptotic, optional parallel MC)
#'  - boot_mean_test(): Bootstrap-based p-values for L2/Supremum (foreach + doRNG)
#'  - boot_conf_bands(): Bootstrap simultaneous confidence bands
#'
#' Erweiterung:
#'  - Rückgabeobjekte sind erweiterte `htest`-Listen mit estimate, null.value,
#'    conf.int, alternative – im Stil von t.test().
#'
#' Numerik-/Stabilitäts-Update:
#'  - Einheitlicher `eps`-Parameter (Default 1e-8) zur Stabilisierung von Divisionen
#'    durch kleine p(t) sowie zur (optionalen) Ridge/Jitter‑Regularisierung in der
#'    KL‑Zerlegung (.tfu_kl_from_cov). Über alle öffentlichen Funktionen propagiert.
#'
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------------------------- Utilities -----------------------------------------

#' @title Convert group labels to logical
#' @description Interne Hilfsfunktion: konvertiert verschiedene Gruppencodierungen
#'   zu einem logischen Vektor (A = TRUE, B = FALSE).
#'   Achtung: welche Gruppe "vollständiger" ist, wird erst in
#'   `.tfu_prepare_inputs()` überprüft und ggf. getauscht.
#' @keywords internal
#' @param groups Gruppierungsvektor (logical, factor, character oder numeric/integer).
#'   Muss genau zwei Ausprägungen enthalten (außer bei logical).
#'   - Bei `logical`: wird direkt übernommen (`TRUE` = Gruppe A, `FALSE` = Gruppe B).
#'   - Bei `factor`: der erste Level nach `droplevels()` wird als Gruppe A genommen.
#'   - Bei `character`/`numeric`/`integer`: die erste im Vektor auftretende Ausprägung
#'     wird als Gruppe A gesetzt.
#'   Die Länge muss `nrow(X_obs)` entsprechen und darf keine fehlenden Werte (`NA`) enthalten.
#' @return Logischer Vektor gleicher Länge: TRUE für Gruppe A, FALSE für Gruppe B.
.tfu_groups_to_logical <- function(groups) {
  if (is.factor(groups)) groups <- droplevels(groups)
  if (is.logical(groups)) return(groups)
  if (is.character(groups) || is.numeric(groups) || is.integer(groups)) {
    u <- unique(groups)
    if (length(u) != 2L) stop("`groups` muss genau 2 verschiedene Werte enthalten.")
    return(groups == u[1])
  }
  stop("`groups` muss logical/character/factor/numeric sein.")
}

#' @title Auto-group from observation ratio
#' @description Interne Hilfsfunktion: erzeugt Gruppenindizes anhand des Anteils
#'   beobachteter Werte pro Kurve.
#' @keywords internal
#' @param O_mat Beobachtungsmatrix (n x m) mit 1 für beobachtet, 0 sonst.
#' @param delta Schwellwert in [0,1] für den minimalen beobachteten Anteil pro Kurve,
#'   um zur Gruppe A zu gehören.
#' @return Logischer Vektor (Länge n): TRUE = Gruppe A (ausreichend beobachtet).
.tfu_group_from_delta <- function(O_mat, delta) rowMeans(O_mat != 0) >= delta

#' @title Coerce tfd to matrix + grid
#' @description Interne Hilfsfunktion: konvertiert ein `tfd`-Objekt in Matrix + Gitter.
#'   Liest (falls möglich) das Grid aus den (numerisch interpretierbaren) Spaltennamen.
#'   Falls nicht möglich, wird auf ein äquidistantes Gitter \code{seq(0, 1, length.out = m)} gefallen.
#'   Ist das Grid nicht aufsteigend, werden Grid und Matrixspalten intern sortiert.
#' @keywords internal
#' @param fd Objekt der Klasse `tfd` (oder kompatibel); wird via \code{as.matrix()} in eine
#'   numerische Matrix (n x m) überführt.
#' @return Liste mit Elementen \code{X_obs} (Matrix n x m mit NAs für fehlende Werte)
#'   und \code{grid} (numerischer Vektor der Länge m).
.tfu_from_fd <- function(fd) {
  X_try <- tryCatch(
    { suppressWarnings(as.matrix(fd)) },
    error = function(e) {
      stop("Konnte `fd` nicht via as.matrix() in eine Matrix überführen: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.matrix(X_try)) stop("`as.matrix(fd)` lieferte kein Matrix-Objekt.", call. = FALSE)
  if (ncol(X_try) < 2L) stop("`fd` enthält weniger als 2 Zeitpunkte (ncol < 2).", call. = FALSE)
  
  cn <- colnames(X_try)
  g  <- suppressWarnings(as.numeric(cn))
  if (is.null(cn) || any(is.na(g))) {
    g <- seq(0, 1, length.out = ncol(X_try))
    warning("Grid konnte nicht sicher aus `fd` gelesen werden – Fallback auf seq(0, 1, length.out = m).")
  } else {
    if (anyDuplicated(g)) {
      warning("Grid aus `fd` enthält doppelte Stützstellen – Reihenfolge wird nach aufsteigendem Grid sortiert.")
    }
    if (is.unsorted(g)) {
      o <- order(g)
      g <- g[o]
      X_try <- X_try[, o, drop = FALSE]
      warning("Grid aus `fd` war nicht aufsteigend – Grid und Matrixspalten wurden intern sortiert.")
    }
  }
  storage.mode(X_try) <- "numeric"
  list(X_obs = X_try, grid = g)
}

#' @title Subdomain selector (strict + overlap fallback)
#' @description Wählt eine Subdomäne (Zeitstellen), auf der beide Gruppen ausreichend
#'   Daten haben. Falls das nicht gelingt, wird auf reine Überlappung (beide Gruppen > 0)
#'   zurückgegriffen. Ist auch das nicht erfüllbar, wird mit Fehler abgebrochen.
#' @keywords internal
#' @param O_mat Beobachtungsmatrix (n x m) mit 1/0.
#' @param group_A Logischer Vektor der Länge n: TRUE = Gruppe A, FALSE = B.
#' @param min_frac Minimaler Anteil (z. B. 0.10), der in jeder Gruppe pro Zeitstelle
#'   vorliegen muss, damit die Zeitstelle in die Subdomäne aufgenommen wird.
#' @return Liste mit \code{idx} (Spaltenindizes), \code{min_frac_used} und \code{fallback}.
.tfu_subdomain_idx_paper <- function(O_mat, group_A, min_frac = 0.10) {
  stopifnot(is.matrix(O_mat), length(group_A) == nrow(O_mat))
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  
  cA <- colSums(O_mat * IA)
  cB <- colSums(O_mat * IB)
  
  idx_strict <- which(pmin(cA, cB) > n * min_frac)
  if (length(idx_strict) >= 2L) {
    return(list(idx = idx_strict, min_frac_used = min_frac, fallback = NULL))
  }
  
  idx_overlap <- which(cA > 0 & cB > 0)
  if (length(idx_overlap) >= 2L) {
    warning("Subdomain: fallback overlap (beide Gruppen mit >0 Beobachtungen).")
    return(list(idx = idx_overlap, min_frac_used = NA_real_, fallback = "overlap"))
  }
  
  stop(
    paste0(
      "Keine geeignete Subdomäne gefunden: ",
      "Weder striktes Kriterium (min_frac = ", format(min_frac), ") ",
      "noch Überlappung (beide Gruppen > 0) sind auf mindestens 2 Zeitstellen erfüllbar. ",
      "Passe ggf. Gruppierung oder 'min_frac' an."
    ),
    call. = FALSE
  )
}

#' @title Prepare inputs (fd or matrix)
#' @description Vereinheitlicht Eingaben, baut Beobachtungsmatrix und Gruppen auf.
#'   Stellt außerdem sicher, dass Gruppe A die vollständigere Gruppe ist.
#' @keywords internal
#' @param fd Optionales `tfd`-Objekt.
#' @param X_obs Optional: numerische Matrix (n x m) mit NAs.
#' @param groups Optional: Gruppierungsvektor (siehe `.tfu_groups_to_logical()`).
#' @param observed_ratio Anteil in [0,1], ab dem eine Kurve automatisch Gruppe A ist,
#'   wenn `groups` fehlt (Default 1).
#' @return Liste mit `X_obs`, `O_mat` (1/0), `group_A` (logical) und `grid`.
.tfu_prepare_inputs <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1) {
  if (!is.null(fd)) {
    conv  <- .tfu_from_fd(fd)
    X_obs <- conv$X_obs
    grid_vec <- conv$grid
  } else {
    if (is.null(X_obs)) stop("Entweder `fd` oder `X_obs` muss angegeben werden.")
    grid_vec <- seq(0, 1, length.out = ncol(X_obs))
  }
  
  n <- nrow(X_obs)
  O_mat <- 1L * (!is.na(X_obs))
  
  if (!is.null(groups)) {
    if (length(groups) != n) stop("`groups` muss Länge nrow(X_obs) haben.")
    if (any(is.na(groups))) stop("`groups` darf keine NAs enthalten.")
    group_A <- .tfu_groups_to_logical(groups)
    
    obs_frac <- rowMeans(O_mat != 0)
    meanA <- mean(obs_frac[group_A], na.rm = TRUE)
    meanB <- mean(obs_frac[!group_A], na.rm = TRUE)
    if (is.finite(meanA) && is.finite(meanB) && meanA < meanB) {
      group_A <- !group_A
      message(sprintf(
        "Hinweis: Gruppenlabels getauscht → Gruppe A ist nun die vollständigere Gruppe (Ø A=%.3f, B=%.3f).",
        meanA, meanB
      ))
    }
  } else {
    group_A <- .tfu_group_from_delta(O_mat, observed_ratio)
    if (sum(group_A) == 0L || sum(!group_A) == 0L) {
      obs_frac <- rowMeans(O_mat != 0)
      medf <- median(obs_frac, na.rm = TRUE)
      group_A <- obs_frac >= medf
      warning("Automatische Gruppierung: observed_ratio angepasst.")
    }
  }
  
  list(
    X_obs   = X_obs,
    O_mat   = O_mat,
    group_A = as.logical(group_A),
    grid    = grid_vec
  )
}

#' @title Trapezoidal integration weights
#' @keywords internal
.tfu_trap_weights <- function(x) {
  m <- length(x); if (m == 1L) return(1)
  w <- numeric(m); w[1] <- (x[2]-x[1])/2; w[m] <- (x[m]-x[m-1])/2
  if (m > 2L) w[2:(m-1)] <- (x[3:m] - x[1:(m-2)])/2
  w
}

#' @title Available-mean estimators by group (fail-fast)
#' @param eps Numerischer Stabilitätsparameter (Default 1e-8). Verhindert Divisionen durch ~0.
#' @keywords internal
.tfu_available_means <- function(X_obs, O_mat, group_A, eps = 1e-8) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  
  pA_hat <- colSums(O_mat * IA) / n
  pB_hat <- colSums(O_mat * IB) / n
  
  # Stabilisierung
  pA_eff <- pmax(pA_hat, eps)
  pB_eff <- pmax(pB_hat, eps)
  
  muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_eff)
  muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_eff)
  
  list(muA = muA_hat, muB = muB_hat, pA = pA_eff, pB = pB_eff)
}

#' @title Corrected covariance under partial observation 
#' @param eps Numerischer Stabilitätsparameter (Default 1e-8). Verhindert Divisionen durch ~0.
#' @keywords internal
.tfu_corrected_cov <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat, eps = 1e-8) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  
  Xtilde <- X_obs
  Xtilde[ group_A, ] <- sweep(X_obs[ group_A, , drop = FALSE], 2, muA_hat, `-`)
  Xtilde[!group_A, ] <- sweep(X_obs[!group_A, , drop = FALSE], 2, muB_hat, `-`)
  Xtilde[is.na(Xtilde)] <- 0
  
  # Stabilisierung p(t)
  pA_eff <- pmax(pA_hat, eps)
  pB_eff <- pmax(pB_hat, eps)
  
  A_resid <- sweep((Xtilde * O_mat) * IA, 2, pA_eff, "/")
  B_resid <- sweep((Xtilde * O_mat) * IB, 2, pB_eff, "/")
  
  K_hat <- (crossprod(A_resid) + crossprod(B_resid)) / n
  (K_hat + t(K_hat)) / 2
}

#' @title KL basis from covariance (strict PSD check + optional ridge)
#' @param eps Optionaler Ridge/Jitter (Default 0): fügt `eps * mean(diag(S))` auf der Diagonalen
#'   der gewichteten Kovarianz hinzu, um numerische PSD‑Probleme zu vermeiden.
#' @keywords internal
.tfu_kl_from_cov <- function(K, grid, tol = sqrt(.Machine$double.eps), eps = 0) {
  w  <- .tfu_trap_weights(grid)
  sw <- sqrt(w)
  S <- (sw %o% sw) * K
  S <- (S + t(S)) / 2
  
  # Optionaler numerischer Jitter zur PSD-Stabilisierung
  if (isTRUE(eps > 0)) {
    scale_S <- mean(diag(S))
    if (!is.finite(scale_S) || scale_S <= 0) scale_S <- 1
    S <- S + diag(eps * scale_S, nrow(S))
  }
  
  ev_test <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev_test < -tol * max(1, abs(ev_test[1])))) {
    stop(sprintf(
      "Gewichtete Kovarianzmatrix ist nicht PSD (min. Eigenwert = %.4g).",
      min(ev_test)
    ), call. = FALSE)
  }
  
  ev  <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, 0)  # num. Stabilisierung
  U   <- ev$vectors
  
  phi   <- sweep(U, 1, sw, "/")
  norms <- sqrt(colSums(phi^2 * w))
  phi   <- sweep(phi, 2, norms, "/")
  
  list(lam = lam, phi = phi, w = w)
}

#' @title Build extended htest object
#' @keywords internal
.tfu_make_htest_ext <- function(stat_name, stat_value, p_value, method, data_name,
                                estimate = NULL, conf.int = NULL, parameter = NULL,
                                null.value = 0, alternative = "two.sided") {
  out <- list(
    statistic  = stats::setNames(as.numeric(stat_value), stat_name),
    parameter  = parameter,
    p.value    = as.numeric(p_value),
    conf.int   = conf.int,
    estimate   = estimate,
    null.value = c("difference in mean functions" = null.value),
    alternative= alternative,
    method     = method,
    data.name  = data_name
  )
  class(out) <- "htest"
  out
}

# ----------------------------- Algorithmus 1: L2-Test --------------------------

#' @title Asymptotischer L2‑Test für MCAR (Mittelwertdifferenzen)
#' @description Testet \eqn{H_0:\ \mu_A=\mu_B} mittels \eqn{T_{µ,L2}=n\|\hat\mu_A-\hat\mu_B\|^2_{L2}}.
#'   Grenzverteilung via KL‑Mischung; **optionale Parallelisierung** der Monte‑Carlo‑Approximation.
#' @param parallel Logisch: Monte‑Carlo‑Mischung parallel berechnen?
#' @param ncpus Anzahl Kerne (nur bei `parallel = TRUE`).
#' @param eps Numerischer Stabilitätsparameter (Default 1e-8). Wird in allen internen
#'   Divisionen sowie in der KL‑Zerlegung (kleiner Ridge) genutzt.
#' @inheritParams .tfu_prepare_inputs
#' @inheritParams .tfu_subdomain_idx_paper
#' @param fve Fraction of variance explained (0–1) zur Wahl der Anzahl KL‑Komponenten q.
#' @param B Anzahl Simulationen für die KL‑Mischung (Monte‑Carlo zur p‑Wert‑Approximation).
#' @param seed Optionaler RNG‑Seed.
#' @param alpha Signifikanzniveau für das (optionale) Intervall in `conf.int`.
#' @return Erweitertes `htest`‑Objekt.
#' @examples
#' # res <- asym_mean_L2_test(X_obs = X, groups = grp, parallel = TRUE)
#' @import foreach doRNG doParallel
asym_mean_L2_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, B = 5000,
                              min_frac = 0.10, seed = NULL, alpha = 0.05,
                              parallel = FALSE, eps = 1e-8,
                              ncpus = parallel::detectCores(logical = FALSE)) {
  
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB, eps = eps)
  KL <- .tfu_kl_from_cov(K, g, eps = eps) ; lam <- KL$lam
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]
  
  # Monte Carlo mixture: parallel oder vektorisiert
  stopClusterOnExit <- FALSE
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(max(1L, ncpus))
    stopClusterOnExit <- TRUE
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed)
    on.exit(if (stopClusterOnExit) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    W <- foreach::foreach(b = 1:B, .combine = c, .inorder = FALSE) %dorng% {
      z <- rnorm(q)
      sum((z^2) * lam_q)
    }
  } else {
    Z <- matrix(rnorm(q * B), nrow = q)
    W <- colSums((Z^2) * lam_q)
  }
  
  p <- (sum(W >= T_L2) + 1) / (length(W) + 1)
  
  ci <- c(mean(diff) - sd(W), mean(diff) + sd(W)); attr(ci,"conf.level") <- 1 - alpha
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  
  .tfu_make_htest_ext("T_{mu,L2}", T_L2, p,
                      if (isTRUE(parallel)) "L2-Test (KL-Mischung; parallel MC)"
                      else "L2-Test (KL-Mischung)",
                      data_name,
                      estimate = c("mean(A)" = mean(muA), "mean(B)" = mean(muB),
                                   "diff(A-B)" = mean(diff)),
                      conf.int = ci,
                      parameter = c(q = q, m = length(g)))
}

# ----------------------------- Algorithmus 2: Supremum-Test --------------------

#' @title Asymptotischer Supremums‑Test für MCAR (Mittelwertdifferenzen)
#' @description Testet \eqn{H_0:\ \mu_A=\mu_B} mittels \eqn{T_{µ,D}=\sqrt{n}\|\hat\mu_A-\hat\mu_B\|_\infty}.
#'   GP‑Approximation über KL‑Basis; **optionale Parallelisierung** der MC‑Supremums‑Simulation.
#' @inheritParams asym_mean_L2_test
#' @param fve Fraction of variance explained (0–1) für die GP‑Approximation (typisch 0.95).
#' @return Erweitertes `htest`‑Objekt.
#' @examples
#' # res <- asym_mean_sup_test(X_obs = X, groups = grp, parallel = TRUE)
#' @import foreach doRNG doParallel
asym_mean_sup_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                               fve = 0.95, B = 5000,
                               min_frac = 0.10, seed = NULL, alpha = 0.05,
                               parallel = FALSE, eps = 1e-8,
                               ncpus = parallel::detectCores(logical = FALSE)) {
  
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  T_D <- sqrt(n) * max(abs(diff))
  
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB, eps = eps)
  KL <- .tfu_kl_from_cov(K, g, eps = eps); lam <- KL$lam; phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]; phi_q <- phi[, seq_len(q), drop = FALSE]
  A <- sweep(phi_q, 2, sqrt(lam_q), "*")  # m x q
  
  stopClusterOnExit <- FALSE
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(max(1L, ncpus))
    stopClusterOnExit <- TRUE
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed)
    on.exit(if (stopClusterOnExit) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    W <- foreach::foreach(b = 1:B, .combine = c, .inorder = FALSE) %dorng% {
      z <- rnorm(q)
      max(abs(A %*% z))
    }
  } else {
    Z <- matrix(rnorm(q * B), nrow = q)
    gp_vals <- A %*% Z  # m x B
    W <- apply(abs(gp_vals), 2, max)
  }
  
  p <- (sum(W >= T_D) + 1) / (length(W) + 1)
  
  ci <- c(mean(diff) - sd(W), mean(diff) + sd(W)); attr(ci,"conf.level") <- 1 - alpha
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  
  .tfu_make_htest_ext("T_{mu,D}", T_D, p,
                      if (isTRUE(parallel)) "Supremum-Test (GP-Approx; parallel MC)"
                      else "Supremum-Test (GP-Approx)",
                      data_name,
                      estimate = c("mean(A)" = mean(muA), "mean(B)" = mean(muB),
                                   "diff(A-B)" = mean(diff)),
                      conf.int = ci,
                      parameter = c(q = q, m = length(g)))
}

# ----------------------------- Algorithmus 3: Konfidenzbänder ------------------

#' @title Simultane Konfidenzbänder (asymptotisch) für \eqn{\mu_A-\mu_B}
#' @description Wie beim Supremums‑Test; **optionale Parallelisierung** der MC‑Supremums‑Simulation.
#' @inheritParams asym_mean_sup_test
#' @return Erweitertes `htest`‑Objekt (mit statistic = q_alpha).
#' @examples
#' # cb <- asym_conf_bands(X_obs = X, groups = grp, parallel = TRUE)
#' @import foreach doRNG doParallel
asym_conf_bands <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                            alpha = 0.05, fve = 0.95, B = 5000,
                            min_frac = 0.10, seed = NULL,
                            parallel = FALSE, eps = 1e-8,
                            ncpus = parallel::detectCores(logical = FALSE)) {
  
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB, eps = eps)
  KL <- .tfu_kl_from_cov(K, g, eps = eps); lam <- KL$lam; phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]; phi_q <- phi[, seq_len(q), drop = FALSE]
  A <- sweep(phi_q, 2, sqrt(lam_q), "*")
  
  stopClusterOnExit <- FALSE
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(max(1L, ncpus))
    stopClusterOnExit <- TRUE
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed)
    on.exit(if (stopClusterOnExit) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    W <- foreach::foreach(b = 1:B, .combine = c, .inorder = FALSE) %dorng% {
      z <- rnorm(q)
      max(abs(A %*% z))
    }
  } else {
    Z <- matrix(rnorm(q * B), nrow = q)
    gp_vals <- A %*% Z
    W <- apply(abs(gp_vals), 2, max)
  }
  
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha))
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth; upper <- diff + halfwidth
  
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB), "diff(A-B)" = mean(diff))
  ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci,"conf.level") <- 1 - alpha
  
  .tfu_make_htest_ext("Band-width", q_alpha, NA,
                      if (isTRUE(parallel)) "Simultane Konfidenzbänder (asympt.; parallel MC)"
                      else "Simultane Konfidenzbänder (asymptotisch)",
                      if (!is.null(fd)) "fd" else "X_obs",
                      estimate = est_vec, conf.int = ci, parameter = c(m = length(g)))
}

# -----------------------------------------------------------------------------
# Bootstrap-Tests (foreach + doRNG) & Bootstrap-Bänder — eps‑stabilisiert
# -----------------------------------------------------------------------------

#' @title Gruppenspezifische Bootstrap‑Mittelwerte (intern)
#' @param eps Numerischer Stabilitätsparameter (Default 1e-8).
#' @keywords internal
if (!exists(".tfu_boot_group_mean", mode = "function")) {
  .tfu_boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
    p_hat <- colSums(O_grp) / n_total
    p_hat <- pmax(p_hat, eps)  # Stabilisierung
    numer <- colSums(replace(X_grp, is.na(X_grp), 0))
    numer / (n_total * p_hat)
  }
}

#' Bootstrap-basierter Mean‑Test (L2 oder Supremum) — foreach + doRNG
#'
#' @description Liefert Bootstrap‑p‑Werte sowie (optional) simultane Bänder.
#' @param fd Optional: `tfd`‑Objekt. Wird gegenüber `X_obs` priorisiert.
#' @param X_obs Optional: Matrix (n x m) mit NAs.
#' @param groups Optional: Gruppierungsvektor (zwei Ausprägungen). Fehlt er, automatische
#'   Gruppierung über `observed_ratio`.
#' @param observed_ratio Anteil in [0,1] für Auto‑Gruppierung (siehe oben).
#' @param B Anzahl Bootstrap‑Zyklen.
#' @param min_frac Minimaler Gruppenbeobachtungsanteil pro Zeitstelle.
#' @param alpha Signifikanzniveau für simultane Bänder (wenn `compute_bands = TRUE`).
#' @param parallel Logisch: Parallelausführung via `foreach`?
#' @param ncpus Anzahl Kerne für Cluster (nur bei `parallel = TRUE` relevant).
#' @param seed Optionaler RNG‑Seed (reproduzierbar via doRNG).
#' @param stat Kennung der Teststatistik: `"L2"` oder `"D"` (Supremum).
#' @param compute_bands Logisch: zusätzlich simultane Bänder berechnen?
#' @param return_boot Logisch: Bootstrap‑Statistiken dem Objekt beilegen?
#' @param eps Numerischer Stabilitätsparameter (Default 1e-8).
#' @return Erweitertes `htest`‑Objekt (inkl. Zusatzfelder).
#' @examples
#' # res <- boot_mean_test(X_obs = X, groups = grp, parallel = TRUE)
#' @import foreach doRNG doParallel
boot_mean_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                           B = 5000,
                           min_frac = 0.10, alpha = 0.05,
                           parallel = TRUE,
                           ncpus = parallel::detectCores(logical = FALSE),
                           seed = NULL,
                           stat = c("L2", "D"),
                           compute_bands = TRUE,
                           return_boot = FALSE,
                           eps = 1e-8) {
  stat <- match.arg(stat)
  prep <- .tfu_prepare_inputs(fd = fd, X_obs = X_obs, groups = groups, observed_ratio = observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  est  <- .tfu_available_means(X, O, group_A, eps = eps)
  muA  <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  w <- .tfu_trap_weights(g)
  T_L2 <- n * sum((diff^2) * w)
  T_D  <- sqrt(n) * max(abs(diff))
  
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  if (!is.null(seed)) set.seed(seed)
  doRNG::registerDoRNG(if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed)
  
  stopClusterOnExit <- FALSE
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(max(1L, ncpus))
    stopClusterOnExit <- TRUE
    doParallel::registerDoParallel(cl)
    on.exit(if (stopClusterOnExit) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }
  
  boot_mat <- foreach::foreach(b = 1:B, .combine = rbind, .inorder = FALSE) %dorng% {
    samp_idx <- sample.int(n, n, replace = TRUE)
    sampA <- samp_idx[group_A[samp_idx] == 1L]
    sampB <- samp_idx[group_A[samp_idx] == 0L]
    
    # Wenn eine Gruppe im Resample komplett leer ist, Abbruch dieses Zugs
    if (length(sampA) == 0L || length(sampB) == 0L) return(c(L2 = NA_real_, D = NA_real_))
    
    if (length(sampA)) {
      OA <- O[sampA,, drop = FALSE]
      XA <- X_cent[sampA,, drop = FALSE]
      pA <- colSums(OA) / n
      # Stabilisierung: verhindert hartes Abbrechen, aber vermeidet Spalten ohne beobachtete Werte in beiden Gruppen
      pA <- pmax(pA, eps)
      numA <- colSums(replace(XA, is.na(XA), 0))
      muA_b <- numA / (n * pA)
    } else muA_b <- rep(0, ncol(X_cent))
    
    if (length(sampB)) {
      OB <- O[sampB,, drop = FALSE]
      XB <- X_cent[sampB,, drop = FALSE]
      pB <- colSums(OB) / n
      pB <- pmax(pB, eps)
      numB <- colSums(replace(XB, is.na(XB), 0))
      muB_b <- numB / (n * pB)
    } else muB_b <- rep(0, ncol(X_cent))
    
    d_b <- muA_b - muB_b
    c(L2 = n * sum((d_b^2) * w),
      D  = sqrt(n) * max(abs(d_b)))
  }
  
  boot_L2 <- boot_mat[, "L2"]
  boot_D  <- boot_mat[, "D"]
  
  n_bad <- sum(is.na(boot_L2) | is.na(boot_D))
  if (n_bad > 0) {
    warning(sprintf("Bootstrap: %d Iterationen übersprungen (%.1f%%).", 
                    n_bad, 100 * n_bad / B))
  }
  boot_L2 <- boot_L2[!is.na(boot_L2)]
  boot_D  <- boot_D[!is.na(boot_D)]
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (length(boot_L2) + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (length(boot_D)  + 1)
  
  if (isTRUE(compute_bands)) {
    q_alpha   <- as.numeric(stats::quantile(boot_D, probs = 1 - alpha, names = FALSE))
    halfwidth <- q_alpha / sqrt(n)
    ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci, "conf.level") <- 1 - alpha
    lower <- diff - halfwidth
    upper <- diff + halfwidth
  } else {
    ci <- NULL; q_alpha <- NA_real_; lower <- upper <- NULL
  }
  
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB), "diff(A-B)" = mean(diff))
  data_name <- if (!is.null(fd)) "fd (tfd/tfd_irreg via tf_gather)" else "X_obs"
  
  if (identical(stat, "L2")) {
    hobj <- .tfu_make_htest_ext("T_{mu,L2}", T_L2, p_L2,
                                "Bootstrap Mean-Test (L2; foreach+doRNG)", data_name,
                                estimate = est_vec, conf.int = ci, parameter = c(m = length(g), B = B - n_bad),
                                null.value = 0, alternative = "two.sided")
  } else {
    hobj <- .tfu_make_htest_ext("T_{mu,D}", T_D, p_D,
                                "Bootstrap Mean-Test (Supremum; foreach+doRNG)", data_name,
                                estimate = est_vec, conf.int = ci, parameter = c(m = length(g), B = B - n_bad),
                                null.value = 0, alternative = "two.sided")
  }
  
  hobj$grid <- g; hobj$muA <- muA; hobj$muB <- muB; hobj$diff <- diff
  hobj$idx  <- idx; hobj$min_frac_used <- sub$min_frac_used; hobj$fallback <- sub$fallback
  hobj$alpha <- alpha; hobj$q_alpha <- q_alpha
  hobj$lower <- lower; hobj$upper <- upper
  if (isTRUE(return_boot)) { hobj$boot_L2 <- boot_L2; hobj$boot_D <- boot_D }
  
  hobj
}

# -----------------------------------------------------------------------------
# Nur Bootstrap-Konfidenzbänder (Wrapper, Algo 7)
# -----------------------------------------------------------------------------

#' @title Bootstrap‑Konfidenzbänder (Wrapper)
#' @description Bequemer Wrapper für `boot_mean_test(stat="D", compute_bands=TRUE)`.
#' @inheritParams boot_mean_test
#' @param with_tests Logisch: komplettes `htest`-Objekt zurückgeben (`TRUE`) oder nur Kernfelder?
#' @return Liste oder `htest`.
#' @examples
#' # bands <- boot_conf_bands(X_obs = X, groups = grp, B = 2000, alpha = 0.05)
boot_conf_bands <- function(fd = NULL, X_obs = NULL, groups = NULL,
                            observed_ratio = 1,
                            B = 10000, min_frac = 0.10,
                            seed = NULL, alpha = 0.05,
                            return_boot = FALSE, with_tests = FALSE,
                            parallel = TRUE, eps = 1e-8,
                            ncpus = parallel::detectCores(logical = FALSE)) {
  res <- boot_mean_test(
    fd = fd, X_obs = X_obs, groups = groups, observed_ratio = observed_ratio,
    B = B, min_frac = min_frac, seed = seed, alpha = alpha,
    parallel = parallel, ncpus = ncpus, eps = eps,
    compute_bands = TRUE, return_boot = return_boot, stat = "D"
  )
  
  if (!isTRUE(with_tests)) {
    keep <- c("grid","muA","muB","diff","lower","upper","q_alpha",
              "idx","min_frac_used","fallback","alpha","estimate","conf.int")
    return(res[keep])
  }
  res
}