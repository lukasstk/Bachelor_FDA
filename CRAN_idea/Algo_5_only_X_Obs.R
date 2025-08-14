#' Tidyfun extension: Algorithm 5 (bootstrap for L2 and Supremum) — X_obs-only
#'
#' Öffentliche Signatur erwartet NUR `X_obs` (numeric Matrix n x m mit NAs).
#' Intern werden O_mat, group_A (delta_A=1 mit automatischer Lockerung) und grid erzeugt.
#' Subdomain folgt der 10%-Paper-Regel mit einmaligen Fallbacks + Warnings.
#'
#' Tuning-Parameter (bei Bedarf hier im Funktionskopf anpassbar):
#'   B        = 5000   (# Bootstrap-Ziehungen)
#'   eps      = 1e-8   (Untergrenze für p̂)
#'   min_frac = 0.10   (Paper-Subdomain)
#'   delta_A  = 1      (Startschwelle für A/B-Split; 1 = vollständig vs. unvollständig)
#'   return_boot = FALSE (bei TRUE werden die Bootstrap-Verteilungen mit zurückgegeben)
#'
#' Erweiterung: *Bootstrap-Konfidenzbänder* (Algorithmus 7)
#'   compute_bands = TRUE  (liefert simultane Bänder für mu_A - mu_B)
#'   alpha         = 0.05   (1 - Konfidenzniveau)
#'
#' @importFrom tf tfd tf_integrate
tfu_algo5_bootstrap <- function(X_obs) {
  # ----------------- feste Defaults (hier anpassbar) -----------------
  B <- 5000L
  eps <- 1e-8
  min_frac <- 0.10
  delta_A <- 1
  return_boot <- FALSE
  compute_bands <- TRUE
  alpha <- 0.05
  # -------------------------------------------------------------------
  
  # ---------- Utilities nur definieren, falls noch nicht vorhanden ----
  if (!exists(".tfu_group_from_delta", mode = "function")) {
    .tfu_group_from_delta <- function(O_mat, delta) {
      rowMeans(O_mat != 0) >= delta
    }
  }
  
  if (!exists(".tfu_subdomain_idx_paper", mode = "function")) {
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
  }
  
  if (!exists(".tfu_prepare_inputs", mode = "function")) {
    # nur Minimalprüfung von X_obs; alles Weitere intern bauen
    .tfu_prepare_inputs <- function(X_obs) {
      if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
      storage.mode(X_obs) <- "numeric"
      n <- nrow(X_obs); m <- ncol(X_obs)
      if (m < 2L) stop("X_obs muss >= 2 Spalten (Zeitpunkte) haben.")
      
      O_mat <- 1L * (!is.na(X_obs))
      
      # automatische Gruppierung: start delta_A=1, falls Gruppe leer -> lockern, sonst Median-Split
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
  }
  
  if (!exists(".tfu_available_means", mode = "function")) {
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
  }
  
  if (!exists(".tfu_boot_group_mean", mode = "function")) {
    .tfu_boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
      p_hat <- colSums(O_grp) / n_total
      p_hat <- pmax(p_hat, eps)
      numer <- colSums(replace(X_grp, is.na(X_grp), 0))
      numer / (n_total * p_hat)
    }
  }
  # -------------------------------------------------------------------
  
  # Inputs vorbereiten
  prep <- .tfu_prepare_inputs(X_obs)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  # Subdomain wählen
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  # Available-data Mittelwerte
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  # Beobachtete Statistiken
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2 <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  T_D  <- sqrt(n) * max(abs(diff))
  
  # Residuen pro Gruppe zentrieren (NAs bleiben; O korrigiert im Mean-Schätzer)
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  # Bootstrap
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
  
  # --- Bootstrap-Konfidenzbänder (Algorithmus 7) ---
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
