# =============================================================================
#  Tidyfun MCAR: Mean-based Tests & Bands (full utilities; fd-or-matrix)
# -----------------------------------------------------------------------------
#  Implements (organized as requested):
#    - (Algo 1)   asym_mean_L2_test():            L2 test (asymptotic, KL-mixture)
#    - (Algo 2+3) asym_sup_and_bands():           Supremum test + asymptotic bands
#    - (Algo 5+7) boot_mean_test_and_bands():     Bootstrap tests + bootstrap bands
# -----------------------------------------------------------------------------
#  All return values are extended `htest` objects (t.test-style) where applicable.
#  The bootstrap combo returns two htests (L2 & D) plus band info.
# =============================================================================

#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------------------------- Utilities (full) ---------------------------------

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

.tfu_group_from_delta <- function(O_mat, delta) rowMeans(O_mat != 0) >= delta

# fd (tfd) -> (X_obs, grid)
.tfu_from_fd <- function(fd) {
  X_try <- try(suppressWarnings(as.matrix(fd)), silent = TRUE)
  if (inherits(X_try, "try-error"))
    stop("Konnte `fd` nicht via as.matrix() in eine Matrix überführen.")
  cn <- colnames(X_try)
  g  <- suppressWarnings(as.numeric(cn))
  if (is.null(cn) || any(is.na(g))) {
    g <- seq(0, 1, length.out = ncol(X_try))
    warning("Grid nicht sicher aus `fd` gelesen -> seq(0,1,len=m).")
  }
  storage.mode(X_try) <- "numeric"
  list(X_obs = X_try, grid = g)
}

# Subdomain (10%) mit Fallback
.tfu_subdomain_idx_paper <- function(O_mat, group_A, min_frac = 0.10, m_min = 2L) {
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  cA <- colSums(O_mat * IA); cB <- colSums(O_mat * IB)
  
  idx <- which(cA >= min_frac * n & cB >= min_frac * n)
  if (length(idx) >= max(2L, m_min)) return(list(idx = idx, min_frac_used = min_frac, fallback = NULL))
  
  idx2 <- which(cA > 0 & cB > 0)
  if (length(idx2) >= max(2L, m_min)) {
    warning("Subdomain (10%): fallback overlap.")
    return(list(idx = idx2, min_frac_used = NA_real_, fallback = "overlap"))
  }
  
  idx3 <- which(colSums(O_mat) > 0)
  if (length(idx3) >= max(2L, m_min)) {
    warning("Subdomain: fallback any-observed.")
    return(list(idx = idx3, min_frac_used = NA_real_, fallback = "any-observed"))
  }
  
  warning("Subdomain: trivial fallback.")
  list(idx = seq_len(min(2L, ncol(O_mat))), min_frac_used = NA_real_, fallback = "trivial")
}

# Inputs vorbereiten
.tfu_prepare_inputs <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1) {
  if (!is.null(fd)) {
    conv  <- .tfu_from_fd(fd)
    X_obs <- conv$X_obs
    grid  <- conv$grid
  } else {
    stopifnot(!is.null(X_obs))
    if (!is.matrix(X_obs)) X_obs <- as.matrix(X_obs)
    storage.mode(X_obs) <- "numeric"
    grid <- seq(0, 1, length.out = ncol(X_obs))
  }
  n <- nrow(X_obs); m <- ncol(X_obs)
  if (m < 2L) stop("Daten enthalten < 2 Zeitpunkte.")
  O_mat <- 1L * (!is.na(X_obs))
  
  if (!is.null(groups)) {
    if (length(groups) != n) stop("`groups` muss Länge nrow(X_obs) haben.")
    if (any(is.na(groups))) stop("`groups` darf keine NAs enthalten.")
    group_A <- .tfu_groups_to_logical(groups)
  } else {
    group_A <- .tfu_group_from_delta(O_mat, observed_ratio)
    if (sum(group_A) == 0L || sum(!group_A) == 0L) {
      obs_frac <- rowMeans(O_mat != 0)
      medf <- median(obs_frac)
      group_A <- obs_frac >= medf
      warning("Automatische Gruppierung: observed_ratio angepasst.")
    }
  }
  list(X_obs = X_obs, O_mat = O_mat, group_A = as.logical(group_A), grid = grid)
}

# Trapezgewichte
.tfu_trap_weights <- function(x) {
  m <- length(x); if (m == 1L) return(1)
  w <- numeric(m); w[1] <- (x[2]-x[1])/2; w[m] <- (x[m]-x[m-1])/2
  if (m > 2L) w[2:(m-1)] <- (x[3:m] - x[1:(m-2)])/2
  w
}

# Mittelwertschätzer (verfügbar)
.tfu_available_means <- function(X_obs, O_mat, group_A, eps = 1e-8) {
  n  <- nrow(X_obs); IA <- as.numeric(group_A); IB <- 1 - IA
  pA_hat <- pmax(colSums(O_mat * IA) / n, eps)
  pB_hat <- pmax(colSums(O_mat * IB) / n, eps)
  muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_hat)
  list(muA = muA_hat, muB = muB_hat, pA = pA_hat, pB = pB_hat)
}

# Kovarianz (korrigiert)
.tfu_corrected_cov <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X_obs); IA <- as.numeric(group_A); IB <- 1 - IA
  Xc <- X_obs
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  if (length(A_idx)) Xc[A_idx, ] <- sweep(X_obs[A_idx, , drop = FALSE], 2, muA_hat, `-`)
  if (length(B_idx)) Xc[B_idx, ] <- sweep(X_obs[B_idx, , drop = FALSE], 2, muB_hat, `-`)
  Xc[O_mat == 0] <- 0; Xc[is.na(Xc)] <- 0
  XA <- Xc[A_idx, , drop = FALSE]; XB <- Xc[B_idx, , drop = FALSE]
  if (nrow(XA)) XA <- sweep(XA, 2, pA_hat, "/")
  if (nrow(XB)) XB <- sweep(XB, 2, pB_hat, "/")
  Sum_A <- if (nrow(XA)) crossprod(XA) else matrix(0, ncol(X_obs), ncol(X_obs))
  Sum_B <- if (nrow(XB)) crossprod(XB) else matrix(0, ncol(X_obs), ncol(X_obs))
  K <- (Sum_A + Sum_B) / n
  (K + t(K)) / 2
}

# KL-Basis aus Kovarianz
.tfu_kl_from_cov <- function(K, grid) {
  w <- .tfu_trap_weights(grid); sw <- sqrt(w)
  S <- (sw * t(sw * K)); S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, 0); U <- ev$vectors
  phi <- sweep(U, 1, sw, "/")
  norms <- sqrt(colSums(phi^2 * w)); phi <- sweep(phi, 2, norms, "/")
  list(lam = lam, phi = phi, w = w)
}

# erweitertes htest erstellen
.tfu_make_htest_ext <- function(stat_name, stat_value, p_value, method, data_name,
                                estimate = NULL, conf.int = NULL, parameter = NULL,
                                null.value = 0, alternative = "two.sided") {
  out <- list(
    statistic  = stats::setNames(as.numeric(stat_value), stat_name),
    parameter  = parameter,
    p.value    = if (!is.null(p_value)) as.numeric(p_value) else NA_real_,
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

# =============================================================================
# (Algo 1)  Asymptotischer L2-Test
# =============================================================================
asym_mean_L2_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, B = 5000, eps = 1e-8,
                              min_frac = 0.10, seed = NULL, alpha = 0.05) {
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  # exakte L2-Integration via tf
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  # KL-Mischung
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .tfu_kl_from_cov(K, g); lam <- KL$lam
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; lam_q <- lam[seq_len(q)]
  Z <- matrix(rnorm(length(lam_q) * B), nrow = length(lam_q))
  W <- colSums((Z^2) * lam_q)
  p <- (sum(W >= T_L2) + 1) / (B + 1)
  
  ci <- c(mean(diff) - sd(W), mean(diff) + sd(W)); attr(ci,"conf.level") <- 1-alpha
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  .tfu_make_htest_ext("T_{mu,L2}", T_L2, p,
                      "L2-Test (KL-Mischung, asymptotisch)", data_name,
                      estimate = c("mean(A)" = mean(muA), "mean(B)" = mean(muB),
                                   "diff(A-B)" = mean(diff)),
                      conf.int = ci,
                      parameter = c(q = q, m = length(g)))
}

# =============================================================================
# (Algo 2 + 3)  Supremum-Test (asymptotisch) + asymptotische simultane Bänder
# =============================================================================
asym_sup_and_bands <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                               alpha = 0.05, fve = 0.95, B = 5000, eps = 1e-8,
                               min_frac = 0.10, seed = NULL) {
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A, eps = eps)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  T_D <- sqrt(n) * max(abs(diff))
  
  # KL-Basis und GP-Approx
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .tfu_kl_from_cov(K, g); lam <- KL$lam; phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]
  lam_q <- lam[seq_len(q)]; phi_q <- phi[, seq_len(q), drop = FALSE]
  Z <- matrix(rnorm(q * B), nrow = q)
  gp_vals <- sweep(phi_q, 2, sqrt(lam_q), "*") %*% Z
  W <- apply(abs(gp_vals), 2, max)               # Supremum der GP-Realisation
  p <- (sum(W >= T_D) + 1) / (B + 1)
  
  # (Algo 3) Bänder aus der gleichen GP-Approx
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha))
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth; upper <- diff + halfwidth
  
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB), "diff(A-B)" = mean(diff))
  ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci,"conf.level") <- 1 - alpha
  
  hobj <- .tfu_make_htest_ext("T_{mu,D}", T_D, p,
                              "Supremum-Test + simultane Konfidenzbänder (asymptotisch, GP-Approx)",
                              if (!is.null(fd)) "fd" else "X_obs",
                              estimate = est_vec,
                              conf.int = ci,
                              parameter = c(q = q, m = length(g)))
  
  # Zusatzfelder (Bands etc.)
  hobj$grid <- g
  hobj$muA  <- muA
  hobj$muB  <- muB
  hobj$diff <- diff
  hobj$lower <- lower
  hobj$upper <- upper
  hobj$alpha <- alpha
  hobj$q_alpha <- q_alpha
  hobj$idx  <- idx
  hobj$min_frac_used <- sub$min_frac_used
  hobj$fallback      <- sub$fallback
  hobj
}

# =============================================================================
# (Algo 5 + 7)  Bootstrap-Tests (L2 & D) + Bootstrap-Bänder
# =============================================================================
boot_mean_test_and_bands <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                                     B = 5000, eps = 1e-8,
                                     min_frac = 0.10, alpha = 0.05,
                                     parallel = TRUE,
                                     ncpus = parallel::detectCores(logical = FALSE),
                                     seed = NULL,
                                     return_boot = FALSE) {
  # --- Inputs vorbereiten ---
  prep <- .tfu_prepare_inputs(fd = fd, X_obs = X_obs, groups = groups, observed_ratio = observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  # --- verfügbare Mittelwerte + beobachtete Statistiken ---
  est  <- .tfu_available_means(X, O, group_A, eps = eps)
  muA  <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  w <- .tfu_trap_weights(g)
  T_L2 <- n * sum((diff^2) * w)
  T_D  <- sqrt(n) * max(abs(diff))
  
  # --- Residuen zentrieren ---
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  # --- foreach + doRNG Setup ---
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
  
  # --- Bootstrap via foreach (Algo 5) ---
  boot_mat <- foreach::foreach(b = 1:B, .combine = rbind, .inorder = FALSE) %dorng% {
    samp_idx <- sample.int(n, n, replace = TRUE)
    sampA <- samp_idx[group_A[samp_idx] == 1L]
    sampB <- samp_idx[group_A[samp_idx] == 0L]
    
    if (length(sampA)) {
      OA <- O[sampA,, drop = FALSE]
      XA <- X_cent[sampA,, drop = FALSE]
      pA <- pmax(colSums(OA) / n, eps)
      numA <- colSums(replace(XA, is.na(XA), 0))
      muA_b <- numA / (n * pA)
    } else muA_b <- rep(0, ncol(X_cent))
    
    if (length(sampB)) {
      OB <- O[sampB,, drop = FALSE]
      XB <- X_cent[sampB,, drop = FALSE]
      pB <- pmax(colSums(OB) / n, eps)
      numB <- colSums(replace(XB, is.na(XB), 0))
      muB_b <- numB / (n * pB)
    } else muB_b <- rep(0, ncol(X_cent))
    
    d_b <- muA_b - muB_b
    c(L2 = n * sum((d_b^2) * w),
      D  = sqrt(n) * max(abs(d_b)))
  }
  
  boot_L2 <- boot_mat[, "L2"]
  boot_D  <- boot_mat[, "D"]
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (B + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (B + 1)
  
  # --- (Algo 7) Konfidenzbänder (Supremum-basiert) ---
  q_alpha   <- as.numeric(stats::quantile(boot_D, probs = 1 - alpha, names = FALSE))
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth; upper <- diff + halfwidth
  
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB),
               "diff(A-B)" = mean(diff))
  ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci, "conf.level") <- 1 - alpha
  data_name <- if (!is.null(fd)) "fd (tfd/tfd_irreg via tf_gather)" else "X_obs"
  
  # Zwei htest-Objekte erzeugen (L2 & D)
  h_L2 <- .tfu_make_htest_ext(
    "T_{mu,L2}", T_L2, p_L2,
    "Bootstrap Mean-Test (L2; foreach+doRNG)", data_name,
    estimate   = est_vec,
    conf.int   = ci,  # CI bezieht sich auf diff-Mean-Bandbreite (gleich für beide Tests)
    parameter  = c(m = length(g), B = B),
    null.value = 0,
    alternative= "two.sided"
  )
  
  h_D <- .tfu_make_htest_ext(
    "T_{mu,D}", T_D, p_D,
    "Bootstrap Mean-Test (Supremum; foreach+doRNG) + simultane Bänder",
    data_name,
    estimate   = est_vec,
    conf.int   = ci,
    parameter  = c(m = length(g), B = B),
    null.value = 0,
    alternative= "two.sided"
  )
  
  # Gemeinsame Zusatzfelder (Bands etc.)
  for (h in list(h_L2, h_D)) {
    h$grid <- g
    h$muA  <- muA
    h$muB  <- muB
    h$diff <- diff
    h$lower <- lower
    h$upper <- upper
    h$alpha <- alpha
    h$q_alpha <- q_alpha
    h$idx  <- idx
    h$min_frac_used <- sub$min_frac_used
    h$fallback      <- sub$fallback
  }
  
  if (isTRUE(return_boot)) {
    h_L2$boot_L2 <- boot_L2
    h_L2$boot_D  <- boot_D
    h_D$boot_L2  <- boot_L2
    h_D$boot_D   <- boot_D
  }
  
  # Rückgabe: zwei htests + Bands (identisch) – klar getrennt
  list(
    test_L2 = h_L2,
    test_D  = h_D
  )
}

