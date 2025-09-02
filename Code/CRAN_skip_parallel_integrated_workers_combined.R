#' Tidyfun extensions: Mean-based MCAR tests & bands (auto inputs; fd-or-matrix)
#'
#' Implements:
#'  - asym_mean_L2_test(): L2 test (asymptotic, KL-mixture)
#'  - asym_mean_sup_test(): Supremum test (asymptotic, KL-based GP approx) with optional bands
#'  - boot_mean_test(): Bootstrap-based p-values for L2/Supremum (+ optional bands)
#'
#' Legacy (kept as thin wrappers, now deprecated):
#'  - asym_conf_bands()  -> calls asym_mean_sup_test(compute_bands=TRUE, bands_only=TRUE)
#'  - boot_conf_bands()  -> calls boot_mean_test(stat="D", compute_bands=TRUE, bands_only=!with_tests)
#'
#' Erweiterung:
#'  - Rückgabeobjekte sind erweiterte `htest`-Listen mit estimate, null.value,
#'    conf.int, alternative – im Stil von t.test().
#'
#' @keywords methods
#' @importFrom tf tfd tf_integrate
NULL

# --------------------------- Utilities -----------------------------------------

#' @title Convert group labels to logical
#' @keywords internal
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
#' @keywords internal
.tfu_group_from_delta <- function(O_mat, delta) rowMeans(O_mat != 0) >= delta

#' @title Coerce tfd to matrix + grid
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
.tfu_available_means <- function(X_obs, O_mat, group_A) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  
  pA_hat <- colSums(O_mat * IA) / n
  pB_hat <- colSums(O_mat * IB) / n
  
  muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_hat)
  
  list(muA = muA_hat, muB = muB_hat, pA = pA_hat, pB = pB_hat)
}

#' @title Corrected covariance under partial observation
#' @keywords internal
.tfu_corrected_cov <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  
  Xtilde <- X_obs
  Xtilde[ group_A, ] <- sweep(X_obs[ group_A, , drop = FALSE], 2, muA_hat, `-`)
  Xtilde[!group_A, ] <- sweep(X_obs[!group_A, , drop = FALSE], 2, muB_hat, `-`)
  Xtilde[is.na(Xtilde)] <- 0
  
  A_resid <- sweep((Xtilde * O_mat) * IA, 2, pA_hat, "/")
  B_resid <- sweep((Xtilde * O_mat) * IB, 2, pB_hat, "/")
  
  K_hat <- (crossprod(A_resid) + crossprod(B_resid)) / n
  (K_hat + t(K_hat)) / 2
}

#' @title KL basis from covariance (strict PSD check)
#' @keywords internal
.tfu_kl_from_cov <- function(K, grid, tol = sqrt(.Machine$double.eps)) {
  w  <- .tfu_trap_weights(grid)
  sw <- sqrt(w)
  S <- (sw %o% sw) * K
  S <- (S + t(S)) / 2
  
  ev_test <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev_test < -tol * max(1, abs(ev_test[1])))) {
    stop(sprintf(
      "Gewichtete Kovarianzmatrix ist nicht PSD (min. Eigenwert = %.4g).",
      min(ev_test)
    ), call. = FALSE)
  }
  
  ev  <- eigen(S, symmetric = TRUE)
  lam <- ev$values
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
    p.value    = if (!is.na(p_value)) as.numeric(p_value) else NA_real_,
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

# ----------------------------- Algorithmus 1: L2-Test ---------

#' @title Asymptotischer L2-Test für MCAR (ohne Parallelisierung)
#' @description Testet H0: µ_A = µ_B mittels T_{µ,L2} = n ||µ̂_A - µ̂_B||^2_{L2}.
#' @inheritParams .tfu_prepare_inputs
#' @inheritParams .tfu_subdomain_idx_paper
#' @param fve Fraction of variance explained (0–1) zur Wahl der Anzahl KL-Komponenten q.
#' @param B Anzahl Simulationen für die KL-Mischung (Monte-Carlo zur p-Wert-Approximation).
#' @param seed Optionaler RNG-Seed.
#' @param alpha Signifikanzniveau (nur für grobe CI-Angabe über MC-Varianz; optional).
#' @return Erweitertes `htest`-Objekt.
asym_mean_L2_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, B = 5000,
                              min_frac = 0.10, seed = 42, alpha = 0.05) {
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  T_L2     <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .tfu_kl_from_cov(K, g); lam <- KL$lam
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]
  
  Z <- matrix(rnorm(q * B), nrow = q)
  W <- colSums((Z^2) * lam_q)
  
  p <- (sum(W >= T_L2) + 1) / (length(W) + 1)
  ci <- c(mean(diff) - sd(W), mean(diff) + sd(W)); attr(ci,"conf.level") <- 1 - alpha
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  
  .tfu_make_htest_ext("T_{mu,L2}", T_L2, p,
                      "L2-Test (KL-Mischung)",
                      data_name,
                      estimate = c("mean(A)" = mean(muA), "mean(B)" = mean(muB),
                                   "diff(A-B)" = mean(diff)),
                      conf.int = ci,
                      parameter = c(q = q, m = length(g)))
}

# ----------------------------- Algorithmus 2: Supremum-Test (+ optionale Bänder) ----

#' @title Asymptotischer Supremums-Test für MCAR (mit optionalen Bändern)
#' @description Testet H0: µ_A = µ_B mittels T_{µ,D} = sqrt(n) ||µ̂_A - µ̂_B||_∞.
#'   GP-Approximation über KL-Basis; die Monte-Carlo-Supremums-Simulation liefert
#'   gleichzeitig die Quantile für **simultane Konfidenzbänder**.
#' @inheritParams asym_mean_L2_test
#' @param compute_bands Logisch: zusätzlich simultane Bänder berechnen?
#' @param bands_only Logisch: nur Band-Informationen (leichtgewichtige Liste) zurückgeben?
#' @return `htest` (mit Zusatzfeldern) oder – bei `bands_only=TRUE` – eine kurze Liste.
asym_mean_sup_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                               fve = 0.99, B = 5000,
                               min_frac = 0.10, seed = 42, alpha = 0.05,
                               compute_bands = TRUE, bands_only = FALSE) {
  prep <- .tfu_prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .tfu_available_means(X, O, group_A)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  T_D <- sqrt(n) * max(abs(diff))
  
  K  <- .tfu_corrected_cov(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .tfu_kl_from_cov(K, g); lam <- KL$lam; phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]; phi_q <- phi[, seq_len(q), drop = FALSE]
  A <- sweep(phi_q, 2, sqrt(lam_q), "*")  # m x q
  
  # Monte-Carlo GP-Supremum (vektorisiert)
  Z <- matrix(rnorm(q * B), nrow = q)     # q x B
  gp_vals <- A %*% Z                      # m x B
  W <- apply(abs(gp_vals), 2, max)
  
  # p-Wert
  p <- (sum(W >= T_D) + 1) / (length(W) + 1)
  
  # Optionale Bänder
  if (isTRUE(compute_bands)) {
    q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
    halfwidth <- q_alpha / sqrt(n)
    lower <- diff - halfwidth
    upper <- diff + halfwidth
    ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci,"conf.level") <- 1 - alpha
  } else {
    q_alpha <- NA_real_; lower <- upper <- NULL; ci <- NULL
  }
  
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB), "diff(A-B)" = mean(diff))
  
  if (isTRUE(bands_only)) {
    return(list(
      grid=g, diff=diff, lower=lower, upper=upper, q_alpha=q_alpha,
      idx=idx, min_frac_used=sub$min_frac_used, fallback=sub$fallback,
      alpha=alpha, estimate=est_vec, conf.int=ci
    ))
  }
  
  h <- .tfu_make_htest_ext("T_{mu,D}", T_D, p,
                           "Supremum-Test (GP-Approx; optional bands)",
                           data_name,
                           estimate = est_vec, conf.int = ci, parameter = c(q = q, m = length(g)))
  h$idx   <- idx
  h$grid  <- g
  h$diff  <- diff
  h$lower <- lower
  h$upper <- upper
  h$q_alpha <- q_alpha
  h$alpha <- alpha
  h$min_frac_used <- sub$min_frac_used
  h$fallback <- sub$fallback
  h
}

# -----------------------------------------------------------------------------
# Bootstrap-Tests (foreach + doRNG) & optionale Bootstrap-Bänder
# -----------------------------------------------------------------------------

#' @title Gruppenspezifische Bootstrap-Mittelwerte (intern)
#' @keywords internal
if (!exists(".tfu_boot_group_mean", mode = "function")) {
  .tfu_boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
    p_hat <- colSums(O_grp) / n_total
    p_hat <- pmax(p_hat, eps)
    numer <- colSums(replace(X_grp, is.na(X_grp), 0))
    numer / (n_total * p_hat)
  }
}

# ---- Persistentes Parallel-Backend (paket-/skript-intern) ---------------------
.tfu_par_env <- local({
  e <- new.env(parent = emptyenv())
  e$cl <- NULL
  e$old_threads <- NULL
  e$finalizer_set <- FALSE
  e
})

shutdown_parallel_tfu <- function() {
  if (!is.null(.tfu_par_env$cl)) {
    try(parallel::stopCluster(.tfu_par_env$cl), silent = TRUE)
    .tfu_par_env$cl <- NULL
  }
  olds <- .tfu_par_env$old_threads
  if (!is.null(olds)) {
    if (!is.na(olds[1])) Sys.setenv(OPENBLAS_NUM_THREADS = olds[1])
    if (!is.na(olds[2])) Sys.setenv(MKL_NUM_THREADS     = olds[2])
    if (!is.na(olds[3])) Sys.setenv(OMP_NUM_THREADS     = olds[3])
    .tfu_par_env$old_threads <- NULL
  }
  invisible(TRUE)
}

.tfu_ensure_backend <- function(manage_backend = c("auto","force_pool","sequential"),
                                ncpus = parallel::detectCores(logical = TRUE),
                                worker_blas_threads = 1L,
                                seed = 42,
                                parallel_flag = TRUE) {
  manage_backend <- match.arg(manage_backend)
  
  if (!isTRUE(parallel_flag) || identical(manage_backend, "sequential")) {
    foreach::registerDoSEQ()
    doRNG::registerDoRNG(seed)
    return(list(nworkers = 1L, used = "sequential"))
  }
  
  existing_workers <- tryCatch(foreach::getDoParWorkers(), error = function(e) 1L)
  existing_backend <- tryCatch(foreach::getDoParName(),    error = function(e) "doSEQ")
  has_external_backend <- isTRUE(existing_workers > 1L && existing_backend != "doSEQ")
  
  if (identical(manage_backend, "auto") && has_external_backend) {
    doRNG::registerDoRNG(seed)
    return(list(nworkers = existing_workers, used = paste0("external:", existing_backend)))
  }
  
  if (!is.null(.tfu_par_env$cl)) {
    doParallel::registerDoParallel(.tfu_par_env$cl)
    doRNG::registerDoRNG(seed)
    nworkers <- tryCatch(foreach::getDoParWorkers(), error = function(e) NA_integer_)
    return(list(nworkers = ifelse(is.na(nworkers), 1L, nworkers), used = "internal-reused"))
  }
  
  ncpus <- max(1L, min(as.integer(ncpus), parallel::detectCores(logical = TRUE)))
  cl <- parallel::makeCluster(ncpus)
  
  .tfu_par_env$old_threads <- Sys.getenv(c("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"), unset = NA)
  parallel::clusterCall(cl, function(k) {
    Sys.setenv(OPENBLAS_NUM_THREADS = as.character(k),
               MKL_NUM_THREADS      = as.character(k),
               OMP_NUM_THREADS      = as.character(k))
    NULL
  }, worker_blas_threads)
  
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed)
  .tfu_par_env$cl <- cl
  
  if (!isTRUE(.tfu_par_env$finalizer_set)) {
    reg.finalizer(.tfu_par_env, function(e) {
      if (!is.null(e$cl)) {
        try(parallel::stopCluster(e$cl), silent = TRUE)
        e$cl <- NULL
      }
      olds <- e$old_threads
      if (!is.null(olds)) {
        if (!is.na(olds[1])) Sys.setenv(OPENBLAS_NUM_THREADS = olds[1])
        if (!is.na(olds[2])) Sys.setenv(MKL_NUM_THREADS     = olds[2])
        if (!is.na(olds[3])) Sys.setenv(OMP_NUM_THREADS     = olds[3])
        e$old_threads <- NULL
      }
    }, onexit = TRUE)
    .tfu_par_env$finalizer_set <- TRUE
  }
  
  list(nworkers = ncpus, used = "internal-new")
}

#' @title Bootstrap-basierter Mean-Test (L2 oder Supremum) — foreach + doRNG (mit optionalen Bändern)
#'
#' @description Liefert Bootstrap-p-Werte sowie (optional) simultane Bänder.
#' @param fd Optional: `tfd`-Objekt. Wird gegenüber `X_obs` priorisiert.
#' @param X_obs Optional: Matrix (n x m) mit NAs.
#' @param groups Optional: Gruppierungsvektor (zwei Ausprägungen). Fehlt er, automatische Gruppierung.
#' @param observed_ratio Anteil in [0,1] für Auto-Gruppierung.
#' @param B Anzahl Bootstrap-Zyklen.
#' @param min_frac Minimaler Gruppenbeobachtungsanteil pro Zeitstelle.
#' @param alpha Signifikanzniveau für simultane Bänder (wenn `compute_bands = TRUE`).
#' @param parallel Logisch: Parallelausführung via `foreach`?
#' @param ncpus Anzahl Kerne für Cluster (nur bei `parallel = TRUE` relevant).
#' @param seed Optionaler RNG-Seed (reproduzierbar via doRNG).
#' @param stat Kennung der Teststatistik: `"L2"` oder `"D"` (Supremum).
#' @param compute_bands Logisch: zusätzlich simultane Bänder berechnen?
#' @param return_boot Logisch: Bootstrap-Statistiken dem Objekt beilegen?
#' @param chunk_size Bootstraps pro Task (Default: automatisch ~ B/(3*Workers), min. 50)
#' @param manage_backend Wie das Parallel-Backend gehandhabt wird: "auto", "force_pool", "sequential".
#' @param worker_blas_threads Anzahl BLAS/OpenMP-Threads **pro Worker** (nur gesetzt, wenn eigener Pool gebaut wird).
#' @param bands_only Logisch: nur Band-Informationen (leichtgewichtige Liste) zurückgeben?
#' @return `htest` (inkl. Zusatzfelder) oder – bei `bands_only=TRUE` – eine kurze Liste.
#' @import foreach doRNG doParallel
boot_mean_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                           B = 5000,
                           min_frac = 0.10, alpha = 0.05,
                           parallel = TRUE,
                           ncpus = parallel::detectCores(logical = TRUE),
                           seed = 42,
                           stat = c("L2", "D"),
                           compute_bands = TRUE,
                           return_boot = FALSE,
                           chunk_size = NULL,
                           manage_backend = c("auto","force_pool","sequential"),
                           worker_blas_threads = 1L,
                           bands_only = FALSE) {
  stat <- match.arg(stat)
  manage_backend <- match.arg(manage_backend)
  
  # ---- Inputs vorbereiten ----
  prep <- .tfu_prepare_inputs(fd = fd, X_obs = X_obs, groups = groups, observed_ratio = observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .tfu_subdomain_idx_paper(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  est  <- .tfu_available_means(X, O, group_A)
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
  
  # ---- Backend: einmalig initialisieren oder vorhandenes nutzen ----
  be <- .tfu_ensure_backend(manage_backend = manage_backend,
                            ncpus = ncpus,
                            worker_blas_threads = worker_blas_threads,
                            seed = seed,
                            parallel_flag = parallel)
  nworkers <- be$nworkers
  
  # ---- Chunking-Heuristik ----
  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size < 1L) {
    chunk_size <- max(50L, ceiling(B / (3L * max(1L, nworkers))))
  } else {
    chunk_size <- as.integer(chunk_size)
  }
  idx_chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
  
  # ---- Foreach über Chunks ----
  boot_mat <- foreach::foreach(
    ch = idx_chunks,
    .combine = rbind,
    .init = NULL,
    .inorder = FALSE
  ) %dorng% {
    out <- matrix(NA_real_, nrow = length(ch), ncol = 3L,
                  dimnames = list(NULL, c("L2","D","redraws")))
    for (ii in seq_along(ch)) {
      redraws <- 0L
      repeat {
        samp_idx <- sample.int(n, n, replace = TRUE)
        gA_samp  <- group_A[samp_idx]
        IA <- as.numeric(gA_samp); IB <- 1 - IA
        OA <- (O[samp_idx, , drop = FALSE] * IA)
        OB <- (O[samp_idx, , drop = FALSE] * IB)
        pA <- colSums(OA) / n
        pB <- colSums(OB) / n
        if (all(pA > 0) && all(pB > 0)) break
        redraws <- redraws + 1L
      }
      
      # zentrierte Werte im Resample
      Xs <- X_cent[samp_idx, , drop = FALSE]
      
      # Gruppe A
      XA <- Xs; XA[IB == 1, ] <- NA
      numA <- colSums(replace(XA, is.na(XA), 0))
      muA_b <- numA / (n * pA)
      
      # Gruppe B
      XB <- Xs; XB[IA == 1, ] <- NA
      numB <- colSums(replace(XB, is.na(XB), 0))
      muB_b <- numB / (n * pB)
      
      d_b <- muA_b - muB_b
      out[ii, ] <- c(
        L2 = n * sum((d_b^2) * w),
        D  = sqrt(n) * max(abs(d_b)),
        redraws = redraws
      )
    }
    out
  }
  
  # ---- Auswertung ----
  boot_L2 <- boot_mat[, "L2"]
  boot_D  <- boot_mat[, "D"]
  total_redraws <- sum(boot_mat[, "redraws"])
  
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
    ci <- c(mean(diff) - halfwidth, mean(diff) + halfwidth); attr(ci,"conf.level") <- 1 - alpha
    lower <- diff - halfwidth
    upper <- diff + halfwidth
  } else {
    ci <- NULL; q_alpha <- NA_real_; lower <- upper <- NULL
  }
  
  est_vec <- c("mean(A)" = mean(muA), "mean(B)" = mean(muB), "diff(A-B)" = mean(diff))
  data_name <- if (!is.null(fd)) "fd (tfd/tfd_irreg via tf_gather)" else "X_obs"
  
  if (isTRUE(bands_only)) {
    return(list(
      grid=g, diff=diff, lower=lower, upper=upper, q_alpha=q_alpha,
      idx=idx, min_frac_used=sub$min_frac_used, fallback=sub$fallback,
      alpha=alpha, estimate=est_vec, conf.int=ci
    ))
  }
  
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
  hobj$total_redraws <- total_redraws
  hobj$skip_rate <- total_redraws / (B + total_redraws)
  if (isTRUE(return_boot)) { hobj$boot_L2 <- boot_L2; hobj$boot_D <- boot_D }
  
  hobj
}