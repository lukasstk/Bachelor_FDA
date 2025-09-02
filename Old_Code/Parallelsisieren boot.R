#' Tidyfun extension: Algorithm 5 (bootstrap for L2 and Supremum) — X_obs-only mit boot::boot
#'
#' Erwartet NUR `X_obs` (numeric Matrix n x m mit NAs).
#' Intern: O_mat, group_A und grid automatisch.
#'
#' Parallelisierung:
#'   - parallel = TRUE (default) → nutzt alle verfügbaren Kerne
#'       * Windows: "snow" (PSOCK)
#'       * Linux/macOS: "multicore"
#'   - parallel = FALSE → seriell
tfu_algo5_bootstrap <- function(X_obs,
                                B = 5000, eps = 1e-8,
                                min_frac = 0.10, alpha = 0.05,
                                parallel = TRUE,
                                ncpus = parallel::detectCores(logical = FALSE),
                                seed = NULL,
                                return_boot = FALSE) {
  # --- Inputs vorbereiten (Hilfsfunktionen gelten als gegeben) ---
  prep <- .tfu_prepare_inputs(X_obs)
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
  
  # schnelle Trapez-Integration (statt tf::tfd + tf_integrate)
  w <- .tfu_trap_weights(g)
  T_L2 <- n * sum((diff^2) * w)
  T_D  <- sqrt(n) * max(abs(diff))
  
  # --- Residuen zentrieren ---
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  # --- Bootstrap-Statistik (ohne externe Helper im Worker) ---
  bootfun <- function(data, indices) {
    # data = X_cent (vom master übergeben); O, group_A, n, w, eps via Closure
    samp_idx <- indices
    sampA <- samp_idx[group_A[samp_idx] == 1L]
    sampB <- samp_idx[group_A[samp_idx] == 0L]
    
    # mu_A^* und mu_B^* per available-data-Formel im Bootstrap
    if (length(sampA)) {
      OA <- O[sampA,, drop = FALSE]
      XA <- data[sampA,, drop = FALSE]
      pA <- pmax(colSums(OA) / n, eps)
      numA <- colSums(replace(XA, is.na(XA), 0))
      muA_b <- numA / (n * pA)
    } else muA_b <- rep(0, ncol(data))
    
    if (length(sampB)) {
      OB <- O[sampB,, drop = FALSE]
      XB <- data[sampB,, drop = FALSE]
      pB <- pmax(colSums(OB) / n, eps)
      numB <- colSums(replace(XB, is.na(XB), 0))
      muB_b <- numB / (n * pB)
    } else muB_b <- rep(0, ncol(data))
    
    d_b <- muA_b - muB_b
    c("L2" = n * sum((d_b^2) * w),
      "D"  = sqrt(n) * max(abs(d_b)))
  }
  
  # --- RNG für Reproduzierbarkeit (parallel-freundlich) ---
  if (!is.null(seed)) {
    set.seed(seed)
    suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
  }
  
  # --- Bootstrapping (plattformabhängig) ---
  os <- .Platform$OS.type
  if (parallel && os == "windows") {
    # Snow/PSOCK
    cl <- parallel::makeCluster(max(1L, ncpus), type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    # 'boot' serialisiert die Funktion + Closure; kein zusätzliches Export nötig
    res_boot <- boot::boot(
      data = X_cent,
      statistic = bootfun,
      R = B,
      parallel = "snow",
      ncpus = max(1L, ncpus),
      cl = cl
    )
  } else if (parallel && os != "windows") {
    # Fork / Multicore
    res_boot <- boot::boot(
      data = X_cent,
      statistic = bootfun,
      R = B,
      parallel = "multicore",
      ncpus = max(1L, ncpus)
    )
  } else {
    # Seriell
    res_boot <- boot::boot(
      data = X_cent,
      statistic = bootfun,
      R = B,
      parallel = "no"
    )
  }
  
  # --- Auswertung ---
  boot_L2 <- if (!is.null(colnames(res_boot$t))) res_boot$t[, "L2"] else res_boot$t[, 1]
  boot_D  <- if (!is.null(colnames(res_boot$t))) res_boot$t[, "D"]  else res_boot$t[, 2]
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (B + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (B + 1)
  
  out <- list(
    grid = g, muA = muA, muB = muB, diff = diff,
    T_L2 = as.numeric(T_L2), T_D = as.numeric(T_D),
    p_L2 = p_L2, p_D = p_D,
    idx = idx, min_frac_used = sub$min_frac_used, fallback = sub$fallback,
    alpha = alpha
  )
  
  # Konfidenzbänder (Algorithmus 7)
  q_alpha   <- as.numeric(stats::quantile(boot_D, probs = 1 - alpha, names = FALSE))
  halfwidth <- q_alpha / sqrt(n)
  out$lower   <- diff - halfwidth
  out$upper   <- diff + halfwidth
  out$q_alpha <- q_alpha
  
  if (isTRUE(return_boot)) {
    out$boot_L2 <- boot_L2
    out$boot_D  <- boot_D
    out$boot    <- res_boot
  }
  
  out
}





# ===============================================================
# Vergleich: Algorithmus 5 (Bootstrap) seriell vs. parallel
# ===============================================================

library(tidyfun)
library(tf)

# ---- Beispiel-Daten simulieren ----
set.seed(123)
n <- 100
grid <- seq(0, 1, length.out = 50)

simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}

bm_mat <- t(replicate(n, simulate_bm(grid)))
O_mat  <- matrix(1L, n, length(grid))  # alles beobachtet
X_obs  <- bm_mat
X_obs[O_mat == 0L] <- NA_real_

# ---- Seriell ----
cat("Starte seriellen Bootstrap...\n")
time_seq <- system.time({
  res_seq <- tfu_algo5_bootstrap(X_obs = X_obs, B = 50000, parallel = FALSE)
})
print(time_seq)

# ---- Parallel ----
cat("Starte parallelen Bootstrap...\n")
time_par <- system.time({
  res_par <- tfu_algo5_bootstrap(X_obs = X_obs, B = 50000, parallel = TRUE)
})
print(time_par)

# ---- Vergleich der Ergebnisse ----
cat("\nVergleich der Ergebnisse:\n")
cat("Seriell:   p_L2 =", res_seq$p_L2, " p_D =", res_seq$p_D, "\n")
cat("Parallel:  p_L2 =", res_par$p_L2, " p_D =", res_par$p_D, "\n")
