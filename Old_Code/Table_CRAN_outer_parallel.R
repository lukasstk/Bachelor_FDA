# ===============================================================
# Type-I-Error Simulation (Case 1: H0, z.B. MCAR) — Äußere Parallelisierung
#   - Tabelle: asym_L2, asym_D, boot_L2, boot_D
#   - für n in {100, 250, 500}
# Benötigt:
#   - asym_mean_L2_test, asym_mean_sup_test, boot_mean_test (aus deinem Paket)
#   - Pakete: tidyfun, tf, tictoc, foreach, doParallel, doRNG
# ===============================================================

suppressPackageStartupMessages({
  library(tidyfun); library(tf); library(tictoc)
  library(foreach); library(doParallel); library(doRNG)
  library(progressr)
})

# ----- Brownian + Missingness (H0: kein Mean-Shift) ---------------------------
simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}
make_O_mcar <- function(grid) {
  if (runif(1) < 0.5) rep(1L, length(grid))
  else { U <- sort(runif(2)); as.integer(grid >= U[1] & grid < U[2]) }
}
make_O_mnar <- function(x_row) {
  if (all(x_row > -1 & x_row < 2)) rep(1L, length(x_row)) else as.integer(x_row > -1 & x_row < 2)
}

simulate_X_obs <- function(n, grid_spacing = 100, mechanism = c("MCAR","MNAR"), p_complete = 0.5) {
  mechanism <- match.arg(mechanism)
  grid <- seq(0, 1, length.out = grid_spacing)
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  I_A <- rbinom(n, 1, p_complete)
  O_mat <- matrix(0L, n, grid_spacing)
  for (i in seq_len(n)) {
    if (I_A[i] == 1L) O_mat[i, ] <- 1L
    else O_mat[i, ] <- if (mechanism == "MCAR") make_O_mcar(grid) else make_O_mnar(bm_mat[i, ])
  }
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  X_obs
}

# ----- 1 Replikat -------------------------------------------------------------
# ACHTUNG: Bootstraps INNEN jetzt SERIELL, um nested parallelism zu vermeiden!
one_rep <- function(n, grid_spacing, mechanism, B_asym, B_boot, alpha, seed_for_rep) {
  set.seed(seed_for_rep)
  X_obs <- simulate_X_obs(n, grid_spacing, mechanism)
  
  # Asymptotische Tests
  res_L2 <- asym_mean_L2_test(X_obs = X_obs, B = B_asym)
  res_D  <- asym_mean_sup_test(X_obs = X_obs, B = B_asym)
  
  # Bootstrap-Tests (innen SERIELL!)
  res_BT_L2 <- boot_mean_test(X_obs = X_obs, B = B_boot, stat = "L2",
                              parallel = FALSE, ncpus = 1)
  res_BT_D  <- boot_mean_test(X_obs = X_obs, B = B_boot, stat = "D",
                              parallel = FALSE, ncpus = 1)
  
  c(
    asym_L2 = as.numeric(res_L2$p.value < alpha),
    asym_D  = as.numeric(res_D$p.value  < alpha),
    boot_L2 = as.numeric(res_BT_L2$p.value < alpha),
    boot_D  = as.numeric(res_BT_D$p.value  < alpha)
  )
}

# ----- Typ-I-Fehler für ein n schätzen (ÄUSSERER PARALLELISMUS mit Progress) -----
run_type1 <- function(n, R = 500, grid_spacing = 100, mechanism = "MCAR",
                      alpha = 0.05, B_asym = 2000, B_boot = 1000,
                      base_seed = 1,
                      ncpus_outer = max(1L, parallel::detectCores(logical = FALSE) - 1L)) {
  
  # Cluster starten
  cl <- parallel::makeCluster(ncpus_outer)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # Pakete auf den Workern laden
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(tf); library(tidyfun)
    })
    NULL
  })
  
  # Alle benötigten Funktionen exportieren
  funcs_to_export <- c(
    "one_rep","simulate_X_obs","simulate_bm","make_O_mcar","make_O_mnar",
    "asym_mean_L2_test","asym_mean_sup_test","asym_conf_bands",
    "boot_mean_test","boot_conf_bands",
    ".tfu_groups_to_logical",".tfu_group_from_delta",".tfu_from_fd",
    ".tfu_subdomain_idx_paper",".tfu_prepare_inputs",
    ".tfu_trap_weights",".tfu_available_means",".tfu_corrected_cov",
    ".tfu_kl_from_cov",".tfu_make_htest_ext",".tfu_boot_group_mean"
  )
  parallel::clusterExport(cl, varlist = funcs_to_export, envir = .GlobalEnv)
  
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(base_seed)
  
  message(sprintf("Starte %d Replikationen auf %d Kernen …", R, ncpus_outer))
  
  # Fortschrittsbalken vorbereiten
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  progress_env <- new.env()
  progress_env$counter <- 0
  
  M <- foreach::foreach(
    r = seq_len(R),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("tf","tidyfun"),
    .export = funcs_to_export
  ) %dorng% {
    res <- one_rep(n, grid_spacing, mechanism, B_asym, B_boot, alpha, base_seed + 1000 * r)
    
    # Fortschritt hochzählen
    progress_env$counter <- progress_env$counter + 1
    setTxtProgressBar(pb, progress_env$counter)
    
    res
  }
  
  close(pb)
  
  est <- colMeans(M)
  data.frame(
    n = n,
    asym_L2 = round(est["asym_L2"], 2),
    asym_D  = round(est["asym_D"],  2),
    boot_L2 = round(est["boot_L2"], 2),
    boot_D  = round(est["boot_D"],  2),
    row.names = NULL
  )
}

# ======= AUSFÜHREN: n = 100, 250, 500 (Case 1: MCAR) =========================
set.seed(42)
ncpus_outer <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

tic("n = 100 (outer parallel)")
res100 <- run_type1(n = 100, R = 100, mechanism = "MCAR",
                    B_asym = 10000, B_boot = 10000,
                    ncpus_outer = ncpus_outer)
toc()

tic("n = 250 (outer parallel)")
res250 <- run_type1(n = 250, R = 100, mechanism = "MCAR",
                    B_asym = 10000, B_boot = 10000,
                    ncpus_outer = ncpus_outer)
toc()

tic("n = 500 (outer parallel)")
res500 <- run_type1(n = 500, R = 100, mechanism = "MCAR",
                    B_asym = 10000, B_boot = 10000,
                    ncpus_outer = ncpus_outer)
toc()

type1_table <- rbind(res100, res250, res500)
print(type1_table, row.names = FALSE)


# n = 100 (outer parallel): 256.35 sec elapsed
# > 
#   > tic("n = 250 (outer parallel)")
# > res250 <- run_type1(n = 250, R = 100, mechanism = "MCAR",
#                       +                     B_asym = 10000, B_boot = 10000,
#                       +                     ncpus_outer = ncpus_outer)
# Starte 100 Replikationen auf 7 Kernen …
# |                                                                                      |   0%
# > toc()
# n = 250 (outer parallel): 520.97 sec elapsed
# > 
#   > tic("n = 500 (outer parallel)")
# > res500 <- run_type1(n = 500, R = 100, mechanism = "MCAR",
#                       +                     B_asym = 10000, B_boot = 10000,
#                       +                     ncpus_outer = ncpus_outer)
# Starte 100 Replikationen auf 7 Kernen …
# |                                                                                      |   0%
# > toc()
# n = 500 (outer parallel): 701.48 sec elapsed
# > 
#   > type1_table <- rbind(res100, res250, res500)
# > print(type1_table, row.names = FALSE)
# n asym_L2 asym_D boot_L2 boot_D
# 100    0.09   0.15    0.06   0.09
# 250    0.04   0.13    0.04   0.10
# 500    0.05   0.03    0.05   0.02
