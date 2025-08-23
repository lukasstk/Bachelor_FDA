# ===============================================================
# Type-I-Error Simulation (Case 1: H0, z.B. MCAR) — INNERE Parallelisierung
#   - Tabelle: asym_L2, asym_D, boot_L2, boot_D
#   - runs = 5000, B_asym = 1000, B_boot = 1000
#   - n in {100, 250, 500}
#   - boot_mean_test() managt Backend & Chunking intern
# ===============================================================

suppressPackageStartupMessages({
  library(tidyfun); library(tf); library(tictoc)
  # doParallel/doRNG müssen NICHT explizit registriert werden – boot_mean_test() kümmert sich
})

## ---------------------- Datengenerator --------------------------
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

## ---------------------- 1 Replikation (innen parallel) ------------------------
# Hinweis: kein 'workers' Parameter mehr nötig
one_rep_inner <- function(n, grid_spacing, mechanism, B_asym, B_boot, alpha, seed_for_rep,
                          ncpus_boot = parallel::detectCores(logical = TRUE) - 1L) {
  set.seed(seed_for_rep)
  X_obs <- simulate_X_obs(n, grid_spacing, mechanism)
  
  # Asymptotische Tests (seriell)
  res_L2 <- asym_mean_L2_test(X_obs = X_obs, B = B_asym, seed = seed_for_rep + 1L)
  res_D  <- asym_mean_sup_test(X_obs = X_obs, B = B_asym, seed = seed_for_rep + 2L)
  
  # Bootstrap-Tests – intern parallel (manage_backend="auto"), mit Chunk-Progress
  res_BT_L2 <- boot_mean_test(
    X_obs = X_obs, B = B_boot, stat = "L2",
    parallel = TRUE, ncpus = ncpus_boot, seed = seed_for_rep + 11L,
    manage_backend = "auto"  # zeigt Chunk-Progress
  )
  res_BT_D  <- boot_mean_test(
    X_obs = X_obs, B = B_boot, stat = "D",
    parallel = TRUE, ncpus = ncpus_boot, seed = seed_for_rep + 22L,
    manage_backend = "auto"
  )
  
  c(
    asym_L2 = as.numeric(res_L2$p.value < alpha),
    asym_D  = as.numeric(res_D$p.value  < alpha),
    boot_L2 = as.numeric(res_BT_L2$p.value < alpha),
    boot_D  = as.numeric(res_BT_D$p.value  < alpha)
  )
}

## ---------------------- R Replikationen (außen seriell) -----------------------
run_type1_inner <- function(n, runs = 5000, grid_spacing = 100, mechanism = "MCAR",
                            alpha = 0.05, B_asym = 1000, B_boot = 1000,
                            base_seed = 1,
                            ncpus_boot = parallel::detectCores(logical = TRUE) - 1L) {
  
  hits <- matrix(0, nrow = runs, ncol = 4)
  colnames(hits) <- c("asym_L2","asym_D","boot_L2","boot_D")
  
  pb <- txtProgressBar(min = 0, max = runs, style = 3)
  for (r in seq_len(runs)) {
    hits[r, ] <- one_rep_inner(
      n, grid_spacing, mechanism, B_asym, B_boot, alpha,
      seed_for_rep = base_seed + 1000L * r,
      ncpus_boot   = ncpus_boot
    )
    setTxtProgressBar(pb, r)  # Run-Progress (jede Iteration)
  }
  close(pb)
  
  colMeans(hits)
}

## ---------------------- Ausführen: n = 100, 250, 500 -------------------------
set.seed(42)

tic("inner-parallel n = 100")
res100 <- run_type1_inner(n = 100, runs = 5000, B_asym = 10000, B_boot = 10000)
toc()

tic("inner-parallel n = 250")
res250 <- run_type1_inner(n = 250, runs = 5000, B_asym = 10000, B_boot = 10000)
toc()

tic("inner-parallel n = 500")
res500 <- run_type1_inner(n = 500, runs = 5000, B_asym = 10000, B_boot = 10000)
toc()

type1_table <- rbind(
  data.frame(n = 100, t(res100)),
  data.frame(n = 250, t(res250)),
  data.frame(n = 500, t(res500))
)
print(type1_table, row.names = FALSE)




# > res100 <- run_type1_inner(n = 100, runs = 2000, B_asym = 2000, B_boot = 2000)
# |======================================================================================| 100%
# Es gab 50 oder mehr Warnungen (Anzeige der ersten 50 mit warnings())
# > toc()
# inner-parallel n = 100: 556.42 sec elapsed
# > 
#   > tic("inner-parallel n = 250")
# > res250 <- run_type1_inner(n = 250, runs = 2000, B_asym = 2000, B_boot = 2000)
# |======================================================================================| 100%
# Es gab 50 oder mehr Warnungen (Anzeige der ersten 50 mit warnings())
# > toc()
# inner-parallel n = 250: 929.98 sec elapsed
# > 
#   > tic("inner-parallel n = 500")
# > res500 <- run_type1_inner(n = 500, runs = 2000, B_asym = 2000, B_boot = 2000)
# |======================================================================================| 100%
# Es gab 50 oder mehr Warnungen (Anzeige der ersten 50 mit warnings())
# > toc()
# inner-parallel n = 500: 1602.88 sec elapsed
# > 
#   > type1_table <- rbind(
#     +   data.frame(n = 100, t(res100)),
#     +   data.frame(n = 250, t(res250)),
#     +   data.frame(n = 500, t(res500))
#     + )
# > print(type1_table, row.names = FALSE)
# n asym_L2 asym_D boot_L2 boot_D
# 100  0.0950 0.1225  0.0835  0.086
# 250  0.0555 0.0765  0.0455  0.054
# 500  0.0610 0.0885  0.0575  0.057
