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
  library(foreach)
  library(doParallel)
  library(doRNG) 
  library(pbapply)
})

## -------------------- Backend EINMAL starten & registrieren -------------------
# workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
# cl <- parallel::makeCluster(workers)
# on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
# doParallel::registerDoParallel(cl)
# doRNG::registerDoRNG(42)
# 
# ## Gegen Oversubscription (BLAS/OpenMP) – wichtig für Performance
# old_threads <- Sys.getenv(c("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"),
#                           unset = NA)
# on.exit({
#   if (!is.na(old_threads[1])) Sys.setenv(OPENBLAS_NUM_THREADS = old_threads[1])
#   if (!is.na(old_threads[2])) Sys.setenv(MKL_NUM_THREADS     = old_threads[2])
#   if (!is.na(old_threads[3])) Sys.setenv(OMP_NUM_THREADS     = old_threads[3])
# }, add = TRUE)
# Sys.setenv(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
# 
# message(sprintf("Inneres Backend: %s mit %d Workern",
#                 foreach::getDoParName(), foreach::getDoParWorkers()))

## ---------------------- Datengenerator (fix) --------------------------
simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}

# MCAR: mit Wkt. 0.5 komplette Kurve, sonst ein Intervall [L, U)
# (wie §5.2 – KEIN zusätzliches p_complete mehr!)
simulate_X_obs <- function(n, grid_spacing = 100,
                           mechanism = c("MCAR","MNAR")) {
  mechanism <- match.arg(mechanism)
  grid <- seq(0, 1, length.out = grid_spacing)
  m <- length(grid)
  
  # Brownian Motion
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  
  if (mechanism == "MCAR") {
    complete <- runif(n) < 0.5                # A: komplett, B: Intervall
    O_mat <- matrix(0L, n, m)
    if (any(complete)) O_mat[complete, ] <- 1L
    if (any(!complete)) {
      k  <- sum(!complete)
      U  <- matrix(runif(k * 2), ncol = 2)
      L  <- pmin(U[,1], U[,2]); R <- pmax(U[,1], U[,2])
      Oi <- outer(L, grid, "<=") & outer(R, grid, ">") 
      O_mat[!complete, ] <- 1L * Oi
    }
  } else {  # MNAR: beobachtet, wenn -1 < X(t) < 2
    within   <- (bm_mat > -1) & (bm_mat < 2)  # n x m
    complete <- rowSums(within) == m
    O_mat <- matrix(0L, n, m)
    if (any(complete))   O_mat[complete, ]   <- 1L
    if (any(!complete))  O_mat[!complete, ]  <- 1L * within[!complete, ]
  }
  
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  X_obs
}

## ---------------------- 1 Replikation (innen parallel) ------------------------
one_rep_inner <- function(n, grid_spacing, mechanism, B_asym, B_boot, alpha,
                          ncpus_boot = parallel::detectCores(logical = TRUE) - 1L) {
  X_obs <- simulate_X_obs(n, grid_spacing, mechanism)
  
  groups <- rowMeans(!is.na(X_obs)) == 1
  res_L2 <- asym_mean_L2_test(X_obs=X_obs, groups=groups, B=B_asym, min_frac=0.10)
  res_D  <- asym_mean_sup_test(X_obs=X_obs, groups=groups, B=B_asym, min_frac=0.10)
  res_BT_L2 <- boot_mean_test(X_obs=X_obs, groups=groups, B=B_boot, stat="L2",
                              parallel=TRUE, ncpus=ncpus_boot, manage_backend="auto",
                              min_frac=0.10)
  res_BT_D  <- boot_mean_test(X_obs=X_obs, groups=groups, B=B_boot, stat="D",
                              parallel=TRUE, ncpus=ncpus_boot, manage_backend="auto",
                              min_frac=0.10)
  
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
                            ncpus_boot = parallel::detectCores(logical = TRUE) - 1L) {
  # pbapply-Progress (Timer-Anzeige)
  pboptions(type = "timer")
  
  # pbreplicate ruft one_rep_inner 'runs'-mal auf und zeigt Progress
  res_arr <- pbreplicate(
    runs,
    one_rep_inner(n, grid_spacing, mechanism, B_asym, B_boot, alpha,
                  ncpus_boot = ncpus_boot)
  )
  # res_arr ist 4 x runs; transponieren zu runs x 4 wie zuvor
  hits <- t(res_arr)
  colnames(hits) <- c("asym_L2","asym_D","boot_L2","boot_D")
  
  colMeans(hits)
}

## ---------------------- Ausführen: n = 100, 250, 500 -------------------------
set.seed(42)

tic("inner-parallel n = 100")
res100 <- run_type1_inner(n = 100, runs = 2000, B_asym = 2000, B_boot = 2000)
toc()

tic("inner-parallel n = 250")
res250 <- run_type1_inner(n = 250, runs = 2000, B_asym = 2000, B_boot = 2000)
toc()

tic("inner-parallel n = 500")
res500 <- run_type1_inner(n = 500, runs = 2000, B_asym = 2000, B_boot = 2000)
toc()

type1_table <- rbind(
  data.frame(n = 100, t(res100)),
  data.frame(n = 250, t(res250)),
  data.frame(n = 500, t(res500))
)
print(type1_table, row.names = FALSE)



# > ## ---------------------- Ausführen: n = 100, 250, 500 -------------------------
#   > tic("inner-parallel n = 100")
# > res100 <- run_type1_inner(n = 100, runs = 2000, B_asym = 2000, B_boot = 2000)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13m 19s
# > toc()
# inner-parallel n = 100: 799.25 sec elapsed
# > 
#   > tic("inner-parallel n = 250")
# > res250 <- run_type1_inner(n = 250, runs = 2000, B_asym = 2000, B_boot = 2000)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=26m 07s
# > toc()
# inner-parallel n = 250: 1567.26 sec elapsed
# > 
#   > tic("inner-parallel n = 500")
# > res500 <- run_type1_inner(n = 500, runs = 2000, B_asym = 2000, B_boot = 2000)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=50m 06s
# > toc()
# inner-parallel n = 500: 3005.72 sec elapsed
# > 
#   > type1_table <- rbind(
#     +   data.frame(n = 100, t(res100)),
#     +   data.frame(n = 250, t(res250)),
#     +   data.frame(n = 500, t(res500))
#     + )
# > print(type1_table, row.names = FALSE)
# n asym_L2 asym_D boot_L2 boot_D
# 100  0.0675 0.0730  0.0600 0.0500
# 250  0.0560 0.0530  0.0505 0.0495
# 500  0.0660 0.0665  0.0685 0.0630
