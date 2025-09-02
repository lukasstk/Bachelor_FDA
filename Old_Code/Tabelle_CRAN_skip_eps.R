# ===============================================================
# Type-I-Error Simulation (Case 1: H0, z.B. MCAR)
#   - liefert Tabelle mit: asym_L2, asym_D, boot_L2, boot_D
#   - für n in {100, 250, 500}
# Benötigt:
#   - asym_mean_L2_test, asym_mean_sup_test, boot_mean_test (aus deinem Paket)
#   - Pakete: tidyfun, tf, tictoc (optional)
# ===============================================================

suppressPackageStartupMessages({
  library(tidyfun); library(tf); library(tictoc)
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
one_rep <- function(n, grid_spacing, mechanism, B_asym, B_boot, alpha, seed_for_rep) {
  set.seed(seed_for_rep)
  X_obs <- simulate_X_obs(n, grid_spacing, mechanism)
  
  # Asymptotische Tests
  res_L2 <- asym_mean_L2_test(X_obs = X_obs, B = B_asym)
  res_D  <- asym_mean_sup_test(X_obs = X_obs, B = B_asym)
  
  # Bootstrap-Tests (parallelisiert innen)
  res_BT_L2 <- boot_mean_test(X_obs = X_obs, B = B_boot, stat = "L2",parallel = TRUE)
  res_BT_D  <- boot_mean_test(X_obs = X_obs, B = B_boot, stat = "D", parallel = TRUE)
  
  c(
    asym_L2 = as.numeric(res_L2$p.value < alpha),
    asym_D  = as.numeric(res_D$p.value  < alpha),
    boot_L2 = as.numeric(res_BT_L2$p.value < alpha),
    boot_D  = as.numeric(res_BT_D$p.value  < alpha)
  )
}

# ----- Typ-I-Fehler für ein n schätzen (mit Progress-Bar) ---------------------
run_type1 <- function(n, R = 500, grid_spacing = 100, mechanism = "MCAR",
                      alpha = 0.05, B_asym = 2000, B_boot = 1000, base_seed = 1) {
  M <- matrix(0, nrow = R, ncol = 4)
  colnames(M) <- c("asym_L2","asym_D","boot_L2","boot_D")
  
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  for (r in seq_len(R)) {
    M[r, ] <- one_rep(n, grid_spacing, mechanism, B_asym, B_boot, alpha, base_seed + 1000*r)
    setTxtProgressBar(pb, r)
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

tic("n = 100")
res100 <- run_type1(n = 100, R = 100, mechanism = "MCAR", B_asym = 10000, B_boot = 10000)
toc()

tic("n = 250")
res250 <- run_type1(n = 250, R = 100, mechanism = "MCAR", B_asym = 10000, B_boot = 10000)
toc()

tic("n = 500")
res500 <- run_type1(n = 500, R = 100, mechanism = "MCAR", B_asym = 10000, B_boot = 10000)
toc()

type1_table <- rbind(res100, res250, res500)
type1_table

#SKIP-SYM with R = 100, B_asym = 2000, B_boot = 2000
# > type1_table <- rbind(res100, res250, res500)
# > type1_table
#     n asym_L2 asym_D boot_L2 boot_D
# 1 100    0.09   0.15    0.06   0.08
# 2 250    0.05   0.13    0.04   0.10
# 3 500    0.05   0.03    0.04   0.03
#         with R = 100, B_asym = 10000, B_boot = 10000
# > type1_table <- rbind(res100, res250, res500)
# > type1_table
#     n asym_L2 asym_D boot_L2 boot_D
# 1 100    0.09   0.15    0.06   0.09
# 2 250    0.05   0.13    0.04   0.10
# 3 500    0.05   0.03    0.05   0.02
