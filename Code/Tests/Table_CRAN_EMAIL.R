# ===============================================================
# Type-I-Error Simulation (Case 1: H0, z.B. MCAR)
#   - liefert Tabelle mit: asym_L2, asym_D, boot_L2, boot_D
#   - für n in {100, 250}
# ===============================================================

library(tidyfun); library(tf)

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
  ids <- paste0("ID", seq_len(n))
  I_A <- rbinom(n, 1, p_complete)  # 1 = komplett, 0 = lückenhaft
  O_mat <- matrix(0L, n, grid_spacing)
  for (i in seq_len(n)) {
    if (I_A[i] == 1L) O_mat[i, ] <- 1L
    else O_mat[i, ] <- if (mechanism == "MCAR") make_O_mcar(grid) else make_O_mnar(bm_mat[i, ])
  }
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  list(X_obs = X_obs, grid = grid, ids = ids)
}

# ----- 1 Replikat: vier Tests (asym L2/D, boot L2/D) -------------------------
one_rep <- function(n, grid_spacing, mechanism, B_asym, B_boot, alpha, base_seed) {
  set.seed(base_seed)
  sim <- simulate_X_obs(n, grid_spacing, mechanism)
  X_obs <- sim$X_obs; grid <- sim$grid; ids <- sim$ids
  
  # fd-Objekt für asymptotische Tests
  fd <- tf::tfd(X_obs, arg = grid, id = ids)
  
  # asymptotische Tests
  res_L2 <- tfu_algo1_L2_test(fd = fd, B = B_asym)
  res_D  <- tfu_algo2_sup_test(fd = fd, B = B_asym)
  
  # Bootstrap-Tests (bleiben matrix-basiert)
  res_BT <- tfu_algo5_bootstrap(X_obs = X_obs, B = B_boot)
  
  c(
    asym_L2 = as.numeric(res_L2$p.value < alpha),
    asym_D  = as.numeric(res_D$p.value  < alpha),
    boot_L2 = as.numeric(res_BT$p_L2    < alpha),
    boot_D  = as.numeric(res_BT$p_D     < alpha)
  )
}

# ----- Typ-I-Fehler für ein n schätzen ---------------------------------------
run_type1 <- function(n, R = 500, grid_spacing = 100, mechanism = "MCAR",
                      alpha = 0.05, B_asym = 2000, B_boot = 1000, base_seed = 1,
                      quiet = TRUE) {
  if (quiet) {
    old_w <- getOption("warn"); options(warn = 1)
    on.exit(options(warn = old_w), add = TRUE)
  }
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  M <- matrix(0, nrow = R, ncol = 4)
  colnames(M) <- c("asym_L2","asym_D","boot_L2","boot_D")
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

# ======= AUSFÜHREN: n = 100 und 250 (Case 1: MCAR) ============================
set.seed(42)
res100 <- suppressWarnings(run_type1(n = 100, R = 100, mechanism = "MCAR", B_asym = 5000, B_boot = 5000))
res250 <- suppressWarnings(run_type1(n = 250, R = 100, mechanism = "MCAR", B_asym = 5000, B_boot = 5000))

type1_table <- rbind(res100, res250)
type1_table
