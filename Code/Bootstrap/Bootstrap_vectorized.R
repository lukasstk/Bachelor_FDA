# ===============================================================
#  Typ-I-Fehler unter H0 für beide Tests (L2 & Supremum), MCAR
#  — Bootstrap gemäß Algorithmus 5 (mit tf_integrate)
# ===============================================================

library(tidyfun)
library(dplyr)
library(tidyr)

# ------------- Grundeinstellungen ------------------------------
set.seed(42)
n            <- 100            # Pfade je Datensatz
grid_spacing <- 100
grid_full    <- seq(0, 1, length.out = grid_spacing)
nsim         <- 100            # Null-Simulationen
Bstar        <- 5000           # Bootstrap-Draws
alpha        <- 0.05
eps_p        <- 1e-8           # Clipping für p-Hüte in Bootstraps

# ------------- Hilfsfunktionen ---------------------------------
simulate_bm <- function(grid) {
  dt <- diff(grid)[1]
  c(0, cumsum(rnorm(length(grid) - 1, sd = sqrt(dt))))
}

make_O_mcar <- function(grid) {
  if (runif(1) < 0.5) {
    rep(1, length(grid))
  } else {
    U <- sort(runif(2))
    as.numeric(grid >= U[1] & grid < U[2])
  }
}

boot_group_mean <- function(X_grp, O_grp, n_total, eps = 1e-8) {
  p_hat <- colSums(O_grp) / n_total
  p_hat <- pmax(p_hat, eps)
  numer <- colSums(replace(X_grp, is.na(X_grp), 0))
  numer / (n_total * p_hat)
}

# Bootstrap (Algo 5) mit tf_integrate für L2
bootstrap_both <- function(X_obs, O_mat, I_A, I_B,
                           muA_hat, muB_hat, grid, n_total, Bstar, eps = 1e-8) {
  
  A_idx <- which(I_A == 1); nA <- length(A_idx)
  B_idx <- which(I_B == 1); nB <- length(B_idx)
  
  # Zentrierung
  X_cent <- X_obs
  X_cent[A_idx, ] <- sweep(X_obs[A_idx, ], 2, muA_hat, `-`)
  X_cent[B_idx, ] <- sweep(X_obs[B_idx, ], 2, muB_hat, `-`)
  
  T_boot_L2 <- numeric(Bstar)
  T_boot_D  <- numeric(Bstar)
  
  for (b in seq_len(Bstar)) {
    sampA <- sample(A_idx, nA, replace = TRUE)
    sampB <- sample(B_idx, nB, replace = TRUE)
    
    muA_b <- boot_group_mean(X_cent[sampA, , drop = FALSE],
                             O_mat[sampA, , drop = FALSE], n_total, eps)
    muB_b <- boot_group_mean(X_cent[sampB, , drop = FALSE],
                             O_mat[sampB, , drop = FALSE], n_total, eps)
    
    diff_b <- muA_b - muB_b
    
    # L2 mit tf_integrate
    diff_tfd <- tfd(matrix(diff_b, nrow = 1), arg = grid)
    T_boot_L2[b] <- n_total * tf_integrate(diff_tfd^2, arg = grid)
    
    # Supremum
    T_boot_D[b] <- sqrt(n_total) * max(abs(diff_b))
  }
  
  list(L2 = T_boot_L2, D = T_boot_D)
}

# ------------- Simulation --------------------------------------
p_vals_L2 <- numeric(nsim)
p_vals_D  <- numeric(nsim)
pb <- txtProgressBar(min = 0, max = nsim, style = 3)

for (sim in seq_len(nsim)) {
  
  ids    <- paste0("ID", seq_len(n))
  bm_mat <- t(replicate(n, simulate_bm(grid_full)))
  O_mat  <- t(sapply(ids, \(.) make_O_mcar(grid_full)))
  rownames(bm_mat) <- ids; colnames(O_mat) <- round(grid_full, 3)
  
  I_A <- as.numeric(rowSums(O_mat) == grid_spacing)
  I_B <- 1 - I_A
  X_obs <- bm_mat
  X_obs[O_mat == 0] <- NA
  
  valid_t <- which(
    colSums(O_mat * I_A) > 0.1 * n &
      colSums(O_mat * I_B) > 0.1 * n
  )
  grid  <- grid_full[valid_t]
  X_obs <- X_obs[, valid_t, drop = FALSE]
  O_mat <- O_mat[, valid_t, drop = FALSE]
  
  pA_hat  <- pmax(colSums(O_mat * I_A) / n, eps_p)
  pB_hat  <- pmax(colSums(O_mat * I_B) / n, eps_p)
  muA_hat <- colSums(X_obs * I_A, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X_obs * I_B, na.rm = TRUE) / (n * pB_hat)
  
  diff <- muA_hat - muB_hat
  
  # L2 mit tf_integrate
  diff_tfd <- tfd(matrix(diff, nrow = 1), arg = grid)
  T_L2 <- n * tf_integrate(diff_tfd^2, arg = grid)
  
  # Supremum
  T_D <- sqrt(n) * max(abs(diff))
  
  boot <- bootstrap_both(X_obs, O_mat, I_A, I_B,
                         muA_hat, muB_hat,
                         grid, n_total = n, Bstar = Bstar, eps = eps_p)
  
  p_vals_L2[sim] <- (sum(boot$L2 >= T_L2) + 1) / (Bstar + 1)
  p_vals_D [sim] <- (sum(boot$D  >= T_D ) + 1) / (Bstar + 1)
  
  setTxtProgressBar(pb, sim)
}
close(pb)

# ------------- Typ-I-Fehler­schätzungen ------------------------
type1_L2 <- mean(p_vals_L2 < alpha)
type1_D  <- mean(p_vals_D  < alpha)

cat("\n-----------------------------------------------------------\n",
    "Typ-I-Fehler (MCAR, n =", n, ", nsim =", nsim, ")\n",
    "  • L2-Test       :", round(type1_L2, 4), "\n",
    "  • Supremums-Test:", round(type1_D,  4), "\n",
    "-----------------------------------------------------------\n")
