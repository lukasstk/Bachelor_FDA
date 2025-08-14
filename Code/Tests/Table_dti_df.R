# --- Einmalige MCAR/MNAR-Auswertung auf dti_df$cca (30% Missing pro INcomplete) ----
suppressPackageStartupMessages({ library(tidyfun); library(tf) })
set.seed(123)

# Daten + Grid
data("dti_df")
cca_all <- dti_df$cca
n        <- 50
grid_len <- 101
grid     <- seq(0, 1, length.out = grid_len)

cca_n <- cca_all[seq_len(min(n, nrow(cca_all))), grid, interpolate = TRUE]
rownames(cca_n) <- paste0("ID", seq_len(nrow(cca_n)))
X_true <- as.matrix(cca_n)
n <- nrow(X_true); m <- ncol(X_true)

# ---- Parameter für Missingness ------------------------------------------------
MISS_FRAC   <- 0.30   # Anteil fehlender Punkte je unvollständiger Kurve
P_COMPLETE  <- 0.50   # Anteil kompletter Kurven

# ---- MCAR: zusammenhängendes Intervall mit genau MISS_FRAC --------------------
algo_mcar_mask <- function(m, miss_frac = 0.30) {
  k <- max(1L, min(m - 1L, round(miss_frac * m)))
  start <- sample.int(m - k + 1L, 1L)
  O <- rep(1L, m); O[start:(start + k - 1L)] <- 0L; O
}

# ---- MNAR: wertsensitiv (Top-k Punkte nach Score) mit exakt MISS_FRAC --------
algo_mnar_mask <- function(x_row, miss_frac = 0.30, beta = 5, jitter_sd = 0.05,
                           block_k = 0L) {
  x <- as.numeric(x_row)
  med <- median(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
  z <- (x - med) / sdv
  
  m <- length(x)
  k_target <- max(1L, min(m - 1L, round(miss_frac * m)))
  
  score <- beta * z + rnorm(m, sd = jitter_sd)  # höherer Score => eher fehlend
  idx   <- order(score, decreasing = TRUE)[1:k_target]
  
  if (block_k > 0L) {
    expanded <- idx
    for (i in idx) expanded <- c(expanded, max(1, i - block_k):min(m, i + block_k))
    expanded <- sort(unique(expanded))
    keep <- order(score[expanded], decreasing = TRUE)[1:k_target]
    idx  <- expanded[keep]
  }
  
  O <- rep(1L, m); O[idx] <- 0L; O
}

# ---- Helfer: X_obs erzeugen ---------------------------------------------------
make_X_obs <- function(X_true, pattern = c("MCAR","MNAR"),
                       p_complete = 0.5, miss_frac = 0.30,
                       beta = 5, block_k = 0L) {
  pattern <- match.arg(pattern)
  n <- nrow(X_true); m <- ncol(X_true)
  I_A <- rbinom(n, 1, p_complete)                   # 1 = komplett, 0 = unvollständig
  if (all(I_A == 1L)) I_A[sample.int(n, 1)] <- 0L
  O <- matrix(1L, n, m)
  for (i in seq_len(n)) if (I_A[i] == 0L) {
    O[i, ] <- if (pattern == "MCAR") {
      algo_mcar_mask(m, miss_frac = miss_frac)
    } else {
      algo_mnar_mask(X_true[i, ], miss_frac = miss_frac, beta = beta, block_k = block_k)
    }
  }
  X_obs <- X_true; X_obs[O == 0L] <- NA_real_
  list(X_obs = X_obs, is_complete = rowSums(is.na(X_obs)) == 0)
}

# Zwei einmalige Datensätze
out_mcar <- make_X_obs(X_true, "MCAR", p_complete = P_COMPLETE, miss_frac = MISS_FRAC)
out_mnar <- make_X_obs(X_true, "MNAR", p_complete = P_COMPLETE, miss_frac = MISS_FRAC, beta = 5, block_k = 0L)

# --- (optional) kurze Kontrolle der realisierten Anteile -----------------------
check_frac <- function(X) c(
  overall = mean(is.na(X)),
  mean_per_curve = mean(rowMeans(is.na(X))),
  mean_incomplete = mean(rowMeans(is.na(X))[rowSums(is.na(X)) > 0], na.rm = TRUE)
)
print(rbind(MCAR = check_frac(out_mcar$X_obs),
            MNAR = check_frac(out_mnar$X_obs)))

# --- Tests je Pattern ----------------------------------------------------------
B_asym <- 2000; B_boot <- 2000
res_mcar_L2 <- tfu_algo1_L2_test(X_obs = out_mcar$X_obs, B = B_asym, seed = 1001)
res_mcar_D  <- tfu_algo2_sup_test(X_obs = out_mcar$X_obs, B = B_asym, seed = 1002)
res_mcar_BT <- tfu_algo5_bootstrap(X_obs = out_mcar$X_obs, B = B_boot, seed = 1003)
res_mnar_L2 <- tfu_algo1_L2_test(X_obs = out_mnar$X_obs, B = B_asym, seed = 2001)
res_mnar_D  <- tfu_algo2_sup_test(X_obs = out_mnar$X_obs, B = B_asym, seed = 2002)
res_mnar_BT <- tfu_algo5_bootstrap(X_obs = out_mnar$X_obs, B = B_boot, seed = 2003)

results <- data.frame(
  pattern   = c("MCAR", "MNAR"),
  asym_L2_p = c(res_mcar_L2$p_value, res_mnar_L2$p_value),
  asym_D_p  = c(res_mcar_D$p_value,  res_mnar_D$p_value),
  boot_L2_p = c(res_mcar_BT$p_L2,    res_mnar_BT$p_L2),
  boot_D_p  = c(res_mcar_BT$p_D,     res_mnar_BT$p_D)
)
print(results, row.names = FALSE)

# ------------------------- Visualisierung --------------------------------------
y_range <- range(c(out_mcar$X_obs, out_mnar$X_obs), na.rm = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot_panel <- function(title, X_obs, is_complete) {
  Xc <- X_obs[ is_complete, , drop = FALSE]
  Xi <- X_obs[!is_complete, , drop = FALSE]
  
  if (nrow(Xc) > 0) {
    matplot(grid, t(Xc), type = "l", lty = 1, lwd = 1,
            col = "grey70", xlab = "t", ylab = "y(t)", ylim = y_range,
            main = title)
    if (nrow(Xi) > 0) matlines(grid, t(Xi), lty = 1, lwd = 1.2, col = "black")
  } else {
    matplot(grid, t(Xi), type = "l", lty = 1, lwd = 1.2,
            col = "black", xlab = "t", ylab = "y(t)", ylim = y_range,
            main = title)
  }
  
  legend("topleft",
         legend = c(sprintf("complete (n=%d)", nrow(Xc)),
                    sprintf("incomplete (n=%d)", nrow(Xi))),
         lty = 1, lwd = c(1, 1.2), bty = "n",
         col = c("grey70", "black"))
}

plot_panel("MCAR: complete=grau, incomplete=schwarz", out_mcar$X_obs, out_mcar$is_complete)
plot_panel("MNAR: complete=grau, incomplete=schwarz", out_mnar$X_obs, out_mnar$is_complete)

par(mfrow = c(1, 1))

# ================================================================
# Type-I-Error @ n = 50 für MCAR & MNAR (30% Missing je INcomplete)
# Voraussetzungen:
# - dti_df geladen (tidyfun, tf)
# - algo_mcar_mask(), algo_mnar_mask(), make_X_obs() definiert (wie oben)
# - tfu_algo1_L2_test(), tfu_algo2_sup_test(), tfu_algo5_bootstrap() verfügbar
# ================================================================

suppressPackageStartupMessages({ library(tidyfun); library(tf) })

# Zufällige Stichprobe von n Kurven aus dti_df$cca auf gleichmäßigem Grid
.get_X_true <- function(n = 50, grid_len = 101) {
  # dti_df laden
  data("dti_df", envir = environment())
  cca_all <- dti_df$cca
  
  # Anzahl verfügbarer Kurven robust bestimmen
  n_all <- length(cca_all)              # <-- statt nrow(cca_all)
  if (is.null(n_all) || n_all < 1) stop("Keine Kurven in dti_df$cca gefunden.")
  
  # Zufallsstichprobe von max. n Kurven
  idx  <- sample.int(n_all, size = min(n, n_all), replace = FALSE)
  
  # gleichmäßiges Grid und Interpolation
  grid <- seq(0, 1, length.out = grid_len)
  cca_n <- cca_all[idx, grid, interpolate = TRUE]
  
  # Matrix der Funktionswerte zurückgeben
  rownames(cca_n) <- paste0("ID", seq_len(length(idx)))
  list(X_true = as.matrix(cca_n), grid = grid)
}

# Ein Replikat: p<alpha-Indikatoren der vier Tests
.one_rep <- function(pattern = c("MCAR","MNAR"),
                     alpha = 0.05, B_asym = 2000, B_boot = 2000,
                     p_complete = 0.5, miss_frac = 0.30,
                     beta = 5, block_k = 0L, grid_len = 101, seed = 1) {
  pattern <- match.arg(pattern)
  set.seed(seed)
  sam <- .get_X_true(n = 50, grid_len = grid_len)
  X_true <- sam$X_true
  
  out <- if (pattern == "MCAR") {
    make_X_obs(X_true, "MCAR", p_complete = p_complete, miss_frac = miss_frac)
  } else {
    make_X_obs(X_true, "MNAR", p_complete = p_complete, miss_frac = miss_frac,
               beta = beta, block_k = block_k)
  }
  
  # Tests: Algo 1 (L2), Algo 2 (Sup), Algo 5 (Bootstrap L2 & D)
  res_L2 <- tfu_algo1_L2_test(X_obs = out$X_obs, B = B_asym, seed = seed + 11)
  res_D  <- tfu_algo2_sup_test(X_obs = out$X_obs, B = B_asym, seed = seed + 12)
  res_BT <- tfu_algo5_bootstrap(X_obs = out$X_obs, B = B_boot, seed = seed + 13)
  
  c(
    asym_L2 = as.numeric(res_L2$p_value < alpha),
    asym_D  = as.numeric(res_D$p_value  < alpha),
    boot_L2 = as.numeric(res_BT$p_L2    < alpha),
    boot_D  = as.numeric(res_BT$p_D     < alpha)
  )
}

# Type-I-Error schätzen (Anteil Zurückweisungen unter H0)
run_type1 <- function(pattern = c("MCAR","MNAR"), n = 50,
                          R = 1000, alpha = 0.05,
                          B_asym = 2000, B_boot = 2000,
                          p_complete = 0.5, miss_frac = 0.30,
                          beta = 5, block_k = 0L, grid_len = 101,
                          base_seed = 42, quiet = TRUE) {
  pattern <- match.arg(pattern)
  if (quiet) {
    old <- getOption("warn"); options(warn = 1); on.exit(options(warn = old), add = TRUE)
  }
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  M <- matrix(0, nrow = R, ncol = 4)
  colnames(M) <- c("asym_L2","asym_D","boot_L2","boot_D")
  for (r in seq_len(R)) {
    M[r, ] <- .one_rep(
      pattern = pattern, alpha = alpha, B_asym = B_asym, B_boot = B_boot,
      p_complete = p_complete, miss_frac = miss_frac,
      beta = beta, block_k = block_k, grid_len = grid_len,
      seed = base_seed + 1000 * r
    )
    setTxtProgressBar(pb, r)
  }
  close(pb)
  est <- colMeans(M)
  data.frame(
    pattern = pattern, n = n,
    asym_L2 = round(est["asym_L2"], 3),
    asym_D  = round(est["asym_D"],  3),
    boot_L2 = round(est["boot_L2"], 3),
    boot_D  = round(est["boot_D"],  3),
    row.names = NULL
  )
}

# ---------------------- AUSFÜHREN (n = 50) ------------------------------------
# Hinweis: R=1000 und B=2000/2000 sind rechenintensiv. Für einen Probelauf:
# R=200 oder B_asym/B_boot kleiner wählen und dann hochdrehen.

set.seed(2025)
res_mcar_50 <- run_type1("MCAR", n = 50, R = 100, B_asym = 2000, B_boot = 2000,
                          p_complete = 0.5, miss_frac = 0.30)

res_mcar_100 <- run_type1("MCAR", n = 100, R = 100, B_asym = 2000, B_boot = 2000,
                          p_complete = 0.5, miss_frac = 0.30)


type1_table_n50 <- rbind(res_mcar_50, res_mcar_100)
type1_table_n50
type1_table_n50 <- rbind(res_mcar_50, res_mcar_100)
# > type1_table_n50
# pattern   n asym_L2 asym_D boot_L2 boot_D
# 1    MCAR  50    0.03   0.03    0.03   0.02
# 2    MCAR 100    0.03   0.03    0.03   0.02
