# ==========================================================
# Helper functions
# ==========================================================
simulate_bm <- function(grid, mean_shift = 0) {
  dt <- diff(grid)[1]
  c(0, cumsum(rnorm(length(grid) - 1, sd = sqrt(dt)))) + mean_shift
}

# MCAR: zufällige Intervalle, immer >50% Coverage möglich
make_O_mcar <- function(grid) {
  if (runif(1) < 0.5) {
    rep(1L, length(grid))
  } else {
    U1 <- runif(1); U2 <- runif(1)
    L  <- min(U1, U2)
    U  <- max(U1, U2)
    as.integer(grid >= L & grid < U)
  }
}

# MNAR: Beobachtung hängt vom Wert selbst ab
make_O_mnar <- function(x_row) {
  if (all(x_row > -1 & x_row < 2)) {
    rep(1, length(x_row))
  } else {
    as.numeric(x_row > -1 & x_row < 2)
  }
}


# ==========================================================
# Type-I Error Simulations for MCAR and MNAR
# (asymptotic tests: L2 and Supremum)
# ==========================================================

test_that("type-I error stays near nominal level under MCAR", {
  set.seed(123)

  grid <- seq(0, 1, length.out = 100)
  n <- 250
  alpha <- 0.05
  runs <- 1000
  n_sim <- 2500

  rejections_L2 <- logical(runs)
  rejections_D  <- logical(runs)

  for (r in seq_len(runs)) {
    # Simulate Brownian motion data under H0
    X <- t(replicate(n, simulate_bm(grid)))
    # Random MCAR censoring with 50% complete curves
    O <- t(replicate(n, make_O_mcar(grid)))
    X[O == 0L] <- NA_real_

    res_L2 <- asym_mean_L2_test(X = X, n_sim = n_sim)
    res_D  <- asym_mean_sup_test(X = X, n_sim = n_sim, compute_bands = FALSE)

    rejections_L2[r] <- res_L2$p.value < alpha
    rejections_D[r]  <- res_D$p.value  < alpha
  }

  rate_L2 <- mean(rejections_L2)
  rate_D  <- mean(rejections_D)

  # Should stay close to nominal 0.05
  expect_lt(rate_L2, 0.07)
  expect_lt(rate_D,  0.07)
})

# ==========================================================
# MNAR Detection Test (Robustness Check)
# ----------------------------------------------------------
# Goal: verify that asymptotic tests systematically reject H0
# under MNAR missingness (p-values < 0.05 most of the time)
# ==========================================================

test_that("asymptotic tests reliably detect MNAR bias (p < 0.05)", {
  set.seed(123)

  grid <- seq(0, 1, length.out = 100)
  n <- 250
  alpha <- 0.05
  runs <- 1000
  n_sim <- 2500

  rejections_L2 <- logical(runs)
  rejections_D  <- logical(runs)

  for (r in seq_len(runs)) {
    # 1️⃣ Simulate Brownian motion trajectories
    X <- t(replicate(n, simulate_bm(grid)))

    # 2️⃣ Apply MNAR mechanism (-1 < X < 2 → observed)
    O <- t(apply(X, 1L, make_O_mnar))
    X[O == 0L] <- NA_real_

    # 3️⃣ Run asymptotic tests
    res_L2 <- asym_mean_L2_test(X = X, n_sim = n_sim)
    res_D  <- asym_mean_sup_test(X = X, n_sim = n_sim, compute_bands = FALSE)

    # 4️⃣ Record detection (reject if p < alpha)
    rejections_L2[r] <- res_L2$p.value < alpha
    rejections_D[r]  <- res_D$p.value  < alpha
  }

  rate_L2 <- mean(rejections_L2)
  rate_D  <- mean(rejections_D)

  # Expect high detection (> 0.8 means the test detects MNAR in 80%+ of cases)
  expect_gt(rate_L2, 0.8)
  expect_gt(rate_D,  0.8)
})

# ==========================================================
# Power / MNAR Sensitivity: rejection probability
# increases with stronger MNAR (smaller b)
# ==========================================================

test_that("asymptotic tests detect stronger MNAR censoring
          (higher rejection with smaller b)", {
  set.seed(42)

  n        <- 100                   # sample size per run
  grid     <- seq(0, 1, length.out = 100)
  alpha    <- 0.05
  b_vals   <- seq(1.0, 2, by = 0.2)   # vary censoring strength (a=-1 fixed)
  n_sims   <- 1000                   # Monte Carlo repetitions
  min_frac <- 0.10

  make_missing_pattern <- function(x_row, a = -1, b = 2) {
    as.integer(x_row > a & x_row < b)
  }

  one_run <- function(b_now) {
    # Simulate n Brownian paths
    bm_mat <- t(replicate(n, simulate_bm(grid)))

    # MNAR censoring: observed only if -1 < X(t) < b
    O <- t(apply(bm_mat, 1L, make_missing_pattern, a = -1, b = b_now))
    bm_mat[O == 0L] <- NA_real_

    # Auto-grouping handles observed_ratio = 1 internally
    p_L2 <- tryCatch(
      asym_mean_L2_test(X = bm_mat, n_sim = 500, min_frac = min_frac)$p.value,
      error = function(e) NA_real_
    )
    p_D <- tryCatch(
      asym_mean_sup_test(X = bm_mat, n_sim = 500, compute_bands = FALSE,
                         min_frac = min_frac)$p.value,
      error = function(e) NA_real_
    )

    c(p_L2, p_D)
  }

  res_list <- lapply(b_vals, function(b_now) {
    ps <- replicate(n_sims, one_run(b_now))
    rej <- rowMeans(ps < alpha, na.rm = TRUE)
    data.frame(b = b_now, L2 = rej[1], D = rej[2])
  })

  res_df <- do.call(rbind, res_list)

  # Expect higher rejection rate
  expect_true(all(diff(res_df$L2) > 0))
  expect_true(all(diff(res_df$D)  > 0))
})
