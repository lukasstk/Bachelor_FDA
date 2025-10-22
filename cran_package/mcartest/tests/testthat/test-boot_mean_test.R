# ==========================================================
# Helper functions (identical to asymptotic tests)
# ==========================================================
simulate_bm <- function(grid, mean_shift = 0) {
  dt <- diff(grid)[1]
  c(0, cumsum(stats::rnorm(length(grid) - 1, sd = sqrt(dt)))) + mean_shift
}

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

make_O_mnar <- function(x_row) {
  if (all(x_row > -1 & x_row < 2)) {
    rep(1, length(x_row))
  } else {
    as.numeric(x_row > -1 & x_row < 2)
  }
}

suppress_bootstrap_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Only [0-9]+ valid bootstrap replicates",
                conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# ==========================================================
# Deterministic behavior when seed is set (L2 and D tests, bootstrap)
# ==========================================================
test_that("bootstrap tests produce deterministic output when seed is set", {
  set.seed(123)

  grid <- seq(0, 1, length.out = 100)
  n <- 100

  # Simulate Brownian motion sample
  X <- t(replicate(n, simulate_bm(grid)))

  # Apply random MCAR censoring
  O <- t(replicate(n, make_O_mcar(grid)))
  X[O == 0L] <- NA_real_

  # Run both tests (L2 and D) together
  suppress_bootstrap_warning({
      res_1 <- boot_mean_test(
      X = X,
      n_boot = 2500,
      stat = c("L2", "D"),
      compute_bands = FALSE,
      manage_backend = "auto",
      ncpus = 2,
      seed = 99
      )

    res_2 <- boot_mean_test(
      X = X,
      n_boot = 2500,
      stat = c("L2", "D"),
      compute_bands = FALSE,
      manage_backend = "auto",
      ncpus = 2,
      seed = 99
      )
  })

  # Expect two outputs (one per stat)
  expect_length(res_1, 2)
  expect_named(res_1, c("L2", "D"))

  # Check reproducibility of p-values
  expect_equal(res_1$L2$p.value, res_2$L2$p.value, tolerance = 1e-12)
  expect_equal(res_1$D$p.value,  res_2$D$p.value,  tolerance = 1e-12)
})


# ==========================================================
# Correct structure when bands_only = TRUE (Supremum test, bootstrap)
# ==========================================================
test_that("boot_mean_test returns correct structure when bands_only = TRUE", {
  set.seed(123)
  grid <- seq(0, 1, length.out = 100)
  n <- 100

  # Simulate Brownian motion trajectories
  X <- t(replicate(n, simulate_bm(grid)))

  # Apply MCAR censoring
  O <- t(replicate(n, make_O_mcar(grid)))
  X[O == 0L] <- NA_real_

  # Run Bootstrap Supremum test with compute_bands = TRUE and bands_only = TRUE
  suppress_bootstrap_warning({
    res_boot <- boot_mean_test(
      X = X,
      n_boot = 2500,
      stat = "D",
      compute_bands = TRUE,
      chunk_size = 50,
      bands_only = TRUE,
      manage_backend = "auto",
      ncpus = 2,
      seed = 42
    )
  })

  # Structural checks
  expect_type(res_boot, "list")
  expect_named(res_boot, c("estimate", "parameter", "bands"))
  expect_s3_class(res_boot$estimate$mean_A, "tfd")
  expect_s3_class(res_boot$estimate$mean_B, "tfd")
  expect_s3_class(res_boot$estimate$mean_diff, "tfd")
  expect_true("n_boot" %in% names(res_boot$parameter))
  expect_true(all(c("lower", "upper", "band", "alpha", "grid")
                  %in% names(res_boot$bands)))
})


# ==========================================================
# Type-I Error under MCAR — Bootstrap test
# ==========================================================
test_that("bootstrap tests maintain nominal Type-I error under MCAR", {
  set.seed(123)

  grid <- seq(0, 1, length.out = 100)
  n <- 250
  alpha <- 0.05
  runs <- 200
  n_boot <- 1000

  rejections_L2 <- logical(runs)
  rejections_D  <- logical(runs)
  suppress_bootstrap_warning({
    for (r in seq_len(runs)) {
      X <- t(replicate(n, simulate_bm(grid)))
      O <- t(replicate(n, make_O_mcar(grid)))
      X[O == 0L] <- NA_real_

      res_L2 <- boot_mean_test(X = X, n_boot = n_boot, stat = "L2",
                               manage_backend = "auto", ncpus = 2)
      res_D  <- boot_mean_test(X = X, n_boot = n_boot, stat = "D",
                               manage_backend = "auto", ncpus = 2)

      rejections_L2[r] <- res_L2$p.value < alpha
      rejections_D[r]  <- res_D$p.value  < alpha
    }
  })

  rate_L2 <- mean(rejections_L2)
  rate_D  <- mean(rejections_D)

  expect_lt(rate_L2, 0.091)
  expect_lt(rate_D,  0.091)
})

# ==========================================================
# MNAR Detection — Bootstrap test
# ==========================================================
test_that("bootstrap tests detect systematic MNAR bias", {
  set.seed(42)

  grid <- seq(0, 1, length.out = 100)
  n <- 250
  alpha <- 0.05
  runs <- 200
  n_boot <- 1000

  rejections_L2 <- logical(runs)
  rejections_D  <- logical(runs)

  suppress_bootstrap_warning({
    for (r in seq_len(runs)) {
      X <- t(replicate(n, simulate_bm(grid)))
      O <- t(apply(X, 1L, make_O_mnar))
      X[O == 0L] <- NA_real_

      res_L2 <- boot_mean_test(X = X, n_boot = n_boot, stat = "L2",
                               manage_backend = "auto", ncpus = 2)
      res_D  <- boot_mean_test(X = X, n_boot = n_boot, stat = "D",
                               manage_backend = "auto", ncpus = 2)

      rejections_L2[r] <- res_L2$p.value < alpha
      rejections_D[r]  <- res_D$p.value  < alpha
    }
  })

  rate_L2 <- mean(rejections_L2)
  rate_D  <- mean(rejections_D)

  expect_gt(rate_L2, 0.8)
  expect_gt(rate_D,  0.8)
})

# ==========================================================
# MNAR Sensitivity / Power Curve — Bootstrap test
# ==========================================================
test_that("bootstrap tests show increasing rejection
          with stronger MNAR (larger b)", {
  set.seed(2025)

  n <- 100
  grid <- seq(0, 1, length.out = 100)
  alpha <- 0.05
  b_vals <- seq(1.0, 2, by = 0.2)
  n_sims <- 200
  n_boot <- 1000

  make_missing_pattern <- function(x_row, a = -1, b = 2) {
    as.integer(x_row > a & x_row < b)
  }

  one_run <- function(b_now) {
    bm_mat <- t(replicate(n, simulate_bm(grid)))
    O <- t(apply(bm_mat, 1L, make_missing_pattern, a = -1, b = b_now))
    bm_mat[O == 0L] <- NA_real_

    p_L2 <- tryCatch(
      boot_mean_test(X = bm_mat, n_boot = n_boot, stat = "L2",
                     manage_backend = "auto", ncpus = 2)$p.value,
      error = function(e) NA_real_
    )
    p_D <- tryCatch(
      boot_mean_test(X = bm_mat, n_boot = n_boot, stat = "D",
                     manage_backend = "auto", ncpus = 2)$p.value,
      error = function(e) NA_real_
    )

    c(p_L2, p_D)
  }

  suppress_bootstrap_warning({
    res_list <- lapply(b_vals, function(b_now) {
      ps <- replicate(n_sims, one_run(b_now))
      rej <- rowMeans(ps < alpha, na.rm = TRUE)
      data.frame(b = b_now, L2 = rej[1], D = rej[2])
    })
  })

  res_df <- do.call(rbind, res_list)

  # Expect higher rejection rate (stronger MNAR → lower b → more rejections)
  expect_true(all(diff(res_df$L2) > 0))
  expect_true(all(diff(res_df$D)  > 0))
})
