# ==========================================================
# .trapezoid_weights
# ==========================================================
test_that("test functionality of .trapezoid_weights", {
  grid <- seq(0, 1, length.out = 5)
  w <- .trapezoid_weights(grid)
  expect_length(w, 5)
  expect_true(all(w > 0))
  expect_equal(sum(w), 1, tolerance = 1e-2)

  expect_error(.trapezoid_weights(rev(grid)))       # sorted = TRUE violated
  expect_error(.trapezoid_weights(c(0, NA, 1)))     # NA in grid
  expect_error(.trapezoid_weights(numeric(0)))      # min.len violated

  # --- Edge case: single grid point ---
  single_grid <- 0.5
  w_single <- .trapezoid_weights(single_grid)
  expect_equal(w_single, 1)         # should return exactly 1
})


# ==========================================================
# .group_mean_estimators
# ==========================================================
test_that("test functionality of .group_mean_estimators", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow = 5)
  O <- matrix(1L, nrow = 5, ncol = 4)
  group_A <- c(TRUE, TRUE, FALSE, FALSE, TRUE)

  res <- .group_mean_estimators(X, O, group_A)
  expect_named(res, c("mean_A", "mean_B", "pA", "pB"))
  expect_length(res$mean_A, 4)

  expected_mean_A <- c(0.02264198, -0.32700300, -0.02730960, 1.06305755)
  expected_mean_B <- c(-0.1630356, -1.2529771, -0.7037129, 0.2781076)
  expected_pA <- rep(0.6, 4)
  expected_pB <- rep(0.4, 4)

  expect_equal(res$mean_A, expected_mean_A, tolerance = 1e-6)
  expect_equal(res$mean_B, expected_mean_B, tolerance = 1e-6)
  expect_equal(res$pA, expected_pA)
  expect_equal(res$pB, expected_pB)
})


# ==========================================================
# .covariance_estimator
# ==========================================================
test_that("test functionality of .covariance_estimator", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow = 5)
  O <- matrix(1L, nrow = 5, ncol = 4)
  group_A <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  gm <- .group_mean_estimators(X, O, group_A)

  K <- .covariance_estimator(X, O, group_A, gm$mean_A, gm$mean_B, gm$pA, gm$pB)
  expect_true(is.matrix(K))
  expect_equal(dim(K), c(4, 4))
  expect_equal(K, t(K))          # symmetric
  expect_false(any(is.na(K)))    # valid numerically

  eigvals <- eigen(K, symmetric = TRUE)$values
  expect_true(all(eigvals >= -1e-8))  # PSD check
})


# ==========================================================
# .kl_decomposition
# ==========================================================
test_that("test functionality of .kl_decomposition", {
  grid <- seq(0, 1, length.out = 5)
  K <- diag(5)
  res <- .kl_decomposition(K, grid)

  expect_named(res, c("eigenvalues", "eigenfunctions", "w"))
  expect_true(all(res$eigenvalues >= 0))
  expect_equal(sum(res$w), 1, tolerance = 1e-2)

  K_bad <- diag(c(1, -1, 1, 1, 1))
  expect_error(.kl_decomposition(K_bad, grid))
})


# ==========================================================
# .confidence_bands
# ==========================================================
test_that("test functionality of .confidence_bands", {
  skip_if_not_installed("tf")

  n <- 10
  grid <- seq(0, 1, length.out = 5)
  diff <- seq(-1, 1, length.out = 5)
  W <- rnorm(100, mean = 0.5)

  res <- .confidence_bands("D", diff, W, n, 0.05, grid)
  expect_named(res, c("type", "band", "lower", "upper", "alpha", "grid"))
  expect_true(all(res$upper >= res$lower))
  expect_equal(res$type, "simultaneous")
  expect_equal(res$alpha, 0.05)
  expect_s3_class(res$band, "tfd")

  expect_error(.confidence_bands("L2", diff, W, n, 0.05, grid))
})


# ==========================================================
# .create_output
# ==========================================================
test_that("test functionality of .create_output", {
  stat_value <- 3.14
  p_value <- 0.05
  res <- .create_output(
    stat_name = "Tµ,L²",
    stat_value = stat_value,
    p_value = p_value,
    method = "Test",
    data_name = "X"
  )
  expect_s3_class(res, "htest")
  expect_equal(res$method, "Test")
  expect_equal(res$data.name, "X")
  expect_true(is.numeric(res$statistic))

  res_na <- .create_output(
    stat_name = "Tµ,L²",
    stat_value = stat_value,
    p_value = NA,
    method = "Test",
    data_name = "X"
  )
  expect_true(is.na(res_na$p.value))
})
