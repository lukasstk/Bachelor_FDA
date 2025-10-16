# ==========================================================
# .groups_to_logical
# ==========================================================
test_that("test functionality of .groups_to_logical", {
  # --- Valid cases ---
  g <- c("A", "B", "A", "B")
  res <- .groups_to_logical(g)
  expect_type(res, "logical")
  expect_length(res, 4)
  expect_true(all(res %in% c(TRUE, FALSE)))

  g_fac <- factor(c("low", "high", "low", "high"))
  res_fac <- .groups_to_logical(g_fac)
  expect_type(res_fac, "logical")

  g_num <- c(1, 2, 1, 2)
  res_num <- .groups_to_logical(g_num)
  expect_type(res_num, "logical")
  expect_equal(sum(res_num), 2)

  g_log <- c(TRUE, FALSE, TRUE)
  res_log <- .groups_to_logical(g_log)
  expect_identical(res_log, g_log)

  # --- Error cases ---
  expect_error(.groups_to_logical(c("A", "B", "C")))     # too many groups
  expect_error(.groups_to_logical(c("A", NA, "B")))      # missing values
  expect_error(.groups_to_logical(list("A", "B", "A")))  # non-atomic
})


# ==========================================================
# .fd_to_matrix
# ==========================================================
test_that("test functionality of .fd_to_matrix", {
  # Regular
  grid <- seq(0, 1, length.out = 5)
  mat <- matrix(1:10, nrow = 2)
  fd <- tf::tfd(mat, arg = grid)
  res <- .fd_to_matrix(fd)
  expect_named(res, c("X", "grid"))
  expect_equal(res$grid, grid)

  # Irregular
  grid_list <- list(seq(0, 0.8, length.out = 4), seq(0, 1, length.out = 6))
  mat_list <- list(rnorm(4), rnorm(6))
  fd_irreg <- tf::tfd(mat_list, arg = grid_list)
  res_irreg <- .fd_to_matrix(fd_irreg)
  expect_true(all(sort(unique(unlist(grid_list))) %in% res_irreg$grid))

  # Error case
  expect_error(.fd_to_matrix(matrix(1:4, 2)))
})


# ==========================================================
# .limit_subdomain
# ==========================================================
test_that("test functionality of .limit_subdomain", {
  O <- matrix(sample(0:1, 20, replace = TRUE), nrow = 4)
  group_A <- c(TRUE, TRUE, FALSE, FALSE)
  res <- .limit_subdomain(O, group_A, min_frac = 0.10)
  expect_named(res, c("idx", "min_frac_used"))
  expect_true(all(res$idx %in% seq_len(ncol(O))))
  expect_equal(res$min_frac_used, 0.10)

  # Edge/error cases
  expect_error(.limit_subdomain(matrix(0, nrow = 4, ncol = 5), group_A))
  expect_error(.limit_subdomain(matrix(1:4, 2), c(TRUE, FALSE, TRUE)))
  expect_error(.limit_subdomain(matrix(1:4, 2), c(TRUE, FALSE), min_frac = -1))
  expect_error(.limit_subdomain(matrix(1:4, 2), c(TRUE, FALSE), min_frac = 1.5))
  expect_error(.limit_subdomain(matrix(c(1, NA, 0, 1), nrow = 2),
                                c(TRUE, FALSE)))
  expect_error(.limit_subdomain(matrix(c("a", "b", "c", "d"), nrow = 2),
                                c(TRUE, FALSE)))
})


# ==========================================================
# .prepare_inputs
# ==========================================================
test_that("test functionality of .prepare_inputs", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow = 5)
  X[sample(length(X), 5)] <- NA

  # Valid automatic grouping
  res <- .prepare_inputs(X = X, observed_ratio = 0.51)
  expect_named(res, c("X", "O", "group_A", "grid"))
  expect_true(is.matrix(res$X))
  expect_true(all(res$O %in% c(0, 1)))

  # Fails when one group empty (threshold too low)
  expect_error(.prepare_inputs(X = X, observed_ratio = 0.5))

  # Manual grouping
  X2 <- matrix(rnorm(12), nrow = 4)
  X2[sample(length(X2), 2)] <- NA
  groups <- c("A", "B", "A", "B")
  res2 <- .prepare_inputs(X = X2, groups = groups)
  expect_length(res2$group_A, 4)

  # Edge cases
  expect_error(.prepare_inputs(X = matrix(numeric(0), nrow = 0, ncol = 0)))
  expect_error(.prepare_inputs(X = X2, observed_ratio = 2))
  expect_error(.prepare_inputs(X = X2, groups = c(TRUE, FALSE)))
})
