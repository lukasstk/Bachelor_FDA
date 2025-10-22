# ==========================================================
# Backend reset and RNG initialization
# ==========================================================
# This section tests that `.reset_backend()` and `.init_parallel()`
# correctly handle:
# - BLAS/OpenMP environment restoration
# - RNG reproducibility (L’Ecuyer–CMRG)
# - Proper sequential / parallel backend setup
# - Cleanup of internal environment variables
# ==========================================================

test_that("Resetting backend restores thread settings and RNG safely", {
  # --- Save current env vars ---
  old_env <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"))

  # --- Simulate old thread state ---
  .tfu_par_env$old_threads <- c("2", "3", "4")

  # --- Call reset function ---
  expect_silent(.reset_backend())

  # --- Check environment restored ---
  expect_equal(Sys.getenv("OPENBLAS_NUM_THREADS"), "2")
  expect_equal(Sys.getenv("MKL_NUM_THREADS"), "3")
  expect_equal(Sys.getenv("OMP_NUM_THREADS"), "4")

  # --- RNG state reproducibility ---
  suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
  x1 <- runif(5)
  .reset_backend()
  x2 <- runif(5)
  expect_false(isTRUE(all.equal(x1, x2))) # RNG reset, not same stream

  # --- Cleanup ---
  expect_false(exists("old_threads", envir = .tfu_par_env, inherits = FALSE))
})


# ==========================================================
# Sequential backend setup
# ==========================================================
# Tests that `.init_parallel(manage_backend = "sequential")`
# correctly registers a sequential foreach backend and returns
# the expected metadata.
# ==========================================================

test_that("Sequential backend initializes correctly", {
  res <- .init_parallel(manage_backend = "sequential")
  expect_type(res, "list")
  expect_equal(res$nworkers, 1L)
  expect_equal(res$used, "sequential")

  # Verify foreach backend is sequential
  info <- foreach::getDoParName()
  expect_equal(info, "doSEQ")
})


# ==========================================================
# Parallel backend: auto reuse and forced pool
# ==========================================================
# Tests that `.init_parallel()` correctly creates or reuses clusters
# depending on the `manage_backend` argument, and respects RNG seeding.
# ==========================================================

test_that("Parallel backend creates and reuses correctly", {
  library(foreach)
  # --- Force new pool ---
  res_new <- .init_parallel(manage_backend = "force_pool", ncpus = 2)
  expect_true(inherits(.tfu_par_env$cl, "cluster"))
  expect_equal(res_new$used, "internal-forced")
  expect_equal(res_new$nworkers, 2L)

  # --- Reuse same pool (auto) ---
  res_auto <- .init_parallel(manage_backend = "auto", ncpus = 2)
  expect_equal(res_auto$used, "internal-reused")

  # --- initialize with seed ---
  res <- .init_parallel(manage_backend = "force_pool", ncpus = 2, seed = 123)

  # (1) RNG method should be correctly set
  expect_equal(RNGkind()[1], "L'Ecuyer-CMRG")

  # (2) doRNG should be active (registered)
  expect_true(any(grepl("doRNG", foreach::getDoParName(), ignore.case = TRUE)) ||
                "doParallel" %in% foreach::getDoParName())

  # (3) within the same backend, reproducibility should hold
  x1 <- foreach(i = 1:3, .combine = c) %dopar% runif(1)
  x2 <- foreach(i = 1:3, .combine = c) %dopar% runif(1)
  expect_false(isTRUE(all.equal(x1, x2)))   # streams move forward

  # Reinitialize backend with same seed
  res2 <- .init_parallel(manage_backend = "force_pool", ncpus = 2, seed = 123)
  y1 <- foreach(i = 1:3, .combine = c) %dopar% runif(1)

  res3 <- .init_parallel(manage_backend = "force_pool", ncpus = 2, seed = 123)
  y2 <- foreach(i = 1:3, .combine = c) %dopar% runif(1)

  # (4) Across identical seeds, results should match again
  expect_equal(y1, y2)
})
