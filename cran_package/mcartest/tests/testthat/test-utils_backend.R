# ==========================================================
# Backend reset and RNG initialization
# ==========================================================

test_that(".reset_backend restores thread settings and RNG safely", {
  # --- Save current env vars ---
  old_env <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", 
                          "OMP_NUM_THREADS"))

  # --- Simulate old thread state ---
  .tfu_par_env$old_threads <- c("2", "3", "4")

  # --- Call reset function ---
  expect_silent(.reset_backend())

  # --- Check environment restored ---
  expect_equal(Sys.getenv("OPENBLAS_NUM_THREADS"), "2")
  expect_equal(Sys.getenv("MKL_NUM_THREADS"), "3")
  expect_equal(Sys.getenv("OMP_NUM_THREADS"), "4")

  # --- RNG kind should be reset to L'Ecuyer-CMRG ---
  expect_equal(RNGkind()[1], "L'Ecuyer-CMRG")

  # --- Internal cleanup: old_threads should be removed ---
  expect_false(exists("old_threads", envir = .tfu_par_env, inherits = FALSE))

})


# ==========================================================
# Sequential/auto/force backend initialization
# ==========================================================

test_that(".init_parallel initializes sequential backend correctly", {
  res_seq <- .init_parallel(manage_backend = "sequential")
  expect_type(res_seq, "list")
  expect_equal(res_seq$nworkers, 1L)
  expect_equal(res_seq$used, "sequential")
  expect_equal(foreach::getDoParName(), "doSEQ")

  # --- Force new pool ---
  res_new <- .init_parallel(manage_backend = "force_pool", ncpus = 2)
  expect_true(inherits(.tfu_par_env$cl, "cluster"))
  expect_equal(res_new$used, "internal-forced")
  expect_equal(res_new$nworkers, 2L)
  expect_type(res_new, "list")
  expect_equal(foreach::getDoParName(), "doParallelSNOW")

  # --- Reuse same pool (auto) ---
  res_auto <- .init_parallel(manage_backend = "auto", ncpus = 2)
  expect_equal(res_auto$used, "internal-reused")
  expect_equal(res_auto$nworkers, 2L)
  expect_type(res_auto, "list")
  expect_equal(foreach::getDoParName(), "doParallelSNOW")
  
  # --- Force new pool with seed ---
  res_new_seed <- .init_parallel(manage_backend = "force_pool", 
                                 ncpus = 2, seed = 1)
  expect_true(inherits(.tfu_par_env$cl, "cluster"))
  expect_equal(res_new_seed$used, "internal-forced")
  expect_equal(res_new_seed$nworkers, 2L)
  expect_type(res_new_seed, "list")
  expect_equal(foreach::getDoParName(), "doRNG")
  
  # --- Reuse same pool (auto) with seed ---
  res_auto_seed <- .init_parallel(manage_backend = "auto", 
                                  ncpus = 2,seed = 1)
  expect_equal(res_auto_seed$used, "internal-reused")
  expect_equal(res_auto_seed$nworkers, 2L)
  expect_type(res_auto_seed, "list")
  expect_equal(foreach::getDoParName(), "doRNG")
})


# ==========================================================
# Parallel backend reuse and RNG seeding
# ==========================================================

test_that(".init_parallel initializes correctly and ensures reproducibility", {
  library(foreach)

  # --- initialize with seed ---
  res <- .init_parallel(manage_backend = "force_pool", ncpus = 2, seed = 123)

  # within the same backend, reproducibility should hold
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
