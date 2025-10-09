#' Fully reset foreach backend and RNG (use after user interrupt)
#' @keywords internal
#' @noRd
.reset_backend <- function() {
  # Reset foreach/doRNG
  try(foreach::registerDoSEQ(), silent = TRUE)
  try(doRNG::registerDoRNG(NULL), silent = TRUE)

  # Stop implicit cluster
  suppressWarnings(try(doParallel::stopImplicitCluster(), silent = TRUE))

  # Stop own pool
  if (exists("cl", envir = .tfu_par_env, inherits = FALSE)) {
    cl <- .tfu_par_env$cl
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
    rm("cl", envir = .tfu_par_env)
  }

  # Reset BLAS/OpenMP threads
  if (exists("old_threads", envir = .tfu_par_env, inherits = FALSE)) {
    olds <- .tfu_par_env$old_threads
    if (!is.null(olds)) {
      if (!is.na(olds[1])) Sys.setenv(OPENBLAS_NUM_THREADS = olds[1])
      if (!is.na(olds[2])) Sys.setenv(MKL_NUM_THREADS = olds[2])
      if (!is.na(olds[3])) Sys.setenv(OMP_NUM_THREADS = olds[3])
    }
    rm("old_threads", envir = .tfu_par_env)
  }

  # Reset RNG to safe reproducible state
  suppressWarnings(RNGkind("L'Ecuyer-CMRG"))

  invisible(TRUE)
}


#' Package load hook: initialize parallel environment
#'
#' Ensures that `.tfu_par_env` exists and is ready when the package
#' is loaded. The actual cluster is only created lazily on demand
#' by `.init_parallel()`.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  assign(".tfu_par_env", new.env(parent = emptyenv()),
    envir = parent.env(environment())
  )
  reg.finalizer(.tfu_par_env, function(e) try(.reset_backend(), silent = TRUE), onexit = TRUE)
}

# # --- Manuell ausfuehren wenn noch kein package ---
# if (!exists(".tfu_par_env", envir = globalenv())) {
#   .tfu_par_env <- new.env(parent = emptyenv())
#   .tfu_par_env$finalizer_set <- FALSE
# }


#' Package unload hook: ensure cluster + RNG cleanup
#' @keywords internal
#' @noRd
.onUnload <- function(libpath) {
  try(.reset_backend(), silent = TRUE)
}


#' Initialize or reuse parallel backend for foreach (robust to interrupts)
#'
#' - Reuse external backends when present
#' - Reuse internal pool if alive; otherwise hard-reset and rebuild it
#' - Set robust RNG (L'Ecuyer-CMRG) and per-worker streams
#' - Pin BLAS/OpenMP threads per worker to avoid oversubscription
#' @keywords internal
#' @noRd
.init_parallel <- function(manage_backend = c("auto", "force_pool", "sequential"),
                           ncpus = parallel::detectCores(logical = TRUE),
                           worker_blas_threads = 1L,
                           seed = 42) {
  checkmate::assert_choice(manage_backend, c("auto", "force_pool", "sequential"))
  checkmate::assert_int(ncpus, lower = 1)
  checkmate::assert_int(worker_blas_threads, lower = 1)
  checkmate::assert_number(seed, null.ok = TRUE)

  # RNG setup
  if (!is.null(seed)) {
    suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
    doRNG::registerDoRNG(seed)
  }

  # Helper: cluster alive?
  .is_alive <- function(cl) {
    if (is.null(cl)) {
      return(FALSE)
    }
    ok <- tryCatch(
      {
        parallel::clusterCall(cl, function() TRUE)
        TRUE
      },
      error = function(e) FALSE
    )
    isTRUE(ok)
  }

  # --- sequential ---
  if (manage_backend == "sequential") {
    foreach::registerDoSEQ()
    return(list(nworkers = 1L, used = "sequential"))
  }

  # --- auto: reuse if possible ---
  if (manage_backend == "auto" && .is_alive(.tfu_par_env$cl)) {
    doParallel::registerDoParallel(.tfu_par_env$cl)
    nworkers <- foreach::getDoParWorkers()
    return(list(nworkers = nworkers, used = "internal-reused"))
  }

  # --- otherwise: create new (force_pool or auto with no alive cluster) ---
  .reset_backend()
  cl <- parallel::makeCluster(ncpus, outfile = "")
  parallel::clusterCall(cl, function(k) {
    Sys.setenv(
      OPENBLAS_NUM_THREADS = k,
      MKL_NUM_THREADS = k,
      OMP_NUM_THREADS = k
    )
    NULL
  }, worker_blas_threads)
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
  doParallel::registerDoParallel(cl)
  .tfu_par_env$cl <- cl

  list(nworkers = ncpus, used = if (manage_backend == "force_pool") "internal-forced" else "internal-new")
}
