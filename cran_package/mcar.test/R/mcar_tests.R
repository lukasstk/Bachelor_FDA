# =============================================================================
# mcar.test - Mean-based MCAR tests & bands
# =============================================================================

#' Internal utilities for mcar.test
#'
#' Helpers used by the exported tests. Not user-facing.
#'
#' @keywords internal
#' @name mcar.test-internal
#'
#' @importFrom tf tfd tf_integrate
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames quantile rnorm median sd
#' @importFrom utils head tail
#'
utils::globalVariables(c("ch"))
NULL

# ---- Shared params topic ---------------------------
#' Common arguments for MCAR mean tests
#'
#' Arguments shared by the exported testing functions.
#' This topic exists only for documentation inheritance.
#'
#' @name mcar_common-params
#' @keywords internal
#'
#' @param fd Optional `tfd`/`tfd_irreg` (tidyfun). If supplied, takes precedence over `X_obs`.
#' @param X_obs Optional numeric matrix (n x m) with NAs (observations).
#' @param groups Optional 2-level grouping vector of length n (logical/character/factor/numeric).
#' @param observed_ratio Threshold in \\[0,1\\] used for auto-grouping when `groups` is missing.
#' @param min_frac Minimum per-time-point coverage per group used to select the subdomain.
NULL

# =============================================================================
# Utilities (internal)  --------------------------------------------------------
# =============================================================================

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
      if (!is.na(olds[2])) Sys.setenv(MKL_NUM_THREADS     = olds[2])
      if (!is.na(olds[3])) Sys.setenv(OMP_NUM_THREADS     = olds[3])
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
         envir = parent.env(environment()))
  .tfu_par_env$finalizer_set <- FALSE
}

# # --- Manuell ausführen wenn noch kein package ---
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

#' Convert group labels to logical
#' @keywords internal
#' @noRd
.groups_to_logical <- function(groups) {
  if (is.factor(groups)) groups <- as.character(droplevels(groups))
  
  if (is.logical(groups)) return(groups)
  
  if (is.character(groups) || is.numeric(groups) || is.integer(groups)) {
    u <- unique(groups)
    if (length(u) != 2L) stop("`groups` must have exactly 2 distinct values.")
    return(groups == u[1])
  }
  
  stop("`groups` must be logical/character/factor/numeric.")
}

#' Coerce `tfd` to dense matrix + grid
#' @keywords internal
#' @noRd
.fd_to_matrix <- function(fd) {
  X_try <- tryCatch(
    { suppressWarnings(as.matrix(fd)) },
    error = function(e) {
      stop("Could not coerce `fd` to matrix via as.matrix(): ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (ncol(X_try) < 2L) stop("`fd` has fewer than 2 support points (ncol < 2).", call. = FALSE)
  
  g <- tf::tf_arg(fd)   
  storage.mode(X_try) <- "numeric"
  
  list(X_obs = X_try, grid = g)
}

#' Subdomain selector (strict + overlap fallback)
#' @keywords internal
#' @noRd
.limit_subdomain <- function(O_mat, group_A, min_frac = 0.10) {
  stopifnot(is.matrix(O_mat), length(group_A) == nrow(O_mat))
  n  <- nrow(O_mat)
  IA <- as.numeric(as.logical(group_A))
  IB <- 1 - IA
  cA <- colSums(O_mat * IA)
  cB <- colSums(O_mat * IB)

  idx_strict <- which(pmin(cA, cB) > n * min_frac)
  if (length(idx_strict) >= 2L) {
    return(list(idx = idx_strict, min_frac_used = min_frac, fallback = NULL))
  }
  stop(
    paste0(
      "No suitable subdomain found: strict criterion (min_frac = ", format(min_frac), 
      ") not satisfied at \u2265 2 time points. Adjust grouping or 'min_frac'."
    ),
    call. = FALSE
  )
}

#' Prepare inputs from `tfd` or matrix
#' @keywords internal
#' @noRd
.prepare_inputs <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1) {
  if (!is.null(fd)) {
    conv  <- .fd_to_matrix(fd)
    X_obs <- conv$X_obs
    grid_vec <- conv$grid
  } else {
    if (is.null(X_obs)) stop("Either `fd` or `X_obs` must be supplied.")
    grid_vec <- seq(0, 1, length.out = ncol(X_obs))
  }

  n <- nrow(X_obs)
  O_mat <- 1L * (!is.na(X_obs))

  if (!is.null(groups)) {
    if (length(groups) != n) stop("`groups` must have length nrow(X_obs).")
    if (any(is.na(groups))) stop("`groups` must not contain NAs.")
    group_A <- .groups_to_logical(groups)

    obs_frac <- rowMeans(O_mat != 0)
    meanA <- mean(obs_frac[group_A], na.rm = TRUE)
    meanB <- mean(obs_frac[!group_A], na.rm = TRUE)
    if (meanA < meanB) {
      group_A <- !group_A
      message(sprintf(
        "Note: swapped labels so Group A is the more complete group (mean A=%.3f, B=%.3f).",
        meanA, meanB
      ))
    }
  } else {
    group_A <- rowMeans(O_mat != 0) >= observed_ratio
    if (sum(group_A) == 0L || sum(!group_A) == 0L) {
      stop("Auto-grouping failed: one group empty. Discard this run or adjust observed_ratio.")
    }
  }

  list(
    X_obs   = X_obs,
    O_mat   = O_mat,
    group_A = as.logical(group_A),
    grid    = grid_vec
  )
}

#' Trapezoidal integration weights
#' @keywords internal
#' @noRd
.trapezoid_weights <- function(x) {
  m <- length(x)
  if (m == 1L) return(1)
  w <- numeric(m)
  w[1] <- (x[2]-x[1])/2
  w[m] <- (x[m]-x[m-1])/2
  if (m > 2L) w[2:(m-1)] <- (x[3:m] - x[1:(m-2)])/2
  w
}

#' Available-mean estimators by group
#' @keywords internal
#' @noRd
.group_mean_estimators <- function(X_obs, O_mat, group_A) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A); IB <- 1 - IA
  pA_hat <- colSums(O_mat * IA) / n
  pB_hat <- colSums(O_mat * IB) / n
  muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_hat)
  list(muA = muA_hat, muB = muB_hat, pA = pA_hat, pB = pB_hat)
}

#' Corrected covariance under partial observation
#' @keywords internal
#' @noRd
.covariance_estimator <- function(X_obs, O_mat, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X_obs)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  Xtilde <- X_obs
  Xtilde[ group_A, ] <- sweep(X_obs[ group_A, , drop = FALSE], 2, muA_hat, `-`)
  Xtilde[!group_A, ] <- sweep(X_obs[!group_A, , drop = FALSE], 2, muB_hat, `-`)
  Xtilde[is.na(Xtilde)] <- 0
  A_resid <- sweep((Xtilde * O_mat) * IA, 2, pA_hat, "/")
  B_resid <- sweep((Xtilde * O_mat) * IB, 2, pB_hat, "/")
  K_hat <- (crossprod(A_resid) + crossprod(B_resid)) / n
  (K_hat + t(K_hat)) / 2
}

#' KL basis from covariance 
#' @keywords internal
#' @noRd
.kl_decomposition <- function(K, grid, tol = sqrt(.Machine$double.eps)) {
  w  <- .trapezoid_weights(grid)
  sw <- sqrt(w)
  S <- (sw %o% sw) * K
  S <- (S + t(S)) / 2
  ev_test <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev_test < -tol * max(1, abs(ev_test[1])))) {
    stop(sprintf("Weighted covariance is not PSD (min eigenvalue = %.4g).", min(ev_test)), call. = FALSE)
  }
  ev  <- eigen(S, symmetric = TRUE)
  lam <- ev$values
  U   <- ev$vectors
  phi   <- sweep(U, 1, sw, "/")
  norms <- sqrt(colSums(phi^2 * w))
  phi   <- sweep(phi, 2, norms, "/")
  list(lam = lam, phi = phi, w = w)
}

#' Functional confidence bands
#'
#' Pointwise bands for L2, simultaneous bands for Sup.
#'
#' @keywords internal
#' @noRd
.confidence_bands <- function(stat, diff, W, n, alpha, grid, method = c("asymptotic","bootstrap")) {
  method <- match.arg(method)
  
  if (identical(stat, "L2")) {
    if (method == "asymptotic") {
      # --- Asymptotic (Normal approximation based on test stat draws) ---
      se <- sd(W)
      halfwidth <- qnorm(1 - alpha/2) * se
      lower <- diff - halfwidth
      upper <- diff + halfwidth
    } else {
      # --- Bootstrap (pointwise bands based on bootstrap resamples) ---
      # W expected to be a matrix: B x m (bootstrap samples × grid points)
      lower <- apply(W, 2, quantile, probs = alpha/2, na.rm = TRUE)
      upper <- apply(W, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
    }
    
    band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
    
    return(list(
      type    = "pointwise",
      lower   = lower,
      upper   = upper,
      band    = band,
      q_alpha = NULL,
      alpha   = alpha
    ))
  }
  
  if (identical(stat, "D")) {
    # --- Simultaneous bands ---
    q_alpha   <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
    halfwidth <- q_alpha / sqrt(n)
    lower <- diff - halfwidth
    upper <- diff + halfwidth
    
    band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
    
    return(list(
      type    = "simultaneous",
      lower   = lower,
      upper   = upper,
      band    = band,
      q_alpha = q_alpha,
      alpha   = alpha
    ))
  }
  
  stop("Unknown stat type in .confidence_bands(): must be 'L2' or 'D'.")
}


#' Build extended htest object
#' @keywords internal
#' @noRd
.create_output <- function(stat_name, stat_value, p_value, method, data_name,
                                estimate = NULL, parameter = NULL,
                                null.value = 0, alternative = "two.sided",
                                alpha = NA_real_) {
  out <- list(
    statistic   = stats::setNames(as.numeric(stat_value), stat_name),
    parameter   = parameter,
    p.value     = if (!is.na(p_value)) as.numeric(p_value) else NA_real_,
    estimate    = estimate,
    null.value  = c("difference in mean functions" = null.value),
    alternative = alternative,
    method      = method,
    data.name   = data_name,
    alpha       = alpha
  )
  class(out) <- "htest"
  out
}


#' Initialize or reuse parallel backend for foreach (robust to interrupts)
#'
#' - Reuse external backends when present
#' - Reuse internal pool if alive; otherwise hard-reset and rebuild it
#' - Set robust RNG (L'Ecuyer-CMRG) and per-worker streams
#' - Pin BLAS/OpenMP threads per worker to avoid oversubscription
#' @keywords internal
#' @noRd
.init_parallel <- function(manage_backend = c("auto","force_pool","sequential"),
                           ncpus = parallel::detectCores(logical = TRUE),
                           worker_blas_threads = 1L,
                           seed = 42,
                           parallel_flag = TRUE) {
  manage_backend <- match.arg(manage_backend)
  ncpus <- max(1L, min(as.integer(ncpus), parallel::detectCores(logical = TRUE)))
  
  # --- RNG setup (session level, once) ---
  if (!is.null(seed)) {
    suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
    doRNG::registerDoRNG(seed)
  }
  
  # --- Helper: check if cluster alive ---
  .is_alive <- function(cl) {
    if (is.null(cl)) return(FALSE)
    ok <- tryCatch({ parallel::clusterCall(cl, function() TRUE); TRUE },
                   error = function(e) FALSE)
    isTRUE(ok)
  }
  
  # --- Sequential path (no cluster) ---
  if (!isTRUE(parallel_flag) || identical(manage_backend, "sequential")) {
    foreach::registerDoSEQ()
    return(list(nworkers = 1L, used = "sequential"))
  }
  
  # --- External backend present? -> reuse ---
  existing_workers <- tryCatch(foreach::getDoParWorkers(), error = function(e) 1L)
  existing_backend <- tryCatch(foreach::getDoParName(),    error = function(e) "doSEQ")
  if (manage_backend == "auto" && existing_workers > 1L && existing_backend != "doSEQ") {
    return(list(nworkers = existing_workers, used = paste0("external:", existing_backend)))
  }
  
  # --- Internal backend: reuse if alive (unless force_pool) ---
  if (.is_alive(.tfu_par_env$cl) && manage_backend != "force_pool") {
    doParallel::registerDoParallel(.tfu_par_env$cl)
    nworkers <- foreach::getDoParWorkers()
    return(list(nworkers = nworkers, used = "internal-reused"))
  }
  
  # --- Hard reset (after interrupts etc.) ---
  .reset_backend()
  
  # --- Create fresh internal cluster ---
  cl <- parallel::makeCluster(ncpus, outfile = "")
  
  # Pin BLAS/OpenMP threads per worker
  .tfu_par_env$old_threads <- Sys.getenv(
    c("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"), unset = NA
  )
  parallel::clusterCall(cl, function(k) {
    Sys.setenv(OPENBLAS_NUM_THREADS = k,
               MKL_NUM_THREADS      = k,
               OMP_NUM_THREADS      = k)
    NULL
  }, worker_blas_threads)
  
  # Worker RNG streams (only if seed provided)
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  
  # Register backend
  doParallel::registerDoParallel(cl)
  .tfu_par_env$cl <- cl
  
  # One-time finalizer for cleanup
  if (!isTRUE(.tfu_par_env$finalizer_set)) {
    reg.finalizer(.tfu_par_env, function(e) {
      if (!is.null(e$cl)) { try(parallel::stopCluster(e$cl), silent = TRUE); e$cl <- NULL }
      olds <- e$old_threads
      if (!is.null(olds)) {
        if (!is.na(olds[1])) Sys.setenv(OPENBLAS_NUM_THREADS = olds[1])
        if (!is.na(olds[2])) Sys.setenv(MKL_NUM_THREADS     = olds[2])
        if (!is.na(olds[3])) Sys.setenv(OMP_NUM_THREADS     = olds[3])
        e$old_threads <- NULL
      }
    }, onexit = TRUE)
    .tfu_par_env$finalizer_set <- TRUE
  }
  
  list(nworkers = ncpus, used = "internal-new")
}


# =============================================================================
# Exported tests  --------------------------------------------------------------
# =============================================================================

#' Asymptotic L2 test for MCAR (optional pointwise bands)
#'
#' Tests \eqn{H_0:\ \mu_A=\mu_B} via \eqn{T_{\mu,L2}=n\lVert \hat\mu_A-\hat\mu_B\rVert^2_{L2}}.
#' p-values are obtained from a KL-mixture; the number of components is chosen by FVE.
#' Optionally, rough pointwise confidence bands for the mean difference can be computed.
#'
#' @inheritParams mcar_common-params
#' @param fve Fraction of variance explained (0–1) to choose \eqn{q}.
#' @param B Number of Monte Carlo draws for the KL-mixture.
#' @param seed RNG seed.
#' @param alpha Significance level (only relevant if \code{compute_bands = TRUE}).
#' @param compute_bands Logical: if TRUE, compute rough pointwise confidence bands.
#' @param bands_only Logical: if TRUE, return only band information instead of a full \code{htest}.
#' @return `htest` (with extras) or a light band list.
#' @examples
#' # Brownian toy example (quick)
#' set.seed(1)
#' m <- 30; n <- 20
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) { d <- diff(g)[1]; c(0, cumsum(rnorm(m-1, sd = sqrt(d)))) }
#' X  <- t(replicate(n, bm(grid)))
#' # Introduce simple censoring (MNAR): observe only when -1 < X < 2
#' O  <- 1L * (X > -1 & X < 2); X[O == 0L] <- NA_real_
#' # Run asymptotic L2 test with bands
#' res <- asym_mean_L2_test(X_obs = X, B = 1000, seed = 1, compute_bands = TRUE)
#' res$p.value
#' @export
asym_mean_L2_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, B = 5000,
                              min_frac = 0.10, seed = NULL, alpha = 0.05,
                              compute_bands = TRUE, bands_only = FALSE) {
  # --- Input preparation ---
  prep <- .prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  # --- Subdomain selection ---
  sub <- .limit_subdomain(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  # --- Group means ---
  est <- .group_mean_estimators(X, O, group_A)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA,  nrow = 1), arg = g)
  muB_tfd  <- tf::tfd(matrix(muB,  nrow = 1), arg = g)
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  
  # --- Test statistic ---
  T_L2 <- n * tf::tf_integrate(diff_tfd^2, arg = g)
  
  # --- Covariance + KL truncation ---
  K  <- .covariance_estimator(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .kl_decomposition(K, g); lam <- KL$lam
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]
  
  # --- Monte Carlo draws ---
  Z <- matrix(rnorm(q * B), nrow = q)
  W <- colSums((Z^2) * lam_q)
  
  # --- p-value ---
  p <- (sum(W >= T_L2) + 1) / (length(W) + 1)
  
  # --- Bands ---
  if (isTRUE(compute_bands)) {
    bands <- .confidence_bands("L2", diff, W, n, alpha, g, method = "asymptotic")
  } else {
    bands <- list(lower = NULL, upper = NULL, band = NULL, q_alpha = NA_real_)
  }
  
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  
  # --- bands_only output ---
  if (isTRUE(bands_only)) {
    return(list(
      grid = g, diff = diff, lower = bands$lower, upper = bands$upper, 
      q_alpha = bands$q_alpha, alpha = alpha,
      estimate = list(muA = muA_tfd, muB = muB_tfd, diff = diff_tfd),
      band_tfd = bands$band
    ))
  }
  
  # --- Full htest output ---
  h <- .create_output("T_{mu,L2}", T_L2, p,
                           "L2 test",
                           data_name,
                           estimate = list(muA  = muA_tfd,
                                           muB  = muB_tfd,
                                           diff = diff_tfd,
                                           band = bands$band),
                           parameter = c(q = q, m = length(g)))
  h$bands <- bands
  h
}

#' Asymptotic supremum test for MCAR (optional simultaneous bands)
#'
#' Tests \eqn{H_0:\ \mu_A=\mu_B} via \eqn{T_{\mu,D}=\sqrt{n}\,\lVert \hat\mu_A-\hat\mu_B\rVert_\infty}.
#' KL-based GP approximation; with `compute_bands=TRUE` also returns simultaneous bands.
#'
#' @inheritParams mcar_common-params
#' @param fve Fraction of variance explained (0-1) to select KL components.
#' @param B Number of Monte Carlo draws for the GP approx.
#' @param seed RNG seed.
#' @param alpha Significance level for bands.
#' @param compute_bands Logical: also compute simultaneous bands?
#' @param bands_only Logical: return only band info?
#' @return `htest` (with extras) or a light band list.
#' @examples
#' set.seed(1)
#' m <- 30; n <- 20
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) { d <- diff(g)[1]; c(0, cumsum(rnorm(m-1, sd = sqrt(d)))) }
#' X  <- t(replicate(n, bm(grid)))
#' O  <- 1L * (X > -1 & X < 2); X[O == 0L] <- NA_real_
#' res <- asym_mean_sup_test(X_obs = X, B = 1000, compute_bands = TRUE, seed = 1)
#' # simple visual:
#' plot(res$grid, res$diff, type = "l")
#' lines(res$grid, res$lower, lty = 2)
#' lines(res$grid, res$upper, lty = 2)
#' @export
asym_mean_sup_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                               fve = 0.99, B = 5000,
                               min_frac = 0.10, seed = NULL, alpha = 0.05,
                               compute_bands = TRUE, bands_only = FALSE) {
  prep <- .prepare_inputs(fd, X_obs, groups, observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .limit_subdomain(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx; g <- grid[idx]
  X  <- X_obs[, idx, drop = FALSE]; O <- O_mat[, idx, drop = FALSE]
  
  est <- .group_mean_estimators(X, O, group_A)
  muA <- est$muA; muB <- est$muB; diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA,  nrow = 1), arg = g)
  muB_tfd  <- tf::tfd(matrix(muB,  nrow = 1), arg = g)
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  
  T_D <- sqrt(n) * max(abs(diff))
  
  K  <- .covariance_estimator(X, O, group_A, muA, muB, est$pA, est$pB)
  KL <- .kl_decomposition(K, g); lam <- KL$lam; phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam); q <- which(cum >= fve)[1]; q <- max(1L, q)
  lam_q <- lam[seq_len(q)]; phi_q <- phi[, seq_len(q), drop = FALSE]
  A <- sweep(phi_q, 2, sqrt(lam_q), "*")  # m x q
  
  Z <- matrix(rnorm(q * B), nrow = q)     # q x B
  gp_vals <- A %*% Z                      # m x B
  W <- apply(abs(gp_vals), 2, max)
  
  p <- (sum(W >= T_D) + 1) / (length(W) + 1)
  
  if (isTRUE(compute_bands)) {
    bands <- .confidence_bands("D", diff, W, n, alpha, g)
  } else {
    bands <- list(lower = NULL, upper = NULL, band = NULL, q_alpha = NA_real_)
  }
  
  data_name <- if (!is.null(fd)) "fd" else "X_obs"
  
  if (isTRUE(bands_only)) {
    return(list(
      grid=g, diff=diff, lower=bands$lower, upper=bands$upper, q_alpha=bands$q_alpha,
      idx=idx, min_frac_used=sub$min_frac_used, fallback=sub$fallback,
      alpha=alpha, estimate=list(muA=muA_tfd, muB=muB_tfd, diff=diff_tfd)
    ))
  }
  
  h <- .create_output("T_{mu,D}", T_D, p,
                           "Supremum test",
                           data_name,
                           estimate = list(muA  = muA_tfd,
                                           muB  = muB_tfd,
                                           diff = diff_tfd,
                                           band = bands$band),
                           parameter = c(q = q, m = length(g)))
  h$bands <- bands
  h
}

#' Bootstrap mean test (L2/Supremum) with optional bands
#'
#' Returns bootstrap p-values for L2 or Supremum and, optionally, confidence
#' bands for the difference curve. Parallelization via foreach/doParallel
#' with reproducible seeding via doRNG.
#'
#' @inheritParams mcar_common-params
#' @param B Number of bootstrap iterations.
#' @param alpha Significance level for bands.
#' @param parallel Logical: run in parallel?
#' @param ncpus Number of workers for an internal cluster (if created).
#' @param seed RNG seed (passed to doRNG).
#' @param stat `"L2"`, `"D"` oder `c("L2","D")`.
#' @param compute_bands Compute confidence bands?
#' @param return_boot Attach bootstrap statistics?
#' @param chunk_size Number of bootstrap replicates per foreach task.
#' @param manage_backend Backend control (`"auto"`, `"force_pool"`, `"sequential"`).
#' @param worker_blas_threads BLAS/OpenMP threads per worker (internal pool only).
#' @param cleanup Logical: if TRUE, detach foreach/doRNG and stop the internal pool on exit.
#' @param bands_only Return only the band payload?
#' @return `htest` (with extras) or a list of both tests.
#' @export
boot_mean_test <- function(fd = NULL, X_obs = NULL, groups = NULL, observed_ratio = 1,
                           B = 5000,
                           min_frac = 0.10, alpha = 0.05,
                           parallel = TRUE,
                           ncpus = parallel::detectCores(logical = TRUE),
                           seed = NULL,
                           stat = c("L2", "D"),
                           compute_bands = TRUE,
                           return_boot = FALSE,
                           chunk_size = NULL,
                           manage_backend = c("auto","force_pool","sequential"),
                           worker_blas_threads = 1L,
                           bands_only = FALSE,
                           cleanup = FALSE) {
  
  if (isTRUE(cleanup)) {
    on.exit({
      try(foreach::registerDoSEQ(), silent = TRUE)
      try(doRNG::registerDoRNG(NULL), silent = TRUE)
      try(.reset_backend(), silent = TRUE)
    }, add = TRUE)
  }
  
  # --- Stat argument prüfen ---
  stat <- unique(stat)
  if (!all(stat %in% c("L2","D"))) {
    stop("Invalid 'stat' argument: must be 'L2', 'D' or c('L2','D').", call. = FALSE)
  }
  
  manage_backend <- match.arg(manage_backend)
  
  prep <- .prepare_inputs(fd = fd, X_obs = X_obs, groups = groups, observed_ratio = observed_ratio)
  X_obs <- prep$X_obs; O_mat <- prep$O_mat; group_A <- prep$group_A; grid <- prep$grid
  n <- nrow(X_obs)
  
  sub <- .limit_subdomain(O_mat, group_A, min_frac = min_frac)
  idx <- sub$idx
  X   <- X_obs[, idx, drop = FALSE]
  O   <- O_mat[, idx, drop = FALSE]
  g   <- grid[idx]
  
  est  <- .group_mean_estimators(X, O, group_A)
  muA  <- est$muA; muB <- est$muB
  diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA,  nrow = 1), arg = g)
  muB_tfd  <- tf::tfd(matrix(muB,  nrow = 1), arg = g)
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = g)
  
  w <- .trapezoid_weights(g)
  T_L2 <- n * sum((diff^2) * w)
  T_D  <- sqrt(n) * max(abs(diff))
  
  IA <- as.numeric(as.logical(group_A)); IB <- 1 - IA
  A_idx <- which(IA == 1); B_idx <- which(IB == 1)
  X_cent <- X
  if (length(A_idx)) X_cent[A_idx, ] <- sweep(X[A_idx, , drop = FALSE], 2, muA, `-`)
  if (length(B_idx)) X_cent[B_idx, ] <- sweep(X[B_idx, , drop = FALSE], 2, muB, `-`)
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- inner bootstrap runner ---
  .run_boot <- function(manage_backend_mode = manage_backend, ncpus_eff = ncpus) {
    be <- .init_parallel(
      manage_backend      = manage_backend_mode,
      ncpus               = ncpus_eff,
      worker_blas_threads = worker_blas_threads,
      seed                = seed,
      parallel_flag       = parallel
    )
    nworkers <- be$nworkers
    
    # Chunking
    cs <- chunk_size
    if (is.null(cs) || !is.finite(cs) || cs < 1L) {
      cs <- max(50L, ceiling(B / (3L * max(1L, nworkers))))
    } else {
      cs <- as.integer(cs)
    }
    idx_chunks <- split(seq_len(B), ceiling(seq_len(B) / cs))
    
    # capture for worker export
    n_       <- n
    X_cent_  <- X_cent
    O_       <- O
    group_A_ <- group_A
    w_       <- w
    
    boot_list <- foreach::foreach(
      ch = idx_chunks, .inorder = FALSE
    ) %dorng% {
      out   <- matrix(NA_real_, nrow = length(ch), ncol = 3L,
                      dimnames = list(NULL, c("L2","D","redraws")))
      diffs <- matrix(NA_real_, nrow = length(ch), ncol = ncol(X_cent_))
      
      for (ii in seq_along(ch)) {
        redraws <- 0L
        repeat {
          samp_idx <- sample.int(n_, n_, replace = TRUE)
          gA_samp  <- group_A_[samp_idx]
          IA <- as.numeric(gA_samp); IB <- 1 - IA
          OA <- (O_[samp_idx, , drop = FALSE] * IA)
          OB <- (O_[samp_idx, , drop = FALSE] * IB)
          pA <- colSums(OA) / n_
          pB <- colSums(OB) / n_
          if (all(pA > 0) && all(pB > 0)) break
          redraws <- redraws + 1L
        }
        Xs <- X_cent_[samp_idx, , drop = FALSE]
        XA <- Xs; XA[IB == 1, ] <- NA
        numA <- colSums(replace(XA, is.na(XA), 0))
        muA_b <- numA / (n_ * pA)
        XB <- Xs; XB[IA == 1, ] <- NA
        numB <- colSums(replace(XB, is.na(XB), 0))
        muB_b <- numB / (n_ * pB)
        d_b <- muA_b - muB_b
        
        out[ii, ] <- c(
          L2 = n_ * sum((d_b^2) * w_),
          D  = sqrt(n_) * max(abs(d_b)),
          redraws = redraws
        )
        diffs[ii, ] <- d_b
      }
      list(stats = out, diffs = diffs)
    }
    
    boot_mat   <- do.call(rbind, lapply(boot_list, `[[`, "stats"))
    boot_diffs <- do.call(rbind, lapply(boot_list, `[[`, "diffs"))
    list(boot_mat = boot_mat, boot_diffs = boot_diffs)
  }
  
  # --- robust run with retry ---
  boot_res <- tryCatch(
    .run_boot(manage_backend_mode = manage_backend, ncpus_eff = ncpus),
    interrupt = function(e) {
      .reset_backend()
      stop("Aborted by user: backend cleaned up; re-execution is possible immediately.", call. = FALSE)
    },
    error = function(e) e
  )
  if (inherits(boot_res, "error")) {
    .reset_backend()
    boot_res <- tryCatch(
      .run_boot(manage_backend_mode = "force_pool",
                ncpus_eff = max(1L, parallel::detectCores(logical = TRUE) - 1L)),
      interrupt = function(e) {
        .reset_backend()
        stop("Aborted by user: backend cleaned up; re-execution is immediately possible.", call. = FALSE)
      },
      error = function(e) e
    )
  }
  if (inherits(boot_res, "error")) stop(boot_res)
  
  boot_mat   <- boot_res$boot_mat
  boot_diffs <- boot_res$boot_diffs
  
  # clean bootstrap stats
  boot_L2 <- boot_mat[, "L2"]; boot_D <- boot_mat[, "D"]
  total_redraws <- sum(boot_mat[, "redraws"])
  n_bad <- sum(is.na(boot_L2) | is.na(boot_D))
  if (n_bad > 0) {
    warning(sprintf("Bootstrap: %d iterations skipped (%.1f%%).", n_bad, 100 * n_bad / B))
  }
  boot_L2 <- boot_L2[!is.na(boot_L2)]
  boot_D  <- boot_D[!is.na(boot_D)]
  
  p_L2 <- (sum(boot_L2 >= T_L2) + 1) / (length(boot_L2) + 1)
  p_D  <- (sum(boot_D  >= T_D ) + 1) / (length(boot_D)  + 1)
  
  results <- list()
  
  if ("L2" %in% stat) {
    bands <- if (isTRUE(compute_bands)) 
      .confidence_bands("L2", diff, boot_diffs, n, alpha, g, method="bootstrap") 
    else list(lower=NULL,upper=NULL,band=NULL,q_alpha=NA_real_)
    
    results$L2 <- .create_output("T_{mu,L2}", T_L2, p_L2,
                                      "Bootstrap mean test (L2)", 
                                      if (!is.null(fd)) "fd (tfd/tfd_irreg)" else "X_obs",
                                      estimate = list(muA=muA_tfd, muB=muB_tfd, diff=diff_tfd, band=bands$band),
                                      parameter = c(m=length(g), B=B-n_bad),
                                      null.value = 0, alternative="two.sided")
    results$L2$q_alpha <- bands$q_alpha; results$L2$lower <- bands$lower; results$L2$upper <- bands$upper
  }
  
  if ("D" %in% stat) {
    bands <- if (isTRUE(compute_bands)) 
      .confidence_bands("D", diff, boot_D, n, alpha, g, method="bootstrap") 
    else list(lower=NULL,upper=NULL,band=NULL,q_alpha=NA_real_)
    
    results$D <- .create_output("T_{mu,D}", T_D, p_D,
                                     "Bootstrap mean test (supremum)", 
                                     if (!is.null(fd)) "fd (tfd/tfd_irreg via tf_gather)" else "X_obs",
                                     estimate = list(muA=muA_tfd, muB=muB_tfd, diff=diff_tfd, band=bands$band),
                                     parameter = c(m=length(g), B=B-n_bad),
                                     null.value = 0, alternative="two.sided")
    results$D$q_alpha <- bands$q_alpha; results$D$lower <- bands$lower; results$D$upper <- bands$upper
  }
  
  # Meta-Infos
  for (nm in names(results)) {
    results[[nm]]$grid <- g; results[[nm]]$muA <- muA; results[[nm]]$muB <- muB; results[[nm]]$diff <- diff
    results[[nm]]$idx  <- idx; results[[nm]]$min_frac_used <- sub$min_frac_used; results[[nm]]$fallback <- sub$fallback
    results[[nm]]$alpha <- alpha
    results[[nm]]$band_tfd <- results[[nm]]$estimate$band
    results[[nm]]$total_redraws <- total_redraws
    results[[nm]]$skip_rate <- total_redraws / (B + total_redraws)
    if (isTRUE(return_boot)) { results[[nm]]$boot_L2 <- boot_L2; results[[nm]]$boot_D <- boot_D }
  }
  
  if (length(results) == 1) return(results[[1]])
  return(results)
}

