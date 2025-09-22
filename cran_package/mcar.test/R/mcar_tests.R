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
#' @importFrom checkmate assert_matrix assert_numeric assert_atomic_vector assert_class assert_logical assert_number assert_true
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
#' @param fd Optional `tfd`/`tfd_irreg` (tidyfun). If supplied, takes precedence over `X`.
#' @param X Optional numeric matrix (n x m) with NAs (observations).
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
  reg.finalizer(.tfu_par_env, function(e) try(.reset_backend(), silent = TRUE), onexit = TRUE)
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
  checkmate::assert_atomic_vector(groups, any.missing = FALSE)
  if (is.factor(groups)) groups <- as.character(droplevels(groups))
  
  # Sicherstellen: genau 2 verschiedene Gruppen
  checkmate::assert(
    length(unique(groups)) == 2,
    .var.name = "groups"
  )
  
  if (is.logical(groups)) {
    return(groups)
  } else if (is.character(groups) || is.numeric(groups) || is.integer(groups)) {
    u <- unique(groups)
    return(groups == u[1])
  }
  
  stop("`groups` must be logical/character/factor/numeric.")
}

#' Coerce `tfd` to dense matrix + grid
#' @keywords internal
#' @noRd
.fd_tOrix <- function(fd) {
  checkmate::assert_class(fd, c("tfd", "tfd_irreg"))
  
  all_grids <- tf::tf_arg(fd)
  if (is.list(all_grids)) {
    g <- sort(unique(unlist(all_grids)))
  } else {
    g <- all_grids
  }
  
  X_try <- suppressWarnings(as.matrix(fd))
  checkmate::assert_matrix(X_try, mode = "numeric", min.cols = 2)
  
  list(X = X_try, grid = g)
}

#' Subdomain selector (strict + overlap fallback)
#' @keywords internal
#' @noRd
.limit_subdomain <- function(O, group_A, min_frac = 0.10) {
  checkmate::assert_matrix(O, any.missing = FALSE)
  checkmate::assert_logical(group_A, len = nrow(O))
  checkmate::assert_number(min_frac, lower = 0, upper = 1)
  
  n  <- nrow(O)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  cA <- colSums(O * IA)
  cB <- colSums(O * IB)
  
  idx_strict <- which(pmin(cA, cB) > n * min_frac)
  if (length(idx_strict) >= 2L) {
    return(list(idx = idx_strict, min_frac_used = min_frac, fallback = NULL))
  }
  
  stop(
    paste0(
      "No suitable subdomain found: strict criterion (min_frac = ", format(min_frac), 
      ") not satisfied at >= 2 time points."
    ),
    call. = FALSE
  )
}


#' Prepare inputs from `tfd` or matrix
#' @keywords internal
#' @noRd
.prepare_inputs <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1) {
  checkmate::assert_number(observed_ratio, lower = 0, upper = 1)
  
  if (!is.null(fd)) {
    conv  <- .fd_tOrix(fd)
    X <- conv$X
    grid_vec <- conv$grid
  } else {
    checkmate::assert_matrix(X, mode = "numeric", min.rows = 1, min.cols = 2)
    grid_vec <- seq(0, 1, length.out = ncol(X))
  }
  
  n <- nrow(X)
  O <- 1L * (!is.na(X))
  
  if (!is.null(groups)) {
    checkmate::assert_atomic_vector(groups, len = n, any.missing = FALSE)
    group_A <- .groups_to_logical(groups)
    
    # Label swap falls nötig
    obs_frac <- rowMeans(O != 0)
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
    group_A <- rowMeans(O != 0) >= observed_ratio
    if (sum(group_A) == 0L || sum(!group_A) == 0L) {
      stop("Auto-grouping failed: one group empty. Adjust `observed_ratio`.")
    }
  }
  
  list(
    X   = X,
    O   = O,
    group_A = as.logical(group_A),
    grid    = grid_vec
  )
}


#' Trapezoidal integration weights
#' @keywords internal
#' @noRd
.trapezoid_weights <- function(x) {
  checkmate::assert_numeric(x, any.missing = FALSE, min.len = 1, sorted = TRUE)
  
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
.group_mean_estimators <- function(X, O, group_A) {
  n  <- nrow(X)
  IA <- as.numeric(group_A); IB <- 1 - IA
  pA_hat <- colSums(O * IA) / n
  pB_hat <- colSums(O * IB) / n
  muA_hat <- colSums(X * IA, na.rm = TRUE) / (n * pA_hat)
  muB_hat <- colSums(X * IB, na.rm = TRUE) / (n * pB_hat)
  list(muA = muA_hat, muB = muB_hat, pA = pA_hat, pB = pB_hat)
}

#' Corrected covariance under partial observation
#' @keywords internal
#' @noRd
.covariance_estimator <- function(X, O, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n  <- nrow(X)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  Xtilde <- X
  Xtilde[ group_A, ] <- sweep(X[ group_A, , drop = FALSE], 2, muA_hat, `-`)
  Xtilde[!group_A, ] <- sweep(X[!group_A, , drop = FALSE], 2, muB_hat, `-`)
  Xtilde[is.na(Xtilde)] <- 0
  A_resid <- sweep((Xtilde * O) * IA, 2, pA_hat, "/")
  B_resid <- sweep((Xtilde * O) * IB, 2, pB_hat, "/")
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
.confidence_bands <- function(stat, diff, W, n, alpha, grid,
                              method = c("asymptotic","bootstrap")) {
  method <- match.arg(method)
  
  if (identical(stat, "L2")) {
    if (method == "asymptotic") {
      se <- sd(W)
      halfwidth <- qnorm(1 - alpha/2) * se
      lower <- diff - halfwidth
      upper <- diff + halfwidth
    } else {
      lower <- apply(W, 2, quantile, probs = alpha/2, na.rm = TRUE)
      upper <- apply(W, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
    }
    band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
    return(list(type = "pointwise", band = band,
                lower = lower, upper = upper,  alpha = alpha, grid = grid))
  }
  
  if (identical(stat, "D")) {
    q_alpha   <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))
    halfwidth <- q_alpha / sqrt(n)
    lower <- diff - halfwidth
    upper <- diff + halfwidth
    band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
    return(list(type = "simultaneous", band = band,
                lower = lower, upper = upper, alpha = alpha, grid = grid))
  }
  
  stop("Unknown stat type in .confidence_bands(): must be 'L2' or 'D'.")
}



#' Build extended htest object
#' @keywords internal
#' @noRd
.create_output <- function(stat_name, stat_value, p_value, method, data_name,
                           estimate = NULL, parameter = NULL,
                           null.value = 0, alternative = "two.sided",
                           alpha = NA_real_, bands = NULL, grid = NULL) {
  out <- list(
    statistic   = stats::setNames(as.numeric(stat_value), stat_name),
    parameter   = parameter,
    p.value     = if (!is.na(p_value)) as.numeric(p_value) else NA_real_,
    estimate    = estimate,
    null.value  = c("difference in mean functions" = null.value),
    alternative = alternative,
    method      = method,
    data.name   = data_name,
    alpha       = alpha,
    bands       = bands,   
    grid        = grid     
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
                           seed = 42) {
  checkmate::assert_choice(manage_backend, c("auto","force_pool","sequential"))
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
    if (is.null(cl)) return(FALSE)
    ok <- tryCatch({ parallel::clusterCall(cl, function() TRUE); TRUE },
                   error = function(e) FALSE)
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
    Sys.setenv(OPENBLAS_NUM_THREADS = k,
               MKL_NUM_THREADS      = k,
               OMP_NUM_THREADS      = k)
    NULL
  }, worker_blas_threads)
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
  doParallel::registerDoParallel(cl)
  .tfu_par_env$cl <- cl
  
  list(nworkers = ncpus, used = if (manage_backend == "force_pool") "internal-forced" else "internal-new")
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
#' @param n_sim Number of Monte Carlo draws for the KL-mixture.
#' @param seed RNG seed.
#' @param alpha Significance level (only relevant if \code{compute_bands = TRUE}).
#' @param compute_bands Logical: if TRUE, compute rough pointwise confidence bands.
#' @param bands_only Logical: if TRUE, return only band information instead of a full \code{htest}.
#' @return `htest` (with extras) or a light band list.
#' @examples
#' set.seed(1)
#' m <- 50; n <- 40
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) { d <- diff(g)[1]; c(0, cumsum(rnorm(length(g)-1, sd = sqrt(d)))) }
#' X  <- t(replicate(n, bm(grid)))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X[O == 0L] <- NA_real_
#'
#' # Asymptotic L2 test
#' res_L2 <- asym_mean_L2_test(
#'   X = X,
#'   n_sim = 2000,   
#'   seed  = 123
#' )
#' res_L2$p.value
#' @export
asym_mean_L2_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                              fve = 0.99, n_sim = 5000,
                              min_frac = 0.10, seed = NULL, alpha = 0.05,
                              compute_bands = TRUE, bands_only = FALSE) {
  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)
  
  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid <- grid[idx]
  X_sub  <- X[, idx, drop = FALSE]
  O_sub <- O[, idx, drop = FALSE]
  
  est <- .group_mean_estimators(X_sub, O_sub, group_A)
  muA <- est$muA
  muB <- est$muB
  diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA,  nrow = 1), arg = subgrid)
  muB_tfd  <- tf::tfd(matrix(muB,  nrow = 1), arg = subgrid)
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = subgrid)
  
  T_L2 <- n * tf::tf_integrate(diff_tfd^2, arg = subgrid)
  
  K  <- .covariance_estimator(X_sub, O_sub, group_A, muA, muB, est$pA, est$pB)
  KL <- .kl_decomposition(K, subgrid); lam <- KL$lam
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam)
  q <- which(cum >= fve)[1]
  q <- max(1L, q)
  lam_q <- lam[seq_len(q)]
  
  Z <- matrix(rnorm(q * n_sim), nrow = q)
  W <- colSums((Z^2) * lam_q)
  
  p <- (sum(W >= T_L2) + 1) / (length(W) + 1)
  
  if (isTRUE(compute_bands)) {
    bands <- .confidence_bands("L2", diff, W, n, alpha, subgrid, method = "asymptotic")
  } else {
    bands <- list(lower = NULL, upper = NULL, band = NULL)
  }
  
  data_name <- if (!is.null(fd)) "fd" else "X"
  
  if (isTRUE(bands_only)) {
    return(list(
      estimate = list(muA = muA_tfd,
                      muB = muB_tfd,
                      diff = diff_tfd),
      parameter = c(q = q, n_sim = n_sim),
      bands = bands
    ))
  }
  
  output <- .create_output("T_{mu,L2}", T_L2, p,
                           "L2 test",
                           data_name,
                           estimate = list(muA  = muA_tfd,
                                           muB  = muB_tfd,
                                           diff = diff_tfd),
                           parameter = c(q = q, n_sim = n_sim),
                      bands = bands)
  output
}

#' Asymptotic supremum test for MCAR (optional simultaneous bands)
#'
#' Tests \eqn{H_0:\ \mu_A=\mu_B} via \eqn{T_{\mu,D}=\sqrt{n}\,\lVert \hat\mu_A-\hat\mu_B\rVert_\infty}.
#' KL-based GP approximation; with `compute_bands=TRUE` also returns simultaneous bands.
#'
#' @inheritParams mcar_common-params
#' @param fve Fraction of variance explained (0-1) to select KL components.
#' @param n_sim Number of Monte Carlo draws for the GP approx.
#' @param seed RNG seed.
#' @param alpha Significance level for bands.
#' @param compute_bands Logical: also compute simultaneous bands?
#' @param bands_only Logical: return only band info?
#' @return `htest` (with extras) or a light band list.
#' @examples
#' set.seed(1)
#' m <- 50; n <- 40
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) { d <- diff(g)[1]; c(0, cumsum(rnorm(length(g)-1, sd = sqrt(d)))) }
#' X  <- t(replicate(n, bm(grid)))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X[O == 0L] <- NA_real_
#'
#' # Asymptotic Supremum test with simultaneous bands
#' res_sup <- asym_mean_sup_test(
#'   X = X,
#'   n_sim = 2000,
#'   compute_bands = TRUE,
#'   seed = 123
#' )
#' res_sup$p.value
#'
#' # Quick plot: mean difference + 95% simultaneous band
#' diff_hat <- tf::tf_evaluate(res_sup$estimate$diff, arg = res_sup$bands$grid)[[1]]
#' plot(res_sup$bands$grid, diff_hat, type = "l",
#'      ylim = range(c(res_sup$bands$lower, res_sup$bands$upper)),
#'      xlab = "t", ylab = "mean difference")
#' abline(h = 0, lty = 3)
#' lines(res_sup$bands$grid, res_sup$bands$lower, lty = 2)
#' lines(res_sup$bands$grid, res_sup$bands$upper, lty = 2)
#' @export
asym_mean_sup_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                               fve = 0.99, n_sim = 5000,
                               min_frac = 0.10, seed = NULL, alpha = 0.05,
                               compute_bands = TRUE, bands_only = FALSE) {
  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)
  
  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid <- grid[idx]
  X_sub  <- X[, idx, drop = FALSE]
  O_sub <- O[, idx, drop = FALSE]
  
  est <- .group_mean_estimators(X_sub, O_sub, group_A)
  muA <- est$muA
  muB <- est$muB
  diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA,  nrow = 1), arg = subgrid)
  muB_tfd  <- tf::tfd(matrix(muB,  nrow = 1), arg = subgrid)
  diff_tfd <- tf::tfd(matrix(diff, nrow = 1), arg = subgrid)
  
  T_D <- sqrt(n) * max(abs(diff))
  
  K  <- .covariance_estimator(X_sub, O_sub, group_A, muA, muB, est$pA, est$pB)
  KL <- .kl_decomposition(K, subgrid)
  lam <- KL$lam
  phi <- KL$phi
  
  if (!is.null(seed)) set.seed(seed)
  cum <- cumsum(lam) / sum(lam)
  q <- which(cum >= fve)[1]
  q <- max(1L, q)
  lam_q <- lam[seq_len(q)]
  phi_q <- phi[, seq_len(q), drop = FALSE]
  A <- sweep(phi_q, 2, sqrt(lam_q), "*")  
  
  Z <- matrix(rnorm(q * n_sim), nrow = q)     
  gp_vals <- A %*% Z                      
  W <- apply(abs(gp_vals), 2, max)
  
  p <- (sum(W >= T_D) + 1) / (length(W) + 1)
  
  if (isTRUE(compute_bands)) {
    bands <- .confidence_bands("D", diff, W, n, alpha, subgrid)
  } else {
    bands <- list(lower = NULL, upper = NULL, band = NULL)
  }
  
  data_name <- if (!is.null(fd)) "fd" else "X"
  
  if (isTRUE(bands_only)) {
    return(list(
      estimate = list(muA = muA_tfd,
                      muB = muB_tfd,
                      diff = diff_tfd),
      parameter = c(q = q, n_sim = n_sim),
      bands    = bands
    ))
  }
  
  output <- .create_output("T_{mu,D}", T_D, p,
                           "Supremum test",
                           data_name,
                           estimate = list(muA  = muA_tfd,
                                           muB  = muB_tfd,
                                           diff = diff_tfd),
                           parameter = c(q = q, n_sim=n_sim),
                           bands = bands)
  output
}

#' Bootstrap mean test (L2/Supremum) with optional bands
#'
#' Returns bootstrap p-values for L2 or Supremum and, optionally, confidence
#' bands for the difference curve. Parallelization via foreach/doParallel
#' with reproducible seeding via doRNG.
#'
#' @inheritParams mcar_common-params
#' @param n_boot Number of bootstrap iterations.
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
#' @param bands_only Return only the band payload?
#' @return `htest` (with extras) or a list of both tests.
#' @examples
#' set.seed(1)
#' m <- 30; n <- 200
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g){ d <- diff(g)[1]; c(0, cumsum(rnorm(length(g)-1, sd = sqrt(d)))) }
#'
#' # Group A: standard BM; Group B: BM with mean shift
#' X_A <- t(replicate(n/2, bm(grid)))
#' X_B <- t(replicate(n/2, bm(grid))) + 0.3
#' X <- rbind(X_A, X_B)
#'
#' # Define groups: FALSE = A, TRUE = B
#' groups <- c(rep(FALSE, n/2), rep(TRUE, n/2))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X <- X; X[O == 0L] <- NA_real_
#'
#' # Bootstrap Supremum test
#' res_boot <- boot_mean_test(
#'   X   = X,
#'   groups  = groups,
#'   n_boot  = 2000,   
#'   stat    = "D",       
#'   alpha   = 0.05,
#'   compute_bands = TRUE,
#'   parallel = FALSE,   
#'   seed    = 1
#' )
#'
#' res_boot$p.value 
#' @export
boot_mean_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                           n_boot = 5000,
                           min_frac = 0.10, alpha = 0.05,
                           parallel = TRUE,
                           ncpus = parallel::detectCores(logical = TRUE),
                           seed = NULL,
                           stat = c("L2", "D"),
                           compute_bands = TRUE,
                           return_boot = FALSE,
                           chunk_size = NULL,
                           manage_backend = "auto",
                           worker_blas_threads = 1L,
                           bands_only = FALSE) {
  
  # --- Argument checks ---
  checkmate::assert_character(stat, any.missing = FALSE, min.len = 1)
  checkmate::assert_subset(stat, c("L2","D"), empty.ok = FALSE)
  stat <- unique(stat)
  checkmate::assert_choice(manage_backend, c("auto","force_pool","sequential"))
  
  # --- Inputs vorbereiten ---
  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)
  
  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid   <- grid[idx]
  X_sub   <- X[, idx, drop = FALSE]
  O_sub   <- O[, idx, drop = FALSE]
  
  est  <- .group_mean_estimators(X_sub, O_sub, group_A)
  muA  <- est$muA
  muB <- est$muB
  diff <- muA - muB
  
  muA_tfd  <- tf::tfd(matrix(muA, 1), arg = subgrid)
  muB_tfd  <- tf::tfd(matrix(muB, 1), arg = subgrid)
  diff_tfd <- tf::tfd(matrix(diff, 1), arg = subgrid)
  
  w <- .trapezoid_weights(subgrid)
  
  # --- Teststatistiken ---
  T_L2 <- if ("L2" %in% stat) n * sum((diff^2) * w) else NULL
  T_D  <- if ("D"  %in% stat) sqrt(n) * max(abs(diff)) else NULL
  
  # --- Daten zentrieren ---
  IA <- as.numeric(group_A); IB <- 1 - IA
  X_cent <- X_sub
  if (any(IA == 1)) X_cent[IA == 1, ] <- sweep(X_sub[IA == 1,,drop=FALSE], 2, muA, `-`)
  if (any(IB == 1)) X_cent[IB == 1, ] <- sweep(X_sub[IB == 1,,drop=FALSE], 2, muB, `-`)
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- Bootstrap Runner ---
  .run_boot <- function(manage_backend_mode = manage_backend, ncpus_eff = ncpus) {
    be <- .init_parallel(manage_backend_mode, ncpus_eff, worker_blas_threads, seed)
    nworkers <- be$nworkers
    
    cs <- chunk_size
    if (is.null(cs) || !is.finite(cs) || cs < 1L) {
      cs <- max(50L, ceiling(n_boot / (3L * max(1L, nworkers))))
    } else cs <- as.integer(cs)
    
    idx_chunks <- split(seq_len(n_boot), ceiling(seq_len(n_boot) / cs))
    
    n_ <- n
    X_cent_ <- X_cent
    O_ <- O_sub
    group_A_ <- group_A
    w_ <- w
    
    do_L2 <- "L2" %in% stat; do_D <- "D" %in% stat
    
    boot_list <- foreach::foreach(ch = idx_chunks, .inorder = FALSE) %dorng% {
      out <- matrix(NA_real_, nrow = length(ch), ncol = sum(c(do_L2, do_D)) + 1L)
      cn  <- c(if (do_L2) "L2", if (do_D) "D", "redraws")
      colnames(out) <- cn
      diffs <- matrix(NA_real_, nrow = length(ch), ncol = ncol(X_cent_))
      
      for (ii in seq_along(ch)) {
        redraws <- 0L
        repeat {
          samp_idx <- sample.int(n_, n_, replace = TRUE)
          gA_samp <- group_A_[samp_idx]
          IA <- as.numeric(gA_samp)
          IB <- 1 - IA
          OA <- (O_[samp_idx,,drop=FALSE] * IA)
          OB <- (O_[samp_idx,,drop=FALSE] * IB)
          if (all(colSums(OA) > 0) && all(colSums(OB) > 0)) break
          redraws <- redraws + 1L
        }
        Xs <- X_cent_[samp_idx,,drop=FALSE]
        XA <- Xs
        XA[IB==1,] <- NA
        muA_b <- colSums(replace(XA, is.na(XA), 0)) / (n_ * colSums(OA)/n_)
        
        XB <- Xs
        XB[IA==1,] <- NA
        muB_b <- colSums(replace(XB, is.na(XB), 0)) / (n_ * colSums(OB)/n_)
        
        d_b <- muA_b - muB_b
        
        vals <- c()
        if (do_L2) vals <- c(vals, n_ * sum((d_b^2) * w_))
        if (do_D)  vals <- c(vals, sqrt(n_) * max(abs(d_b)))
        out[ii, ] <- c(vals, redraws)
        diffs[ii, ] <- d_b
      }
      list(stats = out, diffs = diffs)
    }
    
    boot_mat   <- do.call(rbind, lapply(boot_list, `[[`, "stats"))
    boot_diffs <- do.call(rbind, lapply(boot_list, `[[`, "diffs"))
    list(boot_mat = boot_mat, boot_diffs = boot_diffs)
  }
  
  boot_res <- tryCatch(
    .run_boot(manage_backend, ncpus),
    interrupt = function(e) { .reset_backend(); stop("Aborted by user: backend cleaned up; re-execution is immediately possible.", call.=FALSE) },
    error     = function(e) { .reset_backend(); stop(e) }
  )
  
  boot_mat   <- boot_res$boot_mat
  boot_diffs <- boot_res$boot_diffs
  
  outputs <- list()
  
  for (s in stat) {
    if (s == "L2") {
      boot_vals <- boot_mat[, "L2"]
      boot_vals <- boot_vals[!is.na(boot_vals)]
      T_val <- T_L2
      bands <- if (compute_bands) .confidence_bands("L2", diff, boot_diffs, n, alpha, subgrid, "bootstrap") else NULL
      method <- "Bootstrap mean test (L2)"
      stat_name <- "T_{mu,L2}"
    }
    if (s == "D") {
      boot_vals <- boot_mat[, "D"]
      boot_vals <- boot_vals[!is.na(boot_vals)]
      T_val <- T_D
      bands <- if (compute_bands) .confidence_bands("D", diff, boot_vals, n, alpha, subgrid, "bootstrap") else NULL
      method <- "Bootstrap mean test (supremum)"
      stat_name <- "T_{mu,D}"
    }
    
    p_val <- (sum(boot_vals >= T_val) + 1) / (length(boot_vals) + 1)
    
    outputs[[s]] <- .create_output(
      stat_name   = stat_name,
      stat_value  = T_val,
      p_value     = p_val,
      method      = method,
      data_name   = if (!is.null(fd)) "fd" else "X",
      estimate    = list(muA=muA_tfd, 
                         muB=muB_tfd, 
                         diff=diff_tfd),
      parameter   = c(n_boot=n_boot),
      bands       = bands
    )
  }
  
  if (length(outputs) == 1) return(outputs[[1]])
  outputs
}

