#' Bootstrap mean test (L2/Supremum)
#'
#' Returns bootstrap p-values for L2 or Supremum.
#'
#' @inheritParams mcar_common-params
#' @param n_boot Target number of bootstrap iterations.
#' @param max_redraws Maximum number of redraw attempts for invalid bootstrap samples.
#' @param alpha Significance level for bands (Supremum only).
#' @param ncpus Number of workers for an internal cluster (if created).
#' @param seed RNG seed (passed to doRNG).
#' @param stat `"L2"`, `"D"` oder `c("L2","D")`.
#' @param compute_bands Compute confidence bands (Supremum only)?
#' @param chunk_size Number of bootstrap replicates per foreach task.
#' @param manage_backend Backend control (`"auto"`, `"force_pool"`, `"sequential"`).
#' @param worker_blas_threads BLAS/OpenMP threads per worker (internal pool only).
#' @param bands_only Logical: if TRUE, return only band information instead of a full `htest`.
#' @return `htest` (with extras) or a list of both tests.
#' @examples
#' set.seed(1)
#' m <- 30
#' n <- 200
#' grid <- seq(0, 1, length.out = m)
#' bm <- function(g) {
#'   d <- diff(g)[1]
#'   c(0, cumsum(rnorm(length(g) - 1, sd = sqrt(d))))
#' }
#'
#' # Group A: standard BM; Group B: BM with mean shift
#' X_A <- t(replicate(n / 2, bm(grid)))
#' X_B <- t(replicate(n / 2, bm(grid))) + 0.3
#' X <- rbind(X_A, X_B)
#'
#' # Define groups: FALSE = A, TRUE = B
#' groups <- c(rep(FALSE, n / 2), rep(TRUE, n / 2))
#'
#' # MNAR censoring: observed only if -1 < X(t) < 2
#' O <- 1L * (X > -1 & X < 2)
#' X[O == 0L] <- NA_real_
#'
#' # Bootstrap Supremum test
#' res_boot <- boot_mean_test(
#'   X = X,
#'   groups = groups,
#'   n_boot = 2000,
#'   stat = "D",
#'   alpha = 0.05,
#'   compute_bands = TRUE,
#'   manage_backend = "sequential",
#'   seed = 1
#' )
#'
#' res_boot$p.value
#' @export
boot_mean_test <- function(fd = NULL, X = NULL, groups = NULL, observed_ratio = 1,
                           n_boot = 10000,
                           min_frac = 0.10, alpha = 0.05,
                           ncpus = parallel::detectCores(logical = TRUE),
                           seed = NULL,
                           stat = c("L2", "D"),
                           compute_bands = TRUE,
                           chunk_size = NULL,
                           manage_backend = "auto",
                           worker_blas_threads = 1L,
                           bands_only = FALSE,
                           max_redraws = 20L) {
  checkmate::assert_character(stat, any.missing = FALSE, min.len = 1)
  checkmate::assert_int(chunk_size, lower = 1, null.ok = TRUE)
  checkmate::assert_subset(stat, c("L2", "D"), empty.ok = FALSE)
  stat <- unique(stat)
  checkmate::assert_choice(manage_backend, c("auto", "force_pool", "sequential"))

  prep <- .prepare_inputs(fd, X, groups, observed_ratio)
  X <- prep$X
  O <- prep$O
  group_A <- prep$group_A
  grid <- prep$grid
  n <- nrow(X)

  subdomain <- .limit_subdomain(O, group_A, min_frac = min_frac)
  idx <- subdomain$idx
  subgrid <- grid[idx]
  X_sub <- X[, idx, drop = FALSE]
  O_sub <- O[, idx, drop = FALSE]

  est <- .group_mean_estimators(X_sub, O_sub, group_A)
  mean_A <- est$mean_A
  mean_B <- est$mean_B
  mean_diff <- mean_A - mean_B

  mean_A_tfd <- tf::tfd(matrix(mean_A, 1), arg = subgrid)
  mean_B_tfd <- tf::tfd(matrix(mean_B, 1), arg = subgrid)
  mean_diff_tfd <- tf::tfd(matrix(mean_diff, 1), arg = subgrid)

  w <- .trapezoid_weights(subgrid)

  T_L2 <- if ("L2" %in% stat) n * sum((mean_diff^2) * w) else NULL
  T_D <- if ("D" %in% stat) sqrt(n) * max(abs(mean_diff)) else NULL
  T_vals <- list(L2 = T_L2, D = T_D)

  IA <- as.numeric(group_A)
  IB <- 1 - IA
  X_cent <- X_sub
  X_cent[IA == 1, ] <- sweep(X_sub[IA == 1, , drop = FALSE], 2, mean_A, `-`)
  X_cent[IB == 1, ] <- sweep(X_sub[IB == 1, , drop = FALSE], 2, mean_B, `-`)

  if (!is.null(seed)) set.seed(seed)

  .run_boot <- function(manage_backend_mode = manage_backend, ncpus_eff = ncpus) {
    be <- .init_parallel(manage_backend_mode, ncpus_eff, worker_blas_threads, seed)
    nworkers <- be$nworkers

    cs <- chunk_size
    if (is.null(cs)) {
      cs <- max(50L, ceiling(n_boot / (3L * nworkers)))
    } else {
      cs <- as.integer(cs)
    }

    idx_chunks <- split(seq_len(n_boot), ceiling(seq_len(n_boot) / cs))
    chunk_idx <- NULL

    boot_list <- foreach::foreach(
      chunk_idx = idx_chunks,
      .inorder = FALSE,
      .export = c(
        "stat", "n", "w", "group_A", "O_sub", "X_cent",
        "max_redraws"
      )
    ) %dorng% {
      res <- matrix(NA_real_, nrow = length(chunk_idx), ncol = length(stat))
      colnames(res) <- stat
      diffs <- matrix(NA_real_, nrow = length(chunk_idx), ncol = ncol(X_cent))

      for (i in seq_along(chunk_idx)) {
        redraws <- 0
        draw_successful <- FALSE

        while (!draw_successful && redraws <= max_redraws) {
          samp <- sample.int(n, n, replace = TRUE)
          gA <- group_A[samp]
          IA <- as.numeric(gA)
          IB <- 1 - IA

          OA <- O_sub[samp, , drop = FALSE] * IA
          OB <- O_sub[samp, , drop = FALSE] * IB

          if (any(colSums(OA) == 0) || any(colSums(OB) == 0)) {
            redraws <- redraws + 1
            next
          }

          Xs <- X_cent[samp, , drop = FALSE]
          denomA <- colSums(OA)
          denomB <- colSums(OB)

          mean_A_boot <- colSums(Xs * OA * IA, na.rm = TRUE) / denomA
          mean_B_boot <- colSums(Xs * OB * IB, na.rm = TRUE) / denomB

          diff_boot <- mean_A_boot - mean_B_boot
          diffs[i, ] <- diff_boot

          if ("L2" %in% stat) {
            res[i, "L2"] <- n * sum((diff_boot^2) * w)
          }
          if ("D" %in% stat) {
            res[i, "D"] <- sqrt(n) * max(abs(diff_boot))
          }

          draw_successful <- TRUE
        }

        if (!draw_successful) {
          warning(sprintf(
            "Max redraws (%d) exceeded in one bootstrap replicate, skipping.",
            max_redraws
          ))
        }
      }
      list(stats = res, diffs = diffs)
    }

    boot_mat <- do.call(rbind, lapply(boot_list, `[[`, "stats"))
    boot_diffs <- do.call(rbind, lapply(boot_list, `[[`, "diffs"))

    list(boot_mat = boot_mat, boot_diffs = boot_diffs)
  }

  boot_res <- tryCatch(
    .run_boot(manage_backend, ncpus),
    interrupt = function(e) {
      .reset_backend()
      stop("Aborted by user: backend cleaned up; re-execution is immediately possible.", call. = FALSE)
    },
    error = function(e) {
      .reset_backend()
      stop(e)
    }
  )

  boot_mat <- boot_res$boot_mat
  boot_diffs <- boot_res$boot_diffs

  valid_idx <- rowSums(!is.na(boot_mat)) > 0
  n_valid <- sum(valid_idx)

  n_required <- ceiling(500 / alpha)

  if (n_valid < 0.8 * n_required) {
    warning(sprintf(
      "Only %d valid bootstrap replicates (< 80%% of required %d for alpha=%.3f). Results may not provide enough samples to reliably estimate distributional quantiles.",
      n_valid, n_required, alpha
    ), call. = FALSE)
  }

  outputs <- lapply(stat, function(s) {
    boot_vals <- boot_mat[valid_idx, s]
    T_val <- T_vals[[s]]
    bands <- NULL
    if (s == "D" && compute_bands) {
      bands <- .confidence_bands("D", mean_diff, boot_vals, n, alpha, subgrid)
    }

    if (isTRUE(bands_only)) {
      return(list(
        estimate = list(
          mean_A = tf::tfd(matrix(est$mean_A, 1), arg = subgrid),
          mean_B = tf::tfd(matrix(est$mean_B, 1), arg = subgrid),
          mean_diff = mean_diff_tfd
        ),
        parameter = c(n_boot = n_valid),
        bands = bands
      ))
    }

    .create_output(
      stat_name = switch(s,
        L2 = "T_\u03bc,L\u00B2",
        D = "T_\u03bc,D"
      ),
      stat_value = T_val,
      p_value = (sum(boot_vals >= T_val) + 1) / (length(boot_vals) + 1),
      method = switch(s,
        L2 = "Bootstrap mean test (L\u00B2)",
        D = "Bootstrap mean test (Supremum)"
      ),
      data_name = if (!is.null(fd)) "fd" else "X",
      estimate = list(
        mean_A = tf::tfd(matrix(est$mean_A, 1), arg = subgrid),
        mean_B = tf::tfd(matrix(est$mean_B, 1), arg = subgrid),
        mean_diff = mean_diff_tfd
      ),
      parameter = c(n_boot = n_valid),
      bands = bands
    )
  })

  if (length(outputs) == 1) {
    return(outputs[[1]])
  }
  outputs
}
