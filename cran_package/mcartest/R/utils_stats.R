#' Trapezoidal integration weights
#' @keywords internal
#' @noRd
.trapezoid_weights <- function(x) {
  checkmate::assert_numeric(x, any.missing = FALSE, min.len = 1, sorted = TRUE)

  m <- length(x)
  if (m == 1L) {
    return(1)
  }
  w <- numeric(m)
  w[1] <- (x[2] - x[1]) / 2
  w[m] <- (x[m] - x[m - 1]) / 2
  if (m > 2L) w[2:(m - 1)] <- (x[3:m] - x[1:(m - 2)]) / 2
  w
}


#' Available-mean estimators by group
#' @keywords internal
#' @noRd
.group_mean_estimators <- function(X, O, group_A) {
  n <- nrow(X)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  pA_hat <- colSums(O * IA) / n
  pB_hat <- colSums(O * IB) / n
  mean_A_hat <- colSums(X * IA, na.rm = TRUE) / (n * pA_hat)
  mean_B_hat <- colSums(X * IB, na.rm = TRUE) / (n * pB_hat)
  list(mean_A = mean_A_hat, mean_B = mean_B_hat, pA = pA_hat, pB = pB_hat)
}

#' Corrected covariance under partial observation
#' @keywords internal
#' @noRd
.covariance_estimator <- function(X, O, group_A, muA_hat, muB_hat, pA_hat, pB_hat) {
  n <- nrow(X)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  Xtilde <- X
  Xtilde[group_A, ] <- sweep(X[group_A, , drop = FALSE], 2, muA_hat, `-`)
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
  w <- .trapezoid_weights(grid)
  sw <- sqrt(w)
  S <- (sw %o% sw) * K
  S <- (S + t(S)) / 2
  ev_test <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev_test < -tol * max(1, abs(ev_test[1])))) {
    stop(sprintf("Weighted covariance is not PSD (min eigenvalue = %.4g).", min(ev_test)), call. = FALSE)
  }
  ev <- eigen(S, symmetric = TRUE)
  eigenvalues <- ev$values
  eigenfunctions <- sweep(ev$vectors, 1, sw, "/")
  norms <- sqrt(colSums(eigenfunctions^2 * w))
  eigenfunctions <- sweep(eigenfunctions, 2, norms, "/")
  list(eigenvalues = eigenvalues, eigenfunctions = eigenfunctions, w = w)
}

#' Simultaneous confidence bands
#'
#' @keywords internal
#' @noRd
.confidence_bands <- function(stat, diff, W, n, alpha, grid) {
  if (identical(stat, "D")) {
    q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE, na.rm = TRUE))
    halfwidth <- q_alpha / sqrt(n)
    lower <- diff - halfwidth
    upper <- diff + halfwidth
    band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
    return(list(
      type = "simultaneous", band = band,
      lower = lower, upper = upper, alpha = alpha, grid = grid
    ))
  }

  stop("Unknown stat type in .confidence_bands(): must be 'D'.")
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
