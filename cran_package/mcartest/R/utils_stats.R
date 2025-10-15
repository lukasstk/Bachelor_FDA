#' Trapezoidal integration weights
#' @keywords internal
#' @noRd
.trapezoid_weights <- function(grid) {
  checkmate::assert_numeric(grid, any.missing = FALSE, min.len = 1, sorted = TRUE)
  
  n_points <- length(grid)
  if (n_points == 1L) {
    return(1)
  }
  
  # initialize weights
  w <- numeric(n_points)
  
  # boundary weights
  w[1] <- (grid[2] - grid[1]) / 2
  w[n_points] <- (grid[n_points] - grid[n_points - 1]) / 2
  
  # interior weights
  if (n_points > 2L)
    w[2:(n_points - 1)] <- (grid[3:n_points] - grid[1:(n_points - 2)]) / 2
  
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
  mean_A_hat <- colSums(X * O * IA, na.rm = TRUE) / (n * pA_hat)
  mean_B_hat <- colSums(X * O * IB, na.rm = TRUE) / (n * pB_hat)
  list(mean_A = mean_A_hat, mean_B = mean_B_hat, pA = pA_hat, pB = pB_hat)
}

#' Corrected covariance under partial observation
#' @keywords internal
#' @noRd
.covariance_estimator <- function(X, O, group_A, mean_A, mean_B, pA, pB) {
  n <- nrow(X)
  IA <- as.numeric(group_A)
  IB <- 1 - IA
  X_centered <- X
  X_centered[group_A, ] <- sweep(X[group_A, , drop = FALSE], 2, mean_A, `-`)
  X_centered[!group_A, ] <- sweep(X[!group_A, , drop = FALSE], 2, mean_B, `-`)
  X_centered[is.na(X_centered)] <- 0
  A_resid <- sweep((X_centered * O) * IA, 2, pA, "/")
  B_resid <- sweep((X_centered * O) * IB, 2, pB, "/")
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
  ev <- eigen(S, symmetric = TRUE)
  if (any(ev$values < -tol * max(1, abs(ev$values[1])))) {
    stop(sprintf("Weighted covariance not PSD (min eigenvalue = %.4g).",
                 min(ev$values)), call. = FALSE)
  }
  eigenvalues <- ev$values
  eigenvalues[eigenvalues < 0] <- 0
  eigenfunctions <- sweep(ev$vectors, 1, sw, "/")
  list(eigenvalues = eigenvalues, eigenfunctions = eigenfunctions, w = w)
}

#' Simultaneous confidence bands
#'
#' @keywords internal
#' @noRd
.confidence_bands <- function(stat, diff, W, n, alpha, grid) {
  checkmate::assert_choice(stat, choices = "D", .var.name = "stat")
  q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE, na.rm = TRUE))
  halfwidth <- q_alpha / sqrt(n)
  lower <- diff - halfwidth
  upper <- diff + halfwidth
  band <- tf::tfd(matrix(c(lower, upper), nrow = 2, byrow = TRUE), arg = grid)
  
  list(
    type = "simultaneous",
    band = band,
    lower = lower,
    upper = upper,
    alpha = alpha,
    grid = grid
  )

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
