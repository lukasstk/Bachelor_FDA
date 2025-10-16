#' Convert group labels to logical
#' @keywords internal
#' @noRd
.groups_to_logical <- function(groups) {
  checkmate::assert_atomic_vector(groups, any.missing = FALSE,
                                  .var.name = "groups")

  if (is.factor(groups)) groups <- as.character(droplevels(groups))

  if (length(unique(groups)) != 2) {
    stop("Assertion on 'groups' failed: Must contain exactly two distinct
         values (one for each group).", call. = FALSE)
  }

  if (is.logical(groups)) {
    return(groups)
  } else if (is.character(groups) || is.numeric(groups) || is.integer(groups)) {
    u <- unique(groups)
    return(groups == u[1])
  }

}

#' Coerce `tfd` to dense matrix + grid
#' @keywords internal
#' @noRd
.fd_to_matrix <- function(fd) {
  checkmate::assert_class(fd, c("tfd"))

  all_grids <- tf::tf_arg(fd)

  # For irregular data ("tfd_irreg"), tf_arg(fd) returns a list of grids
  # (one per curve); for regular data ("tfd_reg"), a single numeric grid vector.
  if (is.list(all_grids)) {
    g <- sort(unique(unlist(all_grids)))
  } else {
    g <- all_grids
  }

  X_try <- suppressWarnings(as.matrix(fd))

  list(X = X_try, grid = g)
}

#' Subdomain selector
#'
#' Selects grid points where both groups have sufficient coverage.
#' Each retained time point satisfies:
#' min(cA(t), cB(t)) > n * min_frac
#' If no grid point satisfies this, an error is raised.
#'
#' @keywords internal
#' @noRd
.limit_subdomain <- function(O, group_A, min_frac = 0.10) {
  checkmate::assert_matrix(O, any.missing = FALSE)
  checkmate::assert_logical(group_A, len = nrow(O))
  checkmate::assert_number(min_frac, lower = 0, upper = 1)

  n <- nrow(O)
  IA <- as.numeric(group_A)
  IB <- 1 - IA

  # Observed counts per grid point in each group
  cA <- colSums(O * IA)
  cB <- colSums(O * IB)

  # Select points where both groups exceed threshold
  valid_points <- pmin(cA, cB) > n * min_frac

  # If no grid point satisfies the criterion stop the process
  if (!any(valid_points)) {
    stop(sprintf(
      "No suitable subdomain found: both groups must exceed
      %.0f%% coverage (min_frac = %.2f) at \u2265 1 grid point.",
      100 * min_frac, min_frac
    ), call. = FALSE)
  }

  # Return indices of valid grid points
  idx <- which(valid_points)

  list(idx = idx, min_frac_used = min_frac)
}



#' Prepare inputs from `tfd` or matrix
#' @keywords internal
#' @noRd
.prepare_inputs <- function(fd = NULL, X = NULL,
                            groups = NULL, observed_ratio = 1) {
  checkmate::assert_number(observed_ratio, lower = 0, upper = 1)

  if (!is.null(fd)) {
    conv <- .fd_to_matrix(fd)
    X <- conv$X
    grid <- conv$grid
  } else {
    checkmate::assert_matrix(
      X,
      mode = "numeric",
      any.missing = TRUE,
      all.missing = FALSE,
      null.ok = FALSE
    )
    grid <- seq(0, 1, length.out = ncol(X))
  }

  n <- nrow(X)
  O <- 1L * (!is.na(X))

  if (!is.null(groups)) {
    checkmate::assert_atomic_vector(groups, len = n, any.missing = FALSE)
    group_A <- .groups_to_logical(groups)

    obs_frac <- rowMeans(O != 0)
    meanA <- mean(obs_frac[group_A], na.rm = TRUE)
    meanB <- mean(obs_frac[!group_A], na.rm = TRUE)
    if (meanA < meanB) {
      group_A <- !group_A
      message(sprintf(
        "Note: swapped labels so Group A is the more complete group
        (mean A=%.3f, B=%.3f).",
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
    X = X,
    O = O,
    group_A = as.logical(group_A),
    grid = grid
  )
}
