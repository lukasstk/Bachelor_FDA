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
.fd_to_matrix <- function(fd) {
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

  n <- nrow(O)
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
    conv <- .fd_to_matrix(fd)
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
    X = X,
    O = O,
    group_A = as.logical(group_A),
    grid = grid_vec
  )
}
