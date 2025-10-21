# ===============================================================
# Type-I-Error Simulation — tests whether the empirical Type I error
# approximates the nominal level of 0.05
# ===============================================================

suppressPackageStartupMessages({
  library(pbapply)
})

# ===============================================================
# Data generator 
# ===============================================================

# Generate one Brownian motion sample path
simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}

# Generate observed Brownian data under MCAR or MNAR mechanism
simulate_X_obs <- function(n, grid_spacing = 100,
                           mechanism = c("MCAR","MNAR")) {
  mechanism <- match.arg(mechanism)
  grid <- seq(0, 1, length.out = grid_spacing)
  m <- length(grid)
  
  # Simulate n Brownian motion trajectories
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  
  # Apply missingness mechanism (MCAR or MNAR)
  if (mechanism == "MCAR") {
    complete <- runif(n) < 0.5               
    O_mat <- matrix(0L, n, m)
    if (any(complete)) O_mat[complete, ] <- 1L
    if (any(!complete)) {
      k  <- sum(!complete)
      U  <- matrix(runif(k * 2), ncol = 2)
      L  <- pmin(U[,1], U[,2]); R <- pmax(U[,1], U[,2])
      Oi <- outer(L, grid, "<=") & outer(R, grid, ">") 
      O_mat[!complete, ] <- 1L * Oi
    }
  } else {  # MNAR: observed only when -1 < X(t) < 2
    within   <- (bm_mat > -1) & (bm_mat < 2) 
    complete <- rowSums(within) == m
    O_mat <- matrix(0L, n, m)
    if (any(complete))   O_mat[complete, ]   <- 1L
    if (any(!complete))  O_mat[!complete, ]  <- 1L * within[!complete, ]
  }
  
  # Apply observation mask to the data
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  X_obs
}

# ===============================================================
# One replication
# ===============================================================
one_rep_inner <- function(n, grid_spacing, mechanism, n_sim, n_boot, alpha) {
  
  X_obs <- simulate_X_obs(n, grid_spacing, mechanism)
  
  res_L2 <- asym_mean_L2_test(X=X_obs, n_sim=n_sim)
  res_D  <- asym_mean_sup_test(X=X_obs, n_sim=n_sim)
  res_BT_L2 <- boot_mean_test(X=X_obs, n_boot=n_boot, stat="L2")
  res_BT_D  <- boot_mean_test(X=X_obs, n_boot=n_boot, stat="D")
  
  # Return binary rejection indicators for each test
  c(
    asym_L2 = as.numeric(res_L2$p.value < alpha),
    asym_D  = as.numeric(res_D$p.value  < alpha),
    boot_L2 = as.numeric(res_BT_L2$p.value < alpha),
    boot_D  = as.numeric(res_BT_D$p.value  < alpha)
  )
}

# ===============================================================
# Multiple replications (runs = 5000)
# ===============================================================
run_type1_inner <- function(n, runs = 5000, grid_spacing = 100, mechanism = "MCAR",
                            n_sim = 10000, n_boot = 10000, alpha = 0.05) {
  # Display progress with timer
  pboptions(type = "timer")
  
  # Repeat 'one_rep_inner' for 'runs' iterations 
  res_arr <- pbreplicate(
    runs,
    one_rep_inner(n, grid_spacing, mechanism, n_sim, n_boot, alpha)
  )
  hits <- t(res_arr)
  colnames(hits) <- c("asym_L2","asym_D","boot_L2","boot_D")
  
  # Estimate rejection probabilities (Type-I error rates)
  colMeans(hits)
}

# ===============================================================
# Execute: n = 100, 250, 500
# ===============================================================

# ============================================================================  
# Precomputed Results  
#   - The files below contain precomputed Type I error estimates for n = 100, 250, 500.
# ============================================================================  
# res100 <- readRDS("Data/type_1_error_n100.rds")  
# res250 <- readRDS("Data/type_1_error_n250.rds")  
# res500 <- readRDS("Data/type_1_error_n500.rds")

set.seed(42)

res100 <- run_type1_inner(n = 100, runs = 100, n_sim = 1000, n_boot = 1000)

res250 <- run_type1_inner(n = 250, runs = 100, n_sim = 1000, n_boot = 1000)

res500 <- run_type1_inner(n = 500, runs = 100, n_sim = 1000, n_boot = 1000)

# Combine results for all sample sizes
type1_table <- rbind(
  data.frame(n = 100, t(res100)),
  data.frame(n = 250, t(res250)),
  data.frame(n = 500, t(res500))
)
print(type1_table, row.names = FALSE)
