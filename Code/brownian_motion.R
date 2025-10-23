# ===============================================================
# Brownian-Motion simulation study 
#   - Calculation of rejection probabilities
#   - visualization
# ===============================================================
suppressPackageStartupMessages({
  library(extrafont)
  library(patchwork)
  library(tidyverse)
  library(pbapply)
})

set.seed(42)

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
                           mechanism = c("MCAR", "MNAR")) {
  mechanism <- match.arg(mechanism)
  grid <- seq(0, 1, length.out = grid_spacing)
  m <- length(grid)
  
  # Simulate n Brownian motion trajectories
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  rownames(bm_mat) <- paste0("ID", seq_len(n))
  
  # Apply missingness mechanism
  if (mechanism == "MCAR") {
    complete <- runif(n) < 0.5
    O_mat <- matrix(0L, n, m)
    if (any(complete)) O_mat[complete, ] <- 1L
    if (any(!complete)) {
      k  <- sum(!complete)
      U  <- matrix(runif(k * 2), ncol = 2)
      L  <- pmin(U[, 1], U[, 2]); R <- pmax(U[, 1], U[, 2])
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
  
  # Apply observation mask to data
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  list(X = X_obs, grid = grid)
}

# ===============================================================
# Simulate data and visualize
# ===============================================================

# Generate observed Brownian data
sim_data <- simulate_X_obs(n = 100, grid_spacing = 100, mechanism = "MNAR")
X   <- sim_data$X
grid <- sim_data$grid

# Prepare data in long format
df_obs <- as.data.frame(X) %>%
  mutate(
    id = rownames(X),
    group = ifelse(rowSums(is.na(.)) == 0, "complete", "incomplete")
  ) %>%
  pivot_longer(cols = -c(id, group), names_to = "t_idx", values_to = "x") %>%
  mutate(t = grid[as.integer(sub("V", "", t_idx))])

p_left <- ggplot(df_obs, aes(x = t, y = x, group = id)) +
  geom_line(
    data = subset(df_obs, group == "complete"),
    color = "grey70",
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_line(
    data = subset(df_obs, group == "incomplete"),
    color = "black",
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  labs(x = "time", y = "X(t)") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    legend.position = "none",
    text            = element_text(family = "Times New Roman"),
    axis.title      = element_text(size = 20),
    axis.title.x    = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y    = element_text(size = 20, margin = margin(r = 10)),
    axis.text       = element_text(size = 16),
    axis.line.x     = element_line(color = "black", linewidth = 0.5),
    axis.line.y     = element_line(color = "black", linewidth = 0.5),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt")
  )

# --- Run asymptotic supremum test and extract mean difference with bands ---
res_sup <- asym_mean_sup_test(
  X = X,
  seed = 42,
  bands_only = TRUE
)

# --- Prepare data frame for plotting mean difference and simultaneous bands ---
bands <- data.frame(
  t     = res_sup$bands$grid,
  diff  = tf::tf_evaluate(res_sup$estimate$mean_diff, arg = res_sup$bands$grid)[[1]],
  lower = res_sup$bands$lower,
  upper = res_sup$bands$upper
)

p_right <- ggplot(bands, aes(t, diff)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "grey50") +
  geom_line(linewidth = 1, color = "black") +
  geom_line(aes(y = lower), linetype = 2, color = "black") +
  geom_line(aes(y = upper), linetype = 2, color = "black") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.5),
                     breaks = seq(-1, 1.5, by = 0.5)) +
  labs(x = "time", y = "difference in means") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    text            = element_text(family = "Times New Roman"),
    axis.title      = element_text(size = 18),
    axis.title.x    = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y    = element_text(size = 20, margin = margin(r = 10)),
    axis.text       = element_text(size = 16),
    axis.line.x     = element_line(color = "black", linewidth = 0.5),
    axis.line.y     = element_line(color = "black", linewidth = 0.5),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt")
  )

brownian_motion_plot <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# # Abspeichern
# ggsave(
#   filename = "Plots/brownian_motion_plot.png",
#   plot     = brownian_motion_plot,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )


# ============================================================================
# Plot: Rejection probabilities vs. b (MNAR, a = -1 < b)
# Tests the power of T_{μ,L²} and T_{μ,D} under increasing truncation level b
# ============================================================================

# ---- Simulation settings ----
n        <- 100
grid     <- seq(0, 1, length.out = 100)
alpha    <- 0.05
b_vals   <- seq(1.0, 2.0, by = 0.1)
n_sims   <- 5000

# ---- Helper: MNAR censoring mechanism ----
make_missing_pattern <- function(x_row, a = -1, b = 2) {
  as.integer(x_row > a & x_row < b)
}

# ---- One Monte Carlo run for given b ----
one_run <- function(b_now) {
  bm_mat <- t(replicate(n, simulate_bm(grid)))                 # simulate n Brownian paths
  O  <- t(apply(bm_mat, 1L, make_missing_pattern, a = -1, b = b_now))  # MNAR observation mask
  storage.mode(O) <- "integer"
  
  X <- bm_mat
  X[O == 0L] <- NA_real_                                       # apply missingness
  
  p_L2 <- tryCatch(                                             # asymptotic L2 test
    asym_mean_L2_test(X = X)$p.value,
    error = function(e) NA_real_
  )
  p_D <- tryCatch(                                              # asymptotic supremum test
    asym_mean_sup_test(X = X, compute_bands = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  c(p_L2, p_D)
}

# ---- Monte Carlo loop over b values ----
res_list <- pblapply(b_vals, function(b_now) {
  ps  <- pbreplicate(n_sims, one_run(b_now))                   # replicate simulations
  rej <- rowMeans(ps < alpha, na.rm = TRUE)                    # rejection rates
  
  cat(sprintf("b = %.1f  ->  rej L2 = %.3f,  rej D = %.3f\n",
              b_now, rej[1], rej[2]))
  
  data.frame(
    b    = b_now,
    test = c("T[mu,L^2]", "T[mu,D]"),
    rej  = rej
  )
})


# ============================================================================  
# Precomputed Results:
#   - The file below contains rejection probabilities from 5000 Monte Carlo runs.  
# ============================================================================  
# res_df <- readRDS("Data/brownian_motion_rej_probs_df.rds")

res_df <- dplyr::bind_rows(res_list)


brownian_motion_rej_probs_plot <- ggplot(res_df, aes(x = b, y = rej, shape = test, color = test)) +
  geom_point(size = 5, stroke = 1.5, shape = 4) +
  scale_color_manual(
    values = c("T[mu,L^2]" = "steelblue",
               "T[mu,D]"   = "red"),
    breaks = c("T[mu,L^2]", "T[mu,D]"),
    labels = c(
      "T[mu,L^2]" = expression(T[mu]*","~L^2),
      "T[mu,D]"   = expression(T[mu]*","~D)
    )
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(1, 2), breaks = seq(1, 2, by = 0.2)) +
  labs(x = "b", y = "rejection probability", shape = NULL, color = NULL) +
  theme_classic(base_size = 18, base_family = "Times New Roman") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.895, 0.227),
    legend.margin = margin(7, 7, 7, 7),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 25, margin = margin(t = 30)), 
    axis.title.y = element_text(size = 25, margin = margin(r = 30)),  
    axis.text  = element_text(size = 18),                             
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key.size = unit(1.2, "lines"),             
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.major = element_line(colour = "grey80", size = 0.5),
    panel.grid.minor = element_blank()
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4))
  )

# # Save figure
# ggsave(
#   filename = "Plots/brownian_motion_rej_probs_plot.png",
#   plot     = brownian_motion_rej_probs_plot,
#   width    = 10,    
#   height   = 5,
#   dpi      = 300
# )
