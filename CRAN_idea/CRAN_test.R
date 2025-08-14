# ===============================================================
# Brownian-Motion Testdaten + Anwendung der Algorithmen (1–3,5)
# ===============================================================

library(dplyr)
library(tidyr)
library(tidyfun)
library(tf)
library(ggplot2)
library(patchwork)

set.seed(42)

# ---------------- Parameter ----------------
n            <- 100                    # Anzahl Pfade
grid_spacing <- 100
grid         <- seq(0, 1, length.out = grid_spacing)
ids          <- paste0("ID", seq_len(n))
mechanism    <- "MNAR"                 # "MCAR" oder "MNAR"

# ---------------- Hilfsfunktionen ----------------
simulate_bm <- function(grid, mean_shift = 0) {
  dt <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc)) + mean_shift
}

make_O_mcar <- function(grid) {
  if (runif(1) < 0.5) {
    rep(1, length(grid))  # vollständig
  } else {
    U1 <- runif(1); U2 <- runif(1)
    L  <- min(U1, U2); U <- max(U1, U2)
    as.numeric(grid >= L & grid < U)  # nur im Intervall beobachtet
  }
}

make_O_mnar <- function(x_row) {
  if (all(x_row > -1 & x_row < 2)) {
    rep(1, length(x_row))                  # vollständig beobachtet
  } else {
    as.numeric(x_row > -1 & x_row < 2)     # sonst nur im Bereich beobachtet
  }
}

# ---------------- Daten simulieren ----------------
# Brownian Paths (unter H0 ohne Mean-Shift)
bm_mat <- t(replicate(n, simulate_bm(grid)))
rownames(bm_mat) <- ids

# --- Beobachtungsmatrix rein aus dem Mechanismus erzeugen ---
if (mechanism == "MCAR") {
  O_mat <- t(sapply(seq_len(n), function(i) make_O_mcar(grid)))
} else if (mechanism == "MNAR") {
  O_mat <- t(apply(bm_mat, 1L, make_O_mnar))
} else {
  stop("Unknown mechanism")
}
storage.mode(O_mat) <- "integer"
rownames(O_mat) <- ids
colnames(O_mat) <- NULL

# --- Gruppen erst NACH der Zensur definieren ---
is_complete <- rowSums(O_mat) == length(grid)
A_ids <- ids[is_complete]
B_ids <- ids[!is_complete]

# Beobachtete Werte (NA außerhalb Beobachtung)
X_obs <- bm_mat
X_obs[O_mat == 0L] <- NA_real_

# ---------------- Long-Format bauen (df_long) ----------------
# Spaltennamen als echte t-Werte, damit pivot_longer sauber wird
colnames(X_obs) <- formatC(grid, digits = 12, format = "fg", flag = "#")

df_long <- as.data.frame(X_obs) |>
  mutate(id = ids) |>
  tidyr::pivot_longer(
    cols = -id,
    names_to = "t",
    values_to = "x"
  ) |>
  mutate(t = as.numeric(t)) |>
  arrange(id, t)

# Gruppen-Dataframe
df_groups <- data.frame(
  id    = ids,
  group = ifelse(ids %in% A_ids, "A", "B"),
  stringsAsFactors = FALSE
)

# ---------------- Pipeline (unverändert) ----------------
# 1) Grid & IDs
grid <- sort(unique(df_long$t))
ids  <- unique(df_long$id)

# 2) Wide-Matrix X_obs (n x m)
df_wide <- df_long %>%
  select(id, t, x) %>%
  mutate(t = as.numeric(t)) %>%
  tidyr::pivot_wider(names_from = t, values_from = x) %>%
  arrange(match(id, ids))

ord   <- order(as.numeric(names(df_wide)[-1]))
X_obs <- as.matrix(df_wide[, c(1 + ord)])
rownames(X_obs) <- df_wide$id
colnames(X_obs) <- NULL

# 3) Beobachtungsmatrix O_mat aus NA-Maske
O_mat <- 1L * !is.na(X_obs)

# 4) Gruppencode 0/1
group_A <- df_groups %>%
  arrange(match(id, ids)) %>%
  transmute(is_A = as.integer(group == "A")) %>%
  pull(is_A)

# ---------------- Algorithmen aufrufen ----------------
B_mc <- 2000
alpha <- 0.05

# Algo 1: L2-Test
res_L2 <- tfu_algo1_L2_test(
  X_obs = X_obs
)
cat("Algo1 L2  -> stat:", res_L2$stat, " p:", res_L2$p_value, "\n")

# Algo 2: Supremums-Test
res_D <- tfu_algo2_sup_test(
  X_obs = X_obs
)
cat("Algo2 Sup -> stat:", res_D$stat, " p:", res_D$p_value, "\n")

# Algo 3: Simultane Konfidenzbänder
res_CB <- tfu_algo3_conf_bands(
  X_obs = X_obs
)
bands <- data.frame(
  t = grid[res_CB$idx],
  diff = res_CB$diff,
  lower = res_CB$lower,
  upper = res_CB$upper
)

# Algo 5: Bootstrap-P-Werte (L2 & Supremum)
res_boot <- tfu_algo5_bootstrap(
  X_obs = X_obs
)
cat("Algo5 Boot -> L2: stat", res_boot$T_L2, " p", res_boot$p_L2,
    " | D: stat", res_boot$T_D, " p", res_boot$p_D, "\n")

# ---------------- Plot: links Kurven (unvollständig=schwarz), rechts Diff + 95%-Band ----------------
df_obs <- as.data.frame(X_obs) %>%
  mutate(id = rownames(X_obs),
         gruppe = ifelse(rowSums(is.na(.)) == 0, "vollständig", "unvollständig")) %>%
  pivot_longer(cols = -c(id, gruppe), names_to = "t_idx", values_to = "x") %>%
  mutate(t = grid[as.integer(sub("V", "", t_idx))])

p_left <- ggplot(df_obs, aes(x = t, y = x, group = id)) +
  geom_line(aes(color = gruppe), linewidth = 0.9, na.rm = TRUE) +
  scale_color_manual(values = c("vollständig" = "grey60", "unvollständig" = "black")) +
  labs(title = "Brownian-Pfade nach MNAR-Zensur", x = "Zeit t", y = "X(t)", color = "Gruppe") +
  theme_minimal() + theme(legend.position = "bottom")

p_right <- ggplot(bands, aes(t, diff)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = lower), linetype = 2) +
  geom_line(aes(y = upper), linetype = 2) +
  labs(x = "time", y = "difference in means") +
  theme_minimal()

p_left + p_right

# Plot:Rejection Probs
# ===========================================
# Rejection curves vs. b (MNAR, a = -1 < b)
# Only T_{mu,L2} and T_{mu,D}
# ===========================================

library(dplyr)
library(ggplot2)

set.seed(2025)

# ---- Simulation settings ----
n        <- 100                 # sample size per run (wie im Paper)
m_grid   <- 100                 # grid points
grid     <- seq(0, 1, length.out = m_grid)
alpha    <- 0.05
b_vals   <- seq(1.0, 2.0, by = 0.1)   # vary b
n_sims   <- 1000                       # für Speed; Paper ~5000

# ---- Helpers ----
simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}

# MNAR censoring: observe only when a < X(t) < b  (a = -1)
make_O_mnar_ab <- function(x_row, a = -1, b = 2) {
  as.integer(x_row > a & x_row < b)
}

one_run <- function(b_now) {
  # 1) simulate n Brownian paths on grid
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  # 2) build observation mask under MNAR censoring
  O_mat  <- t(apply(bm_mat, 1L, make_O_mnar_ab, a = -1, b = b_now))
  storage.mode(O_mat) <- "integer"
  
  # 3) observed matrix with NAs
  X_obs <- bm_mat
  X_obs[O_mat == 0L] <- NA_real_
  
  # 4) groups: complete vs incomplete AFTER censoring
  complete   <- rowSums(O_mat) == ncol(O_mat)
  nA <- sum(complete); nB <- sum(!complete)
  if (nA == 0 || nB == 0) return(c(NA_real_, NA_real_))  # skip run if only one group
  
  # 5) run tests (asymptotic)
  p_L2 <- tryCatch(tfu_algo1_L2_test(X_obs)$p_value, error = function(e) NA_real_)
  p_D  <- tryCatch(tfu_algo2_sup_test(X_obs)$p_value, error = function(e) NA_real_)
  
  c(p_L2, p_D)
}

# ---- Monte Carlo over b ----
res_list <- vector("list", length(b_vals))
names(res_list) <- sprintf("b=%.1f", b_vals)

for (j in seq_along(b_vals)) {
  b_now <- b_vals[j]
  ps <- replicate(n_sims, one_run(b_now))
  # ps is 2 x n_sims: [1,] = p_L2, [2,] = p_D
  rej_L2 <- mean(ps[1, ] < alpha, na.rm = TRUE)
  rej_D  <- mean(ps[2, ] < alpha, na.rm = TRUE)
  
  res_list[[j]] <- data.frame(
    b   = b_now,
    test = c("T[mu,L^2]", "T[mu,D]"),
    rej = c(rej_L2, rej_D)
  )
  cat(sprintf("b=%.1f  ->  rej L2=%.3f,  rej D=%.3f\n", b_now, rej_L2, rej_D))
}

res_df <- dplyr::bind_rows(res_list)

# ---- Plot (Figure-4-Style, ohne TF) ----
ggplot(res_df, aes(x = b, y = rej, shape = test)) +
  geom_point(size = 3, stroke = 1.1) +
  scale_shape_manual(values = c("T[mu,L^2]" = 1, "T[mu,D]" = 2)) + # 1: o, 2: △
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(
    limits = c(1, 2),
    breaks = seq(1, 2, by = 0.2))+ 
  labs(x = "b", y = "rejection probability",
       shape = NULL,
       title = "Rejection probabilities under MNAR censoring (a = -1, n = 100)") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(.88, .18))
