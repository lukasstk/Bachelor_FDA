# ===============================================================
# Brownian-Motion Testdaten + Anwendung der Algorithmen (1–3,5)
# -> angepasst: Input als <tfd>
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
mechanism    <- "MCAR"                 # "MCAR" oder "MNAR"

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
bm_mat <- t(replicate(n, simulate_bm(grid)))
rownames(bm_mat) <- ids

# Beobachtungsmatrix nach Mechanismus
if (mechanism == "MCAR") {
  O_mat <- t(sapply(seq_len(n), function(i) make_O_mcar(grid)))
} else if (mechanism == "MNAR") {
  O_mat <- t(apply(bm_mat, 1L, make_O_mnar))
} else stop("Unknown mechanism")

storage.mode(O_mat) <- "integer"
rownames(O_mat) <- ids

# Gruppen nach Zensur
is_complete <- rowSums(O_mat) == length(grid)
A_ids <- ids[is_complete]
B_ids <- ids[!is_complete]

# Beobachtete Werte (NA außerhalb Beobachtung)
X_obs <- bm_mat
X_obs[O_mat == 0L] <- NA_real_

# ---------------- <tfd> Objekt ----------------
fd <- tf::tfd(X_obs, arg = grid, id = ids)

# Gruppenvektor für Tests
groups_vec <- ifelse(ids %in% A_ids, "A", "B")

# ---------------- Algorithmen aufrufen ----------------
B_mc  <- 2000
alpha <- 0.05

# Algo 1: L2-Test
res_L2 <- tfu_algo1_L2_test(fd = fd, groups = groups_vec, B = B_mc)
print(res_L2)

# Algo 2: Supremum-Test
res_D <- tfu_algo2_sup_test(fd = fd, groups = groups_vec, B = B_mc)
print(res_D)

# Algo 3: Simultane Konfidenzbänder
res_CB <- tfu_algo3_conf_bands(fd = fd, groups = groups_vec, alpha = alpha, B = B_mc)
bands <- data.frame(
  t     = res_CB$grid,
  diff  = res_CB$diff,
  lower = res_CB$lower,
  upper = res_CB$upper
)

# Algo 5: Bootstrap (Matrix-Input)
res_boot <- tfu_algo5_bootstrap(X_obs = X_obs)
cat("Algo5 Boot -> L2: stat", res_boot$T_L2, " p", res_boot$p_L2,
    " | D: stat", res_boot$T_D, " p", res_boot$p_D, "\n")

# ---------------- Plot: links Kurven, rechts Konfidenzbänder ----------------
df_obs <- as.data.frame(X_obs) %>%
  mutate(id = rownames(X_obs),
         gruppe = ifelse(rowSums(is.na(.)) == 0, "vollständig", "unvollständig")) %>%
  pivot_longer(cols = -c(id, gruppe), names_to = "t_idx", values_to = "x") %>%
  mutate(t = grid[as.integer(sub("V", "", t_idx))])

p_left <- ggplot(df_obs, aes(x = t, y = x, group = id)) +
  geom_line(aes(color = gruppe), linewidth = 0.9, na.rm = TRUE) +
  scale_color_manual(values = c("vollständig" = "grey60", "unvollständig" = "black")) +
  labs(title = "Brownian-Pfade nach Zensur", x = "Zeit t", y = "X(t)", color = "Gruppe") +
  theme_minimal() + theme(legend.position = "bottom")

p_right <- ggplot(bands, aes(t, diff)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = lower), linetype = 2) +
  geom_line(aes(y = upper), linetype = 2) +
  labs(x = "Zeit", y = "Differenz der Mittelwerte") +
  theme_minimal()

p_left + p_right
