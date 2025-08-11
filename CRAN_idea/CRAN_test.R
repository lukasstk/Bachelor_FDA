# ===============================================================
# Brownian-Motion Testdaten + Anwendung der Algorithmen (1–3,5)
# ===============================================================

library(dplyr)
library(tidyr)
library(tidyfun)
library(tf)

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
    rep(1L, length(grid))                       # komplett beobachtet
  } else {
    U <- sort(runif(2)); as.integer(grid >= U[1] & grid < U[2])  # Intervall
  }
}

make_O_mnar <- function(x_row) {
  # beobachtet nur, wenn X in (-1, 2); voll wenn gesamte Kurve dort liegt
  if (all(x_row > -1 & x_row < 2)) rep(1L, length(x_row))
  else as.integer(x_row > -1 & x_row < 2)
}

# ---------------- Daten simulieren ----------------
# Brownian Paths (unter H0 ohne Mean-Shift)
bm_mat <- t(replicate(n, simulate_bm(grid)))
rownames(bm_mat) <- ids

# Gruppen: A = komplett, B = unvollständig (päckchenweise Hälfte/Hälfte)
A_ids <- sample(ids, size = n/2)
B_ids <- setdiff(ids, A_ids)

# Beobachtungsmatrix je nach Mechanismus erzeugen
O_mat <- matrix(0L, nrow = n, ncol = length(grid), dimnames = list(ids, NULL))
for (i in seq_len(n)) {
  if (ids[i] %in% A_ids) {
    # Gruppe A: komplett beobachtet
    O_mat[i, ] <- 1L
  } else {
    # Gruppe B: je nach Mechanismus
    if (mechanism == "MCAR") {
      O_mat[i, ] <- make_O_mcar(grid)
    } else if (mechanism == "MNAR") {
      O_mat[i, ] <- make_O_mnar(bm_mat[i, ])
    } else stop("Unknown mechanism")
  }
}

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

# ---------------- Deine Pipeline (unverändert) ----------------
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

# (NAs in X_obs bleiben wie sie sind — so erwartet von den Funktionen)

# ---------------- Algorithmen aufrufen ----------------
# Hinweis: B in Tests ruhig kleiner (z. B. 1000); final 5000+
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
head(bands)

# Algo 5: Bootstrap-P-Werte (L2 & Supremum)
res_boot <- tfu_algo5_bootstrap(
  X_obs = X_obs
)
cat("Algo5 Boot -> L2: stat", res_boot$T_L2, " p", res_boot$p_L2,
    " | D: stat", res_boot$T_D, " p", res_boot$p_D, "\n")
