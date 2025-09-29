# ===============================================================
# Brownian-Motion Testdaten + Anwendung der neuen Funktionen
# (asym_mean_L2_test, asym_mean_sup_test, boot_mean_test)
# → Keine eigenen Gruppen: Auto-Gruppierung via observed_ratio = 1 (Default)
# ===============================================================

library(dplyr)
library(tidyr)
library(tidyfun)
library(tf)
library(ggplot2)
library(patchwork)
library(foreach)
library(doRNG)
library(pbapply)
library(checkmate)

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
  O <- t(sapply(seq_len(n), function(i) make_O_mcar(grid)))
} else if (mechanism == "MNAR") {
  O <- t(apply(bm_mat, 1L, make_O_mnar))
} else {
  stop("Unknown mechanism")
}
storage.mode(O) <- "integer"
rownames(O) <- ids
colnames(O) <- NULL

# Beobachtete Werte (NA außerhalb Beobachtung)
X <- bm_mat
X[O == 0L] <- NA_real_

# ---------------- Long-Format bauen (df_long) ----------------
# Spaltennamen als echte t-Werte, damit pivot_longer sauber wird
colnames(X) <- formatC(grid, digits = 12, format = "fg", flag = "#")

df_long <- as.data.frame(X) |>
  mutate(id = ids) |>
  tidyr::pivot_longer(
    cols = -id,
    names_to = "t",
    values_to = "x"
  ) |>
  mutate(t = as.numeric(t)) |>
  arrange(id, t)

# ---------------- Pipeline (bis zur Matrix) ----------------
# 1) Grid & IDs
grid <- sort(unique(df_long$t))
ids  <- unique(df_long$id)

# 2) Wide-Matrix X (n x m)
df_wide <- df_long %>%
  select(id, t, x) %>%
  mutate(t = as.numeric(t)) %>%
  tidyr::pivot_wider(names_from = t, values_from = x) %>%
  arrange(match(id, ids))

ord   <- order(as.numeric(names(df_wide)[-1]))
X <- as.matrix(df_wide[, c(1 + ord)])
rownames(X) <- df_wide$id
colnames(X) <- NULL

# 3) Beobachtungsmatrix O aus NA-Maske
O <- 1L * !is.na(X)

# ---------------- Neue Algorithmen aufrufen (Auto-Gruppierung) ----------------
# Algo 1 (neu): L2-Test (asymptotisch)
res_L2 <- asym_mean_L2_test(
  X   = X,
  # keine groups → Auto-Gruppierung (observed_ratio = 1)
  fve     = 0.99,
  n_sim   = 10000,       # MC für KL-Mischung (p-Wert)
  min_frac= 0.10,
  seed    = 42
)
cat("Algo1 L2  -> stat:", unname(res_L2$statistic), " p:", res_L2$p.value, "\n")

# Algo 2 (neu): Supremums-Test (asymptotisch) + simultane Bänder
res_sup <- asym_mean_sup_test(
  X        = X,
  # keine groups → Auto-Gruppierung (observed_ratio = 1)
  fve          = 0.99,
  n_sim        = 10000,
  min_frac     = 0.10,
  seed         = 42,
  alpha        = alpha,
  compute_bands= TRUE,
  bands_only   = FALSE
)
cat("Algo2 Sup -> stat:", unname(res_sup$statistic), " p:", res_sup$p.value, "\n")

# Bänder (kompakte Liste)
bands_list <- asym_mean_sup_test(
  X        = X,
  # keine groups → Auto-Gruppierung (observed_ratio = 1)
  fve          = 0.99,
  n_sim        = 10000,
  min_frac     = 0.10,
  seed         = 42,
  alpha        = alpha,
  compute_bands= TRUE,
  bands_only   = TRUE
)

bands <- data.frame(
  t     = bands_list$bands$grid,
  diff  = tf::tf_evaluate(bands_list$estimate$mean_diff, arg = bands_list$bands$grid)[[1]],
  lower = bands_list$bands$lower,
  upper = bands_list$bands$upper
)

# Algo 5 (neu): Bootstrap-P-Werte – getrennte Aufrufe für L2 und D
res_boot_L2 <- boot_mean_test(
  X        = X,
  # keine groups → Auto-Gruppierung (observed_ratio = 1)
  n_boot       = 10000,
  min_frac     = 0.10,
  alpha        = alpha,
  stat         = "L2",
  compute_bands= TRUE,
  manage_backend = "auto",
  seed         = 123
)

res_boot_D <- boot_mean_test(
  X        = X,
  # keine groups → Auto-Gruppierung (observed_ratio = 1)
  n_boot        = 10000,
  min_frac     = 0.10,
  alpha        = alpha,
  stat         = "D",
  compute_bands= TRUE,   # hier optional zugleich Bandbreiten-Quantil
  manage_backend = "auto",
  seed         = 123
)

cat(
  "Algo5 Boot -> L2: stat", unname(res_boot_L2$statistic), " p", res_boot_L2$p.value,
  " | D: stat", unname(res_boot_D$statistic), " p", res_boot_D$p.value, "\n"
)

# ---------------- Plot: links Kurven (vollständig vs. unvollständig nur für Anzeige),
# rechts Diff + 95%-Band ----------------
library(ggplot2)
library(patchwork)

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
  labs(x = "time", y = "X(t)", color = "Group") +
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

Figure_5 <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# Abspeichern
ggsave(
  filename = "Plots/Figure_5.png",
  plot     = Figure_5,
  width    = 10,
  height   = 4,
  dpi      = 300
)



# ============================================================================
# Plot: Rejection-Probs vs. b (MNAR, a = -1 < b) – nur T_{mu,L2} und T_{mu,D}
# (asymptotische Varianten der neuen Funktionen; ohne eigene Gruppen)
# ============================================================================

library(dplyr)
library(ggplot2)


# ---- Simulation settings ----
n        <- 100                 # sample size per run (wie im Paper)
m_grid   <- 100                 # grid points
grid     <- seq(0, 1, length.out = m_grid)
alpha    <- 0.05
b_vals   <- seq(1.0, 2.0, by = 0.1)   # vary b
n_sims   <- 1000                    # für Speed; Paper ~5000

# ---- Helpers ----
simulate_bm <- function(grid) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc))
}

# MNAR censoring: observe only when a < X(t) < b  (a = -1)
make_missing_pattern <- function(x_row, a = -1, b = 2) {
  as.integer(x_row > a & x_row < b)
}

one_run <- function(b_now) {
  # 1) simulate n Brownian paths on grid
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  # 2) build observation mask under MNAR censoring
  O  <- t(apply(bm_mat, 1L, make_missing_pattern, a = -1, b = b_now))
  storage.mode(O) <- "integer"
  
  # 3) observed matrix with NAs
  X <- bm_mat
  X[O == 0L] <- NA_real_
  
  # 4) keine expliziten Gruppen – Auto-Gruppierung (observed_ratio=1)
  p_L2 <- tryCatch(
    asym_mean_L2_test(X = X, n_sim = 10000)$p.value,
    error = function(e) NA_real_
  )
  p_D  <- tryCatch(
    asym_mean_sup_test(X = X, n_sim = 10000, compute_bands = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  c(p_L2, p_D)
}

# ---- Monte Carlo over b ----
res_list <- pblapply(b_vals, function(b_now) {
  # 2 x n_sims: [1,] = p_L2, [2,] = p_D
  ps  <- pbreplicate(n_sims, one_run(b_now))
  rej <- rowMeans(ps < alpha, na.rm = TRUE)  # c(L2, D)
  
  cat(sprintf("b=%.1f  ->  rej L2=%.3f,  rej D=%.3f\n",
              b_now, rej[1], rej[2]))
  
  data.frame(
    b    = b_now,
    test = c("T[mu,L^2]", "T[mu,D]"),
    rej  = rej,
    row.names = NULL
  )
})

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

library(ggplot2)
library(extrafont)

#Figure 4: Rejection Probabilities
rejection_proba <- readRDS("C:/LMU/Bachelor/Bacherlor_FDA/results/Figure_4_5/res_df_20250913-193109.rds")

Figure_4 <- ggplot(rejection_proba, aes(x = b, y = rej, shape = test, color = test)) +
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
    axis.title.x = element_text(size = 25, margin = margin(t = 30)),  # größer + mehr Abstand
    axis.title.y = element_text(size = 25, margin = margin(r = 30)),  # größer + mehr Abstand
    axis.text  = element_text(size = 18),                             # Ticklabels größer
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




ggsave(
  filename = "Plots/Figure_4.png",
  plot     = Figure_4,
  width    = 10,    
  height   = 5,
  dpi      = 300
)


