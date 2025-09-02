# ===================== Setup & Packages =====================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tf)          # für tf::tfd / tf::tf_integrate
  library(tidyfun)
  library(patchwork)
})

# ---- Funktionen: asym_mean_L2_test(), asym_mean_sup_test() ----
# müssen im Environment vorhanden sein

# ===================== 1) Daten laden ======================
csv_path <- "Data/logbidcurves.csv"
M <- as.matrix(read.csv(csv_path, header = TRUE, check.names = FALSE))
storage.mode(M) <- "numeric"

grid_full <- seq(1700, 2500, length.out = 161)

# Orientierung: Zeilen = Kurven, Spalten = Gridpunkte
X_full <- if (ncol(M) == length(grid_full)) M else t(M)
rm(M)

# ===================== 2) Tests (nur deine Funktionen) ======================
set.seed(2025)
fd_full <- tf::tfd(X_full, arg = grid_full)

res_L2 <- asym_mean_L2_test(
  fd = fd_full,
  observed_ratio = 0.43,   # Gruppenbildung intern
  min_frac = 0.10,
  fve = 0.99,
  B = 5000,
  seed = 2025
)

res_sup <- asym_mean_sup_test(
  fd = fd_full,
  observed_ratio = 0.43,   # Gruppenbildung intern
  min_frac = 0.10,
  fve = 0.99,
  B = 5000,
  compute_bands = TRUE,
  seed = 2025
)

cat(sprintf("\n=== Ergebnisse auf I = [%.0f, %.0f] MW ===\n",
            min(res_sup$grid), max(res_sup$grid)))
print(res_L2)
print(res_sup)

# ===================== 3) Gruppenlabels nur für Plot ========================
observed_ratio <- 0.43
group_A <- rowMeans(!is.na(X_full)) >= observed_ratio
group_labs <- ifelse(group_A, "Large extent", "Small extent")

# ===================== 4) Plots =============================================
# (a) Linkes Panel: Kurven farbig nach Gruppen
df_long <- as.data.frame(X_full)
colnames(df_long) <- as.character(grid_full)
df_long$id <- sprintf("curve_%05d", seq_len(nrow(df_long)))
df_long$group <- group_labs

df_long <- tidyr::pivot_longer(df_long, cols = -c(id, group),
                               names_to = "demand", values_to = "logprice") |>
  mutate(demand = as.numeric(demand))

p_left <- ggplot(df_long, aes(x = demand, y = logprice, group = id, colour = group)) +
  geom_line(linewidth = 0.4, alpha = 0.85, na.rm = TRUE) +
  scale_colour_manual(values = c("Large extent" = "grey70", "Small extent" = "black")) +
  scale_x_continuous(limits = c(1700, 2500),
                     breaks = c(1800, 2000, 2200, 2400),
                     expand = c(0,0)) +
  labs(x = "demand (MW)", y = "log price") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# (b) Rechts: Mean-Diff + simultane Bänder
df_band <- data.frame(
  demand = res_sup$grid,
  diff   = as.numeric(res_sup$diff),
  lower  = as.numeric(res_sup$lower),
  upper  = as.numeric(res_sup$upper)
) |> tidyr::drop_na() |> arrange(demand)

p_right <- ggplot(df_band, aes(x = demand)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_line(aes(y = lower), linetype = "dashed", na.rm = TRUE) +
  geom_line(aes(y = upper), linetype = "dashed", na.rm = TRUE) +
  geom_line(aes(y = diff),  linewidth = 0.9, na.rm = TRUE) +
  scale_x_continuous(breaks = c(1800, 2000, 2200, 2400), expand = c(0,0)) +
  coord_cartesian(xlim = c(1700, 2500), ylim = c(-4, 4)) +
  labs(x = "demand (MW)", y = "difference in means") +
  theme_classic(base_size = 14)

# Anzeigen
p_left + p_right
