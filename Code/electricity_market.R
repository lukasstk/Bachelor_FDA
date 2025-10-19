# ===================== Setup & Packages =====================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tf)          # für tf::tfd / tf::tf_integrate
  library(tidyfun)
  library(patchwork)
  library(extrafont)
  loadfonts(device = "win")
})

# ---- Funktionen: asym_mean_L2_test(), asym_mean_sup_test() ----
# müssen im Environment vorhanden sein

# ===================== 1) Daten laden ======================
data <- as.matrix(read.csv("Data/logbidcurves.csv", header = TRUE, check.names = FALSE))

grid_full <- seq(1700, 2500, length.out = 161)

# Orientierung: Zeilen = Kurven, Spalten = Gridpunkte
X_full <- if (ncol(data) == length(grid_full)) data else t(data)

# ===================== 2) Tests (nur deine Funktionen) ======================
set.seed(2025)
fd_full <- tf::tfd(X_full, arg = grid_full)
observed_ratio <- 0.43

res_L2 <- asym_mean_L2_test(
  fd = fd_full,
  observed_ratio = observed_ratio,  
  min_frac = 0.10,
  fve = 0.99,
  n_sim = 10000,
  seed = 42
)

res_sup <- asym_mean_sup_test(
  fd = fd_full,
  observed_ratio = observed_ratio,   
  min_frac = 0.10,
  fve = 0.99,
  n_sim = 10000,
  compute_bands = TRUE,
  seed = 42
)

cat(sprintf("\n=== Ergebnisse auf I = [%.0f, %.0f] MW ===\n",
            min(res_sup$grid), max(res_sup$grid)))
print(res_L2)
print(res_sup)

# ===== Labels passend zu δ =====
lab_A <- sprintf("obs \u2265 %.0f%%", observed_ratio * 100)  # "obs ≥ 43%"
lab_B <- sprintf("obs < %.0f%%",  observed_ratio * 100)      # "obs < 43%"

# Falls noch nicht gesetzt:
group_A    <- rowMeans(!is.na(X_full)) >= observed_ratio
group_labs <- ifelse(group_A, lab_A, lab_B)

# ===================== 4) Plots =============================================
# (a) Linkes Panel: Kurven farbig nach Gruppen
df_long <- as.data.frame(X_full)
colnames(df_long) <- as.character(grid_full)
df_long$id    <- sprintf("curve_%05d", seq_len(nrow(df_long)))
df_long$group <- factor(group_labs, levels = c(lab_A, lab_B))  # feste Reihenfolge

df_long <- tidyr::pivot_longer(
  df_long, cols = -c(id, group),
  names_to = "demand", values_to = "logprice"
) |> dplyr::mutate(demand = as.numeric(demand))

df_band <- tidyfun::tf_unnest(res_sup$estimate$mean_diff) %>%
  dplyr::select(demand = arg, diff = value) %>%
  dplyr::mutate(
    lower = res_sup$bands$lower,
    upper = res_sup$bands$upper
  ) %>%
  tidyr::drop_na() %>%
  dplyr::arrange(demand)

# Links
p_left <- ggplot(df_long, aes(x = demand, y = logprice, group = id)) +
  geom_line(
    data = subset(df_long, group == lab_A),
    color = "grey70", linewidth = 0.6, na.rm = TRUE
  ) +
  geom_line(
    data = subset(df_long, group == lab_B),
    color = "black", linewidth = 0.7, na.rm = TRUE
  ) +
  scale_x_continuous(
    limits  = c(1700, 2500),
    breaks  = c(1800, 2000, 2200, 2400),
    expand  = c(0,0)
  ) +
  labs(x = "demand (MW)", y = "log price") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    legend.position   = "none",
    text              = element_text(family = "Times New Roman"),
    axis.title        = element_text(size = 20),
    axis.text         = element_text(size = 16),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# Rechts
p_right <- ggplot(df_band, aes(x = demand)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, fill = "grey50") +
  geom_line(aes(y = diff), linewidth = 1, color = "black") +
  geom_line(aes(y = lower), linetype = 2, color = "black") +
  geom_line(aes(y = upper), linetype = 2, color = "black") +
  scale_x_continuous(
    breaks = c(1800, 2000, 2200, 2400),
    expand = c(0,0)
  ) +
  coord_cartesian(xlim = c(1700, 2500), ylim = c(-4, 4)) +
  labs(x = "demand (MW)", y = "difference in means") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    text              = element_text(family = "Times New Roman"),
    axis.title        = element_text(size = 18),
    axis.text         = element_text(size = 16),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# Kombi
electricity_market_plot <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# # Abspeichern
# ggsave(
#   filename = "Plots/electricity_market_plot.png",
#   plot     = electricity_market_plot,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )
