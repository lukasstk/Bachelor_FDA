# ===================== Setup & Packages =====================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyfun)   # tfd / tf_integrate
  library(patchwork)
  # Optional für Bootstrap:
  library(doParallel); library(doRNG); library(foreach)
  library(extrafont)
  loadfonts(device = "win")
})

# ---- Deine neuen Funktionen müssen im Environment sein ----
# asym_mean_L2_test(), asym_mean_sup_test(), boot_mean_test(), ...

# ===================== 1) Daten laden ======================
load("Data/heart_rate.RData")
# Erwartet: rate_dt mit Spalten id, ti (Stunden 20–26), y0 (bpm)

# ===================== 2) tfd-Objekt & Wide-Matrix ============================
rate_dt <- rate_dt %>% mutate(ti = as.numeric(ti))

fd_hr <- tf::tfd(rate_dt, arg = "ti", id = "id", value = "y0")

# Wide: Zeilen = id, Spalten = Zeitpunkte (direkt mit Pipe)
wide_hr <- rate_dt |>
  tidyr::pivot_wider(
    id_cols   = id,
    names_from  = ti,
    values_from = y0
  )

# Matrix X (ohne id)
X <- wide_hr |>
  dplyr::select(-id) |>
  as.matrix()

row.names(X) <- as.character(wide_hr$id)  

# ===================== 3) Tests (asymptotisch + Bänder) ======================
set.seed(2025)

res_L2 <- asym_mean_L2_test(
  fd = fd_hr, n_sim = 10000
)

res_sup <- asym_mean_sup_test(
  fd = fd_hr, n_sim = 10000
)

res_sup_bands <- asym_mean_sup_test(
  fd = fd_hr, n_sim = 10000, bands_only = TRUE
)

cat("\n=== HEART RATE: Ergebnisse (complete vs incomplete) ===\n")
print(res_L2)
print(res_sup_bands)

# ===================== 4) Plot: Kurvenschar (links) ===========================
# Gruppierung wie in den Tests: complete (alle Zeitpunkte belegt) vs incomplete
group_A <- rowMeans(!is.na(X)) == 1
grp_df  <- tibble(
  id    = row.names(X),
  group = ifelse(group_A, "Complete", "Incomplete")
)

# >>> Fix: gleiche ID-Typen (character) vor dem Join
plot_df <- rate_dt %>%
  select(id, ti, y0) %>%
  mutate(id = as.character(id)) %>%
  inner_join(grp_df, by = "id")

# Gemeinsame Breaks wie im Paper (alle 2 Stunden)
hr_breaks <- c(20, 22, 24, 26)
hr_labels <- c("20:00","22:00","24:00","02:00")

# ===================== 5) Plot: Mean-Diff + Bänder (rechts) ===================
df_band <- data.frame(
  ti    = res_sup$bands$grid,
  diff  = tf::tf_evaluate(res_sup$estimate$mean_diff)[[1]], 
  lower = res_sup$bands$lower,
  upper = res_sup$bands$upper
) %>%
  arrange(ti)


# Lücken in der Subdomäne als getrennte Linienblöcke zeichnen
unique_ti <- sort(unique(df_band$ti))
base_step <- min(diff(unique_ti))
block <- c(0, cumsum(diff(df_band$ti) > 1.5 * base_step))
df_band$block <- block

# ===================== Heart Rate Plots in Publikationsstil =====================

# Links (Kurvenschar)
p_left <- ggplot(plot_df, aes(x = ti, y = y0, group = id)) +
  geom_line(
    data = subset(plot_df, group == "Complete"),
    color = "grey70", linewidth = 0.6, na.rm = TRUE
  ) +
  geom_line(
    data = subset(plot_df, group == "Incomplete"),
    color = "black", linewidth = 0.7, na.rm = TRUE
  ) +
  scale_x_continuous(
    breaks = hr_breaks,
    labels = hr_labels,
    expand = c(0,0)
  ) +
  labs(x = "time", y = "heart rate (bpm)") +
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

# Rechts (Mean-Diff + Bänder)
p_right <- ggplot(df_band, aes(x = ti)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, fill = "grey50") +
  geom_line(aes(y = diff), linewidth = 1, color = "black") +
  geom_line(aes(y = lower), linetype = 2, color = "black") +
  geom_line(aes(y = upper), linetype = 2, color = "black") +
  scale_x_continuous(
    breaks = hr_breaks,
    labels = hr_labels,
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(-40, 40)) +
  labs(x = "time", y = "difference in means") +
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

# Kombinieren mit Spacer
Figure_6_Heart_rate <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# # Abspeichern
# ggsave(
#   filename = "Plots/Figure_6_Heart_rate.png",
#   plot     = Figure_6_Heart_rate,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )
