# ===================== Setup & Packages =====================
# Uses the NEW API (2025-08) for:
#   boot_mean_test()                            [bootstrap]
#   asym_mean_L2_test(), asym_mean_sup_test()  [optional/asymptotic, not used for bands]
# Konfidenzbänder werden unten explizit per Bootstrap erzeugt.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(tidyfun)   # tfd / tf_integrate
  library(patchwork)
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(pbapply)   # pblapply()
  # library(tf)      # optional; wir nutzen tf::tfd() qualifiziert
  library(extrafont)
  loadfonts(device = "win")
})

# ===================== 1) Tageskurven aus temp_graz ============================
# temp_graz: data.frame mit Spalten time (POSIXct/char), temp_east, temp_west
load("Data/temp_graz.rda")

# Ensure POSIXct and build date + slot (HH:MM)
df <- temp_graz %>%
  mutate(
    time = as.POSIXct(time, tz = "UTC"),
    date = as.Date(time),
    slot = sprintf("%02d:%02d", hour(time), minute(time))
  ) %>%
  arrange(time)

# 48 Halbstunden-Slots pro Tag: 00:00..23:30
slots48 <- sprintf("%02d:%02d", rep(0:23, each = 2), rep(c(0, 30), times = 24))

# Wide-Matrix: Zeilen = Tage, Spalten = Halbstunden-Slots, Werte = temp_east
wide <- df %>%
  filter(slot %in% slots48) %>%
  mutate(slot = factor(slot, levels = slots48)) %>%
  select(date, slot, temp = temp_east) %>%
  distinct() %>%
  pivot_wider(names_from = slot, values_from = temp) %>%
  arrange(date)

X <- wide %>% select(all_of(slots48)) %>% as.matrix()

# kurzer Check
obs_per_day <- rowMeans(!is.na(X))
cat("Tage gesamt:", nrow(X),
    " | Median beobachtet/Tag (Anteil):", signif(median(obs_per_day), 3), "\n")

# ===================== 2) tfd-Objekt (irregulär) ===============================
slot_to_hours <- function(s) {
  hh <- as.numeric(substr(s, 1, 2))
  mm <- as.numeric(substr(s, 4, 5))
  hh + mm / 60
}

df_long <- wide %>%
  pivot_longer(cols = all_of(slots48), names_to = "slot", values_to = "temp") %>%
  mutate(hour_num = slot_to_hours(slot))

fd <- tf::tfd(df_long, arg = "hour_num", id = "date", value = "temp")

# ===================== 3) Tests & (BOOTSTRAP-)Bänder ===========================
set.seed(2025)

# (A) L2-Test (optional, ohne Bänder)
res_L2 <- boot_mean_test(
  fd = fd, observed_ratio = 1, min_frac = 0.10,
  n_boot = 10000, alpha = 0.05, stat = "L2",
  compute_bands = FALSE
)

# (B) Supremums-/D-Test MIT Bootstrap-Bändern (Konfidenzbänder via Bootstrap)
res_sup <- boot_mean_test(
  fd = fd, observed_ratio = 1, min_frac = 0.10,
  n_boot = 10000, alpha = 0.05, stat = "D",
  compute_bands = TRUE
)

# ===================== 4) Ergebnisse ===========================================
n_total <- nrow(X)

# Einfache Heuristik zur Visualisierung: complete vs incomplete
group_A <- rowMeans(!is.na(X)) == 1
n_A <- sum(group_A); n_B <- n_total - n_A

m_band <- length(res_sup$idx)  # Subdomain-Punkte aus den Bootstrap-Bändern

cat("\n==== Ergebnisse (observed_ratio = 1 – complete vs incomplete) ====\n")
cat("n =", n_total, "\n")
cat("Größe Gruppe A/B:", n_A, "/", n_B, "\n")
cat("Subdomain-Punkte (Bands) =", m_band, "\n")

print(res_L2)
print(res_sup)

# ===================== 5) Plots ================================================
# Einheitliche Slots (00:00 .. 23:30)
slots48 <- sprintf("%02d:%02d", rep(0:23, each = 2), rep(c(0,30), 24))
tick_idx <- seq(1, length(slots48), by = 6)  # Ticks alle 3h

wide$group <- ifelse(group_A, "Complete", "Incomplete")

df_long_plot <- wide %>%
  tidyr::pivot_longer(cols = all_of(slots48), names_to = "slot", values_to = "temp") %>%
  mutate(slot = factor(slot, levels = slots48))

# ---------- Rechte Seite: Diff + Bootstrap-Bänder MIT LÜCKE --------------------
grid_sub  <- res_sup$bands$grid
grid_full <- slot_to_hours(slots48)

# Finde Positionen des Subgrids im vollen Grid
idx_out <- vapply(grid_sub, function(g) which.min(abs(grid_full - g)), integer(1))

# Sortieren + Blöcke bilden
o          <- order(idx_out)
idx_sorted <- idx_out[o]
block      <- c(0, cumsum(diff(idx_sorted) > 1))

# Extrahiere Differenz und Bänder
diff_vec  <- as.numeric(as.matrix(res_sup$estimate$mean_diff))[o]
lower_vec <- as.numeric(res_sup$bands$lower[o])
upper_vec <- as.numeric(res_sup$bands$upper[o])

# Dataframe bauen
df_cb <- data.frame(
  time  = factor(slots48[idx_sorted], levels = slots48),
  diff  = diff_vec,
  lower = lower_vec,
  upper = upper_vec,
  block = block
)

# ===================== Final Plot: Figure_6_temp_data ===========================

# Links (Tageskurven, grau vs schwarz)
p_left <- ggplot(df_long_plot, aes(x = slot, y = temp, group = date)) +
  geom_line(
    data = subset(df_long_plot, group == "Complete"),
    color = "grey70", linewidth = 0.6, alpha = 0.9, na.rm = TRUE
  ) +
  geom_line(
    data = subset(df_long_plot, group == "Incomplete"),
    color = "black", linewidth = 0.7, alpha = 0.9, na.rm = TRUE
  ) +
  scale_x_discrete(
    breaks = c("00:00", "06:00", "12:00", "18:00", "23:30"),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
    drop   = FALSE,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = function(x) setdiff(pretty(x), 10)  
  ) + 
  labs(x = "time", y = expression("temperature ("*degree*C*")")) +
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

# Rechts (Differenzen + Bootstrap-Bänder, schwarz/grau)
p_right <- ggplot(df_cb, aes(x = time)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper, group = block),
              alpha = 0.15, fill = "grey50") +
  geom_line(aes(y = diff, group = block), linewidth = 1, color = "black") +
  geom_line(aes(y = lower, group = block), linetype = 2, color = "black") +
  geom_line(aes(y = upper, group = block), linetype = 2, color = "black") +
  scale_x_discrete(
    breaks = c("00:00", "06:00", "12:00", "18:00", "23:30"),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
    drop   = FALSE,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(x = "time", y = expression("Difference in means ("*degree*C*")")) +
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

# Kombination
Figure_6_temp_data <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# Abspeichern
# ggsave(
#   filename = "Plots/Figure_6_temp_data.png",
#   plot     = Figure_6_temp_data,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )


# ==== Figure-7-Style (Bootstrap-p-Werte über m) – FIX: Gruppen manuell/fix ====

set.seed(2025)
B <- 10000
m_grid <- c(24,26,28,31,33,36,38,40,43,45)

# Wichtig: Gruppen (complete vs incomplete) werden EINMAL auf dem vollen Gitter
# bestimmt und eingefroren. Andernfalls könnten Kurven durch Subgrids ihre
# Gruppenzugehörigkeit wechseln

# 1) Gruppen EINFRIEREN auf dem vollen Gitter (48 Slots):
groups_fixed <- rowMeans(!is.na(X)) == 1    

# 2) Slots nach Gesamtverfügbarkeit sortieren (Paper-Ansatz)
avail    <- colSums(!is.na(X))
slot_idx <- seq_len(ncol(X))

# 3) p-Werte über m berechnen – mit fixen Gruppen + min_frac = 0
res_list <- pblapply(m_grid, function(m) {
  # Top-m Slots nach Gesamt-Verfügbarkeit
  idx_top <- order(-avail, slot_idx)[1:m] |> sort()
  X_sub   <- X[, idx_top, drop = FALSE]
  
  # Bootstrap-Tests auf diesem Subgrid
  #  - groups = groups_fixed (fixe Gruppenzugehörigkeit)
  #  - min_frac = 0 (paper-konform, keine zusätzliche 10%-Hürde), opt. hab die Hürde drin gelassen
  rb_L2 <- boot_mean_test(
    X = X_sub, groups = groups_fixed,
    n_boot = B, stat = "L2",
    compute_bands = FALSE)
  
  rb_D <- boot_mean_test(
    X = X_sub, groups = groups_fixed,
    n_boot = B, stat = "D",
    compute_bands = FALSE)
  
  data.frame(
    m    = m,
    p_L2 = as.numeric(rb_L2$p.value),
    p_D  = as.numeric(rb_D$p.value)
  )
})

# 4) Aufbereiten & Tabelle
df_p <- dplyr::bind_rows(res_list) |>
  tidyr::pivot_longer(cols = c(p_L2, p_D), names_to = "test", values_to = "p")

df_table <- df_p |>
  tidyr::pivot_wider(names_from = test, values_from = p) |>
  dplyr::rename(`T[mu,L2]` = p_L2, `T[mu,D]` = p_D)

print(df_table, n = Inf)

# Figure 7 im Figure-4-Stil
Figure_7 <- ggplot(df_p, aes(x = m, y = p, shape = test, color = test)) +
  geom_point(size = 5, stroke = 1.5, shape = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_color_manual(
    values = c("p_L2" = "steelblue",
               "p_D"  = "red"),
    breaks = c("p_L2", "p_D"),
    labels = c(
      "p_L2" = expression(T[mu]*","~L^2),
      "p_D"  = expression(T[mu]*","~D)
    )
  ) +
  scale_x_continuous(
    limits = c(min(m_grid)-1, max(m_grid)+1),
    breaks = seq(25, 45, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "number of grid points in I", y = "p-values",
       shape = NULL, color = NULL) +
  theme_classic(base_size = 18, base_family = "Times New Roman") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.91, 0.3),
    legend.margin = margin(7, 7, 7, 7),
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    axis.title.x = element_text(size = 25, margin = margin(t = 30)),
    axis.title.y = element_text(size = 25, margin = margin(r = 30)),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key.size = unit(1, "lines"),             
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.major = element_line(colour = "grey80", size = 0.5),
    panel.grid.minor = element_blank()
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4))
  )

print(Figure_7)

# Speichern
# ggsave(
#   filename = "Plots/Figure_7.png",
#   plot     = Figure_7,
#   width    = 10,
#   height   = 5,
#   dpi      = 300
# )

