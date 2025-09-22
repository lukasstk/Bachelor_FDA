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
  n_boot = 5000, alpha = 0.05,
  parallel = TRUE, stat = "L2",
  compute_bands = FALSE, return_boot = FALSE
)

# (B) Supremums-/D-Test MIT Bootstrap-Bändern (Konfidenzbänder via Bootstrap)
res_sup <- boot_mean_test(
  fd = fd, observed_ratio = 1, min_frac = 0.10,
  n_boot = 5000, alpha = 0.05,
  parallel = TRUE, stat = "D",
  compute_bands = TRUE,    # <-- Bänder per Bootstrap
  return_boot = FALSE
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

p_groups <- ggplot(df_long_plot, aes(x = slot, y = temp, group = date, colour = group)) +
  geom_line(linewidth = 0.7, alpha = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c("Complete" = "grey60", "Incomplete" = "red")) +
  scale_x_discrete(breaks = slots48[tick_idx], labels = slots48[tick_idx], drop = FALSE) +
  labs(title = "Tageskurven (temp_east) – Complete (grau) vs Incomplete (rot)",
       x = "Uhrzeit", y = "Temperatur (°C)", colour = "Gruppe") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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

p_left  <- p_groups
p_right <- ggplot(df_cb, aes(x = time)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_line(aes(y = lower, group = block), linetype = "dashed") +
  geom_line(aes(y = upper, group = block), linetype = "dashed") +
  geom_line(aes(y = diff,  group = block), linewidth = 1) +
  scale_x_discrete(
    breaks = slots48[tick_idx],
    labels = slots48[tick_idx],
    drop   = FALSE
  ) +
  scale_y_continuous(
    limits = c(-10, 10),
    breaks = seq(-10, 10, by = 2)
  ) +
  labs(x = "time", y = "difference in means") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

# Doppel-Plot anzeigen
p_left + p_right

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
    n_boot = B,parallel = TRUE, stat = "L2",
    compute_bands = FALSE)
  
  rb_D <- boot_mean_test(
    X = X_sub, groups = groups_fixed,
    n_boot = B, parallel = TRUE, stat = "D",
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

# 5) Plot im Figure-7-Stil
p_fig7 <- ggplot(df_p, aes(x = m, y = p, shape = test)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_manual(
    values = c(p_L2 = 1, p_D = 2),
    labels = c(p_L2 = expression(T[mu*","*L^2]),
               p_D  = expression(T[mu*","*D]))
  ) +
  scale_x_continuous(limits = c(min(m_grid)-1, max(m_grid)+1),
                     breaks = seq(25,45, by = 5),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "number of grid points in I", y = "p-values", shape = NULL) +
  theme_classic(base_size = 14)

print(p_fig7)
