# ===================== Setup & Packages =====================
# Uses the NEW API (2025-08) for:
#   asym_mean_L2_test(), asym_mean_sup_test()  [asymptotic]
#   boot_mean_test()                            [bootstrap]
# Legacy wrappers asym_conf_bands()/boot_conf_bands() are NOT used.

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
})

# ===================== 1) Tageskurven aus temp_graz ============================
# temp_graz: data.frame mit Spalten time (POSIXct/char), temp_east, temp_west
load("C:/LMU/Bachelor/Bacherlor_FDA/Code/Tests/temp_graz.rda")

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

X_obs <- wide %>% select(all_of(slots48)) %>% as.matrix()

# kurzer Check
obs_per_day <- rowMeans(!is.na(X_obs))
cat("Tage gesamt:", nrow(X_obs),
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

# ===================== 3) Tests & (asymptotische) Bänder =======================
set.seed(2025)

# (A) Asymptotischer L2‑Test (KL‑Mixture)
res_L2 <- asym_mean_L2_test(
  fd = fd, observed_ratio = 1, min_frac = 0.10, seed = 2025, fve = 0.99, B = 5000
)

# (B) Asymptotischer Supremums‑Test; hier OHNE Bänder (leichtgewichtiger)
res_sup <- asym_mean_sup_test(
  fd = fd, observed_ratio = 1, min_frac = 0.10, seed = 2025, fve = 0.99,
  B = 5000, compute_bands = TRUE
)

# # (C) Nur Bänder (leichtgewichtiges Objekt) – ersetzt asym_conf_bands()
# res_cb <- asym_mean_sup_test(
#   fd = fd, observed_ratio = 1, min_frac = 0.10, seed = 2025, alpha = 0.05,
#   compute_bands = TRUE, bands_only = TRUE
# )

# ===================== 4) Ergebnisse ===========================================
n_total <- nrow(X_obs)
# Info: Die Tests gruppieren intern automatisch über observed_ratio=1 (vollständig vs. unvollständig)
# Für die Anzeige behalten wir die einfache Heuristik bei:
group_A <- rowMeans(!is.na(X_obs)) == 1
n_A <- sum(group_A); n_B <- n_total - n_A

m_asym <- as.integer(res_sup$parameter["m"])          # Anzahl Subdomain-Punkte im asym‑Test
m_band <- length(res_sup$idx)                             # dito aus dem Band‑Objekt

cat("\n==== Ergebnisse (observed_ratio = 1 – complete vs incomplete) ====\n")
cat("n =", n_total, "\n")
cat("Größe Gruppe A/B:", n_A, "/", n_B, "\n")
cat("Subdomain-Punkte: asym=", m_asym, "; bands=", m_band,
    " | Fallback:", ifelse(is.null(res_sup$fallback), "none", res_sup$fallback), "\n")

print(res_L2)
print(res_sup)

# ===================== 5) Plots ================================================
# Einheitliche Slots (00:00 .. 23:30) – falls nicht mehr im Scope:
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

# ---------- Rechte Seite: Diff + Bänder MIT LÜCKE -----------------------------
# Sortiere idx und bilde Blöcke zusammenhängender Slots (Lücke trennt Blöcke)
o          <- order(res_sup$idx)
idx_sorted <- res_sup$idx[o]
block      <- c(0, cumsum(diff(idx_sorted) > 1))  # neue Gruppe, wenn Sprung > 1

df_cb <- data.frame(
  time  = factor(slots48[idx_sorted], levels = slots48),  # alle 48 Levels
  diff  = res_sup$diff[o],
  lower = res_sup$lower[o],
  upper = res_sup$upper[o],
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



# ===================== 6) Figure-7-Style (Bootstrap-p-Werte über m) =============
set.seed(2025)
B <- 10000
m_grid <- c(24, 27, 30, 32, 33, 36, 38, 40, 43, 45)

avail <- colSums(!is.na(X_obs))
slot_idx <- seq_len(ncol(X_obs))

res_list <- lapply(m_grid, function(m) {
  # wähle die m Slots mit den meisten Beobachtungen (Paper-Ansatz)
  idx_top <- order(-avail, slot_idx)[1:m] |> sort()
  X_sub   <- X_obs[, idx_top, drop = FALSE]
  
  # Bootstrap-Tests auf diesem Subgrid
  rb_L2 <- boot_mean_test(
    X_obs = X_sub, B = B, min_frac = 0, seed = 2025,  # min_frac=0: keine 10%-Hürde
    parallel = TRUE, stat = "L2", compute_bands = FALSE, return_boot = FALSE
  )
  rb_D  <- boot_mean_test(
    X_obs = X_sub, B = B, min_frac = 0, seed = 2025,
    parallel = TRUE, stat = "D", compute_bands = FALSE, return_boot = FALSE
  )
  
  data.frame(
    m    = m,
    p_L2 = as.numeric(rb_L2$p.value),
    p_D  = as.numeric(rb_D$p.value)
  )
})

df_p <- dplyr::bind_rows(res_list) |>
  tidyr::pivot_longer(cols = c(p_L2, p_D), names_to = "test", values_to = "p")

# Plot im Figure-7-Stil
p_fig7 <- ggplot(df_p, aes(x = m, y = p, shape = test)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_manual(
    values = c(p_L2 = 1, p_D = 2),
    labels = c(p_L2 = expression(T[mu*","*L^2]),
               p_D  = expression(T[mu*","*D]))
  ) +
  scale_x_continuous(limits = c(min(m_grid), max(m_grid)),
                     breaks = m_grid,
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "number of grid points in I", y = "p-values", shape = NULL) +
  theme_classic(base_size = 14)

print(p_fig7)


