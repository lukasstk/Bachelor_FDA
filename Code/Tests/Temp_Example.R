# Pakete
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyfun)   # deine Algorithmen nutzen tf::tfd / tf_integrate

# relativer Pfad von deinem aktuellen Working Directory
file <- "Code/Tests/luis-daten.csv"

# Missing nach Paper: Sensors’ sensitivity to high temperature values or 
# caused by a power problem in the transmission unit
# Kann ich nicht so bestätigen siehe Plot

# Headerzeile finden
hdr_line <- which(grepl("^Datum;Uhrzeit;Wert\\s*$", readLines(file, warn = FALSE)))[1]

dat <- read_delim(
  file, delim = ";", skip = hdr_line - 1,
  locale = locale(decimal_mark = ","),
  col_types = cols(
    Datum   = col_character(),
    Uhrzeit = col_character(),
    Wert    = col_number()
  )
)

# --- 2) Tageskurven bauen: 24 Halbstunden-Slots (00:30..23:30) ------------------
slots <- sprintf("%02d:30", 0:23)

wide <- dat |>
  filter(Uhrzeit %in% slots) |>
  mutate(slot = factor(Uhrzeit, levels = slots)) |>
  select(Datum, slot, Wert) |>
  distinct() |>
  pivot_wider(names_from = slot, values_from = Wert) |>
  arrange(Datum)

X_obs <- wide |> select(all_of(slots)) |> as.matrix()

# kurzer Check
obs_per_day <- rowSums(!is.na(X_obs))
cat("Tage gesamt:", nrow(X_obs),
    " | Median beobachtet/Tag:", median(obs_per_day), "\n")

# --- 3) Analyse mit δ = 1 (complete vs incomplete) ----------------
set.seed(2025)
res_L2  <- tfu_algo1_L2_test(X_obs, delta_A = 1, min_frac = 0.10, seed = 2025)
res_sup <- tfu_algo2_sup_test(X_obs, delta_A = 1, min_frac = 0.10, seed = 2025)
res_cb  <- tfu_algo3_conf_bands(X_obs, delta_A = 1, min_frac = 0.10, seed = 2025)
res_boot <- tfu_algo5_bootstrap(
  X_obs,
  delta_A   = 1,       # complete vs incomplete
  min_frac  = 0.10,    # 10%-Subdomain
  B         = 2000,    # Anzahl Bootstrap-Wiederholungen
  seed      = 2025,
  return_boot = TRUE   # falls du die Bootstraps behalten willst
)

# --- 4) Ergebnisse ausgeben -----------------------------------------------------
cat("\n==== Ergebnisse (δ = 1 – complete vs incomplete) ====\n")
cat("n =", nrow(X_obs), "\n")
cat("Größe Gruppe A/B:",
    sum(rowMeans(!is.na(X_obs)) >= 1), "/",
    sum(rowMeans(!is.na(X_obs)) < 1), "\n")
cat("Subdomain-Punkte:", length(res_L2$idx),
    " | Fallback:", ifelse(is.null(res_L2$fallback), "none", res_L2$fallback), "\n")
cat("L2-Test  p-value:", signif(res_L2$p_value, 4),
    " | Stat:", signif(res_L2$stat, 4), "\n")
cat("Sup-Test p-value:", signif(res_sup$p_value, 4),
    " | Stat:", signif(res_sup$stat, 4), "\n")
# Ergebnisse ausgeben
cat("\n==== Bootstrap-Ergebnisse (δ = 1) ====\n")
cat("n =", nrow(X_obs), "\n")
cat("Größe Gruppe A/B:",
    sum(rowMeans(!is.na(X_obs)) >= 1), "/",
    sum(rowMeans(!is.na(X_obs)) < 1), "\n")
cat("Subdomain-Punkte:", length(res_boot$idx),
    " | Fallback:", ifelse(is.null(res_boot$fallback), "none", res_boot$fallback), "\n")
cat("L2-Test  p-value (boot):", signif(res_boot$p_L2, 4),
    " | Stat:", signif(res_boot$T_L2, 4), "\n")
cat("Sup-Test p-value (boot):", signif(res_boot$p_D, 4),
    " | Stat:", signif(res_boot$T_D, 4), "\n")

# --- 5) Gruppen-Plot: complete (grau) vs incomplete (rot) -----------------------
group_A <- rowMeans(!is.na(X_obs)) >= 1
wide$group <- ifelse(group_A, "Complete", "Incomplete")

df_long <- wide |>
  pivot_longer(cols = all_of(slots), names_to = "slot", values_to = "temp") |>
  mutate(slot = factor(slot, levels = slots))

ggplot(df_long, aes(x = slot, y = temp, group = Datum, colour = group)) +
  geom_line(na.rm = TRUE) +
  scale_colour_manual(values = c("Complete" = "grey", "Incomplete" = "black")) +
  labs(x = "Uhrzeit", y = "Temperatur (°C)",
       title = "Tageskurven – Complete (grau) vs. Incomplete (rot)",
       colour = "Gruppe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# --- 5) Plots -------------------------------------------------------------------
# 5a) Mittelkurven μ_A und μ_B auf der Subdomain
times_sub <- slots[res_L2$idx]
df_mu <- data.frame(
  time = factor(times_sub, levels = slots),
  muA = res_L2$muA,
  muB = res_L2$muB
)

ggplot(df_mu, aes(time, group = 1)) +
  geom_line(aes(y = muA)) +
  geom_line(aes(y = muB), linetype = "dashed") +
  labs(x = "Uhrzeit", y = "Temperatur (°C)",
       title = "Mittelkurven (A: durchgezogen, B: gestrichelt) – δ = 0.43") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# 5b) Simultanes 95%-Konfidenzband für μ_A − μ_B
df_cb <- data.frame(
  time  = factor(slots[res_cb$idx], levels = slots),
  diff  = res_cb$diff,
  lower = res_cb$lower,
  upper = res_cb$upper
)

ggplot(df_cb, aes(time, diff, group = 1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Uhrzeit", y = "μ_A − μ_B (°C)",
       title = "Simultanes 95%-Band – δ = 0.43") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- Doppel-Plot: links Kurven, rechts Diff + 95%-Band ------------------------
library(ggplot2)
library(patchwork)

# Linkes Panel: beobachtete Tageskurven (Complete vs Incomplete)
p_left <- ggplot(df_long, aes(x = slot, y = temp, group = Datum, colour = group)) +
  geom_line(linewidth = 0.7, alpha = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c("Complete" = "grey60", "Incomplete" = "black")) +
  labs(
    title = "Tageskurven – Complete (grau) vs Incomplete (schwarz)",
    x = "Uhrzeit", y = "Temperatur (°C)", colour = "Gruppe"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Rechtes Panel: Differenz μ_A - μ_B mit simultanem 95%-Konfidenzband (aus res_cb)
df_cb <- data.frame(
  time  = factor(slots[res_cb$idx], levels = slots),
  diff  = res_cb$diff,
  lower = res_cb$lower,
  upper = res_cb$upper
)

p_right <- ggplot(df_cb, aes(time, diff, group = 1)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = lower), linetype = 2) +
  geom_line(aes(y = upper), linetype = 2) +
  labs(
    title = "Differenz der Mittelwerte (μA − μB) mit 95%-Band",
    x = "Uhrzeit", y = expression(mu[A] - mu[B] ~ "(°C)")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Nebeneinander anzeigen
p_left + p_right

# ===== Figure-7-Style: p-values vs number of grid points in I =====
# bin mir nicht sicher ob das so passt
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(2025)
B <- 2000                 # Bootstrap-Wiederholungen (für Paper: 10000)
m_grid <- 24:45           # wie im Paper

# 1) Verfügbarkeit je Slot zählen
avail <- colSums(!is.na(X_obs))
slot_idx <- seq_len(ncol(X_obs))      # 1..24 (00:30..23:30)

# 2) Für jedes m: Top-m Slots wählen, Tests auf Subdomain fahren
res_list <- lapply(m_grid, function(m) {
  idx_top <- order(-avail, slot_idx)[1:m] |> sort()     # Zeitreihen-Reihenfolge
  X_sub   <- X_obs[, idx_top, drop = FALSE]
  
  rb <- tfu_algo5_bootstrap(
    X_sub,
    delta_A   = 1,         # complete vs incomplete
    min_frac  = 0.10,      # Subdomain schon gewählt; Wert wird ignoriert/ok
    B         = B,
    seed      = 2025
  )
  
  data.frame(
    m    = m,
    p_L2 = rb$p_L2,
    p_D  = rb$p_D
  )
})

df_p <- bind_rows(res_list)

# 3) Long-Format + Plot
df_long <- df_p |>
  pivot_longer(cols = c(p_L2, p_D), names_to = "test", values_to = "p")

ggplot(df_long, aes(x = m, y = p, shape = test)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_manual(
    values = c(p_L2 = 1, p_D = 2),   # 1: offener Kreis, 2: offenes Dreieck
    labels = c(p_L2 = expression(T[mu*","*L^2]),
               p_D  = expression(T[mu*","*D]))
  ) +
  scale_x_continuous(breaks = m_grid) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "number of grid points in I", y = "p-values", shape = NULL) +
  theme_minimal(base_size = 12)
