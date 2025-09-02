# Pakete
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyfun)   # nutzt tf::tfd / tf_integrate
library(patchwork)
# Für parallelen Bootstrap (die Funktionen nutzen foreach + doRNG intern)
library(doParallel)
library(doRNG)
library(foreach)

# -------- 1) Einlesen ----------------------------------------------------------
file <- "Code/Tests/luis-daten.csv"

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

# -------- 2) Tageskurven (24 Halbstunden-Slots: 00:30..23:30) ------------------
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
obs_per_day <- rowMeans(!is.na(X_obs))
cat("Tage gesamt:", nrow(X_obs),
    " | Median beobachtet/Tag (Anteil):", signif(median(obs_per_day), 3), "\n")

# -------- 2b) fd-Objekt (tfd_irreg) --------------------------------------------
slot_to_hour <- function(s) as.numeric(substr(s, 1, 2)) + 0.5

df_long <- wide |>
  pivot_longer(cols = all_of(slots), names_to = "slot", values_to = "temp") |>
  mutate(hour_num = slot_to_hour(slot))

fd <- tf::tfd(df_long, arg = "hour_num", id = "Datum", value = "temp")

# -------- 3) Analyse mit observed_ratio = 1 -----------------------------------
set.seed(2025)

# Algo 1: asymptotischer L2-Test (unverändert)
res_L2 <- asym_mean_L2_test(fd = fd, observed_ratio = 1, min_frac = 0.10, seed = 2025)

# Algo 2+3: Supremum-Test + asymptotische simultane Bänder (NEU: kombiniert)
res_asym <- asym_sup_and_bands(fd = fd, observed_ratio = 1, min_frac = 0.10, seed = 2025)
res_sup  <- res_asym          # für Kompatibilität der Variablennamen unten
res_cb   <- res_asym          # Bänder kommen nun aus dem Supremum-Objekt

# Algo 5+7: Bootstrap-Tests (L2 & D) + Bootstrap-Bänder (NEU: kombiniert)
res_boot_all <- boot_mean_test_and_bands(
  X_obs       = X_obs,
  min_frac    = 0.10,
  B           = 2000,
  seed        = 2025,
  return_boot = TRUE,
  parallel    = TRUE
)
# einzelne htest-Objekte für bequemes Printen
res_boot_L2 <- res_boot_all$test_L2
res_boot_D  <- res_boot_all$test_D

# -------- 4) Ergebnisse ausgeben -----------------------------------------------
n_total <- nrow(X_obs)
group_A <- rowMeans(!is.na(X_obs)) == 1
n_A <- sum(group_A); n_B <- n_total - n_A

cat("\n==== Ergebnisse (observed_ratio = 1 – complete vs incomplete) ====\n")
cat("n =", n_total, "\n")
cat("Größe Gruppe A/B:", n_A, "/", n_B, "\n")

cat("Subdomain-Punkte (Bänder/asym):", as.integer(res_cb$parameter["m"]),
    " | Fallback:", ifelse(is.null(res_cb$fallback), "none", res_cb$fallback), "\n")

# direkte Ausgabe von htest-kompatiblen Objekten
print(res_L2)
print(res_sup)

cat("\n==== Bootstrap-Ergebnisse (foreach + doRNG) ====\n")
cat("L2:\n"); print(res_boot_L2)
cat("\nSupremum (D):\n"); print(res_boot_D)

# -------- 5) Gruppen-Plot ------------------------------------------------------
wide$group <- ifelse(group_A, "Complete", "Incomplete")

df_long_plot <- wide |>
  pivot_longer(cols = all_of(slots), names_to = "slot", values_to = "temp") |>
  mutate(slot = factor(slot, levels = slots))

p_groups <- ggplot(df_long_plot, aes(x = slot, y = temp, group = Datum, colour = group)) +
  geom_line(linewidth = 0.7, alpha = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c("Complete" = "grey60", "Incomplete" = "red")) +
  labs(title = "Tageskurven – Complete (grau) vs Incomplete (rot)",
       x = "Uhrzeit", y = "Temperatur (°C)", colour = "Gruppe") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# -------- 5a) Mittelkurven -----------------------------------------------------
times_sub <- slots[res_cb$idx]
df_mu <- data.frame(time = factor(times_sub, levels = slots),
                    muA = res_cb$muA, muB = res_cb$muB)

p_mu <- ggplot(df_mu, aes(time, group = 1)) +
  geom_line(aes(y = muA)) +
  geom_line(aes(y = muB), linetype = "dashed") +
  labs(x = "Uhrzeit", y = "Temperatur (°C)",
       title = "Mittelkurven (A: durchgezogen, B: gestrichelt)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# -------- 5b) Konfidenzbänder -------------------------------------------------
df_cb <- data.frame(time  = factor(slots[res_cb$idx], levels = slots),
                    diff  = res_cb$diff,
                    lower = res_cb$lower,
                    upper = res_cb$upper)

p_band <- ggplot(df_cb, aes(time, diff, group = 1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Uhrzeit", y = expression(mu[A] - mu[B]),
       title = "Simultanes 95%-Band") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# -------- Doppel-Plot ----------------------------------------------------------
p_left  <- p_groups
p_right <- ggplot(df_cb, aes(time, diff, group = 1)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = lower), linetype = 2) +
  geom_line(aes(y = upper), linetype = 2) +
  labs(title = "Differenz der Mittelwerte (μA − μB) mit 95%-Band",
       x = "Uhrzeit", y = expression(mu[A] - mu[B])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_left + p_right

# -------- 6) Figure-7-Style: p-values vs number of grid points ----------------
set.seed(2025)
B <- 5000
m_grid <- 6:24

avail <- colSums(!is.na(X_obs))
slot_idx <- seq_len(ncol(X_obs))

res_list <- lapply(m_grid, function(m) {
  idx_top <- order(-avail, slot_idx)[1:m] |> sort()
  X_sub   <- X_obs[, idx_top, drop = FALSE]
  
  rb_all <- boot_mean_test_and_bands(
    X_obs       = X_sub,
    B           = B,
    min_frac    = 0.10,
    seed        = 2025,
    return_boot = FALSE,
    parallel    = TRUE
  )
  
  data.frame(
    m    = m,
    p_L2 = as.numeric(rb_all$test_L2$p.value),
    p_D  = as.numeric(rb_all$test_D$p.value)
  )
})

df_p <- bind_rows(res_list)

df_long_p <- df_p |>
  pivot_longer(cols = c(p_L2, p_D), names_to = "test", values_to = "p")

ggplot(df_long_p, aes(x = m, y = p, shape = test)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_manual(
    values = c(p_L2 = 1, p_D = 2),
    labels = c(p_L2 = expression(T[mu*","*L^2]),
               p_D  = expression(T[mu*","*D]))
  ) +
  scale_x_continuous(breaks = m_grid) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "number of grid points in I", y = "p-values", shape = NULL)
