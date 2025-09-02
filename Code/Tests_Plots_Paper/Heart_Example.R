# ===================== Setup & Packages =====================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyfun)   # tfd / tf_integrate
  library(patchwork)
  # Optional für Bootstrap:
  library(doParallel); library(doRNG); library(foreach)
})

# ---- Deine neuen Funktionen müssen im Environment sein ----
# asym_mean_L2_test(), asym_mean_sup_test(), boot_mean_test(), ...

# ===================== 1) Daten laden ======================
load("Data/heart_rate.RData")
# Erwartet: rate_dt mit Spalten id, ti (Stunden 20–26), y0 (bpm)

# ===================== 2) tfd-Objekt & Wide-Matrix ============================
stopifnot(all(c("id","ti","y0") %in% names(rate_dt)))
rate_dt <- rate_dt %>% mutate(ti = as.numeric(ti))

fd_hr <- tf::tfd(rate_dt, arg = "ti", id = "id", value = "y0")

# Wide: Zeilen = id, Spalten = Zeitpunkte
wide_hr <- rate_dt %>%
  select(id, ti, y0) %>%
  distinct() %>%
  arrange(id, ti) %>%
  tidyr::pivot_wider(names_from = ti, values_from = y0)

# Matrix X_obs (ohne id), Spalten numerisch sortieren
X_obs <- wide_hr %>% select(-id) %>% as.matrix()
col_order <- order(as.numeric(colnames(X_obs)))
X_obs <- X_obs[, col_order, drop = FALSE]
colnames(X_obs) <- as.character(sort(as.numeric(colnames(X_obs))))
row.names(X_obs) <- as.character(wide_hr$id)  # IDs als character

# ===================== 3) Tests (asymptotisch + Bänder) ======================
set.seed(2025)

res_L2 <- asym_mean_L2_test(
  fd = fd_hr,
  observed_ratio = 1,
  min_frac = 0.10,
  fve = 0.99,
  B = 5000,
  seed = 2025
)

res_sup <- asym_mean_sup_test(
  fd = fd_hr,
  observed_ratio = 1,
  min_frac = 0.10,
  fve = 0.99,
  B = 5000,
  compute_bands = TRUE,
  seed = 2025
)

cat("\n=== HEART RATE: Ergebnisse (complete vs incomplete) ===\n")
print(res_L2)
print(res_sup)

# ===================== 4) Plot: Kurvenschar (links) ===========================
# Gruppierung wie in den Tests: complete (alle Zeitpunkte belegt) vs incomplete
group_A <- rowMeans(!is.na(X_obs)) == 1
grp_df  <- tibble(
  id    = row.names(X_obs),
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

p_left <- ggplot(plot_df, aes(x = ti, y = y0, group = id, colour = group)) +
  geom_line(linewidth = 0.5, alpha = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c("Complete" = "grey70", "Incomplete" = "black")) +
  scale_x_continuous(breaks = hr_breaks, labels = hr_labels, expand = expansion(0)) +
  labs(x = "time", y = "heart rate (bpm)", colour = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# ===================== 5) Plot: Mean-Diff + Bänder (rechts) ===================
df_band <- data.frame(
  ti    = res_sup$grid,
  diff  = res_sup$diff,
  lower = res_sup$lower,
  upper = res_sup$upper
) %>% arrange(ti)

# Lücken in der Subdomäne als getrennte Linienblöcke zeichnen
unique_ti <- sort(unique(df_band$ti))
base_step <- min(diff(unique_ti))
block <- c(0, cumsum(diff(df_band$ti) > 1.5 * base_step))
df_band$block <- block

p_right <- ggplot(df_band, aes(x = ti)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_line(aes(y = lower, group = block), linetype = "dashed") +
  geom_line(aes(y = upper, group = block), linetype = "dashed") +
  geom_line(aes(y = diff,  group = block), linewidth = 0.9) +
  scale_x_continuous(breaks = hr_breaks, labels = hr_labels, expand = expansion(0)) +
  scale_y_continuous(limits = c(-40, 40), breaks = seq(-40, 40, by = 20)) +
  labs(x = "time", y = "difference in means") +
  theme_classic(base_size = 14)

# ===================== 6) Anzeigen (Seite an Seite) ===========================
p_left + p_right

# ===================== (Optional) Bootstrap-Variante ===========================
# res_boot <- boot_mean_test(
#   fd = fd_hr, observed_ratio = 1, min_frac = 0.10, alpha = 0.05,
#   stat = "D", B = 10000, seed = 2025,
#   parallel = TRUE, ncpus = max(1L, parallel::detectCores() - 1L),
#   manage_backend = "auto", compute_bands = TRUE
# )
# print(res_boot)
# # Plot analog mit res_boot$grid, $diff, $lower, $upper
