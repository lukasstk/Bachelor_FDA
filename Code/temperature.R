# ===============================================================
# Temperature Data (Graz)
# Bootstrap-based L2 and Supremum tests with simultaneous bands
# ===============================================================

library(extrafont)
library(lubridate)

# ===============================================================
# Load and prepare data
# ===============================================================
load("Data/temp_graz.rda")

# Convert to proper datetime, extract date and half-hour slots
df <- temp_graz %>%
  mutate(
    time = as.POSIXct(time, tz = "UTC"),
    date = as.Date(time),
    slot = sprintf("%02d:%02d", hour(time), minute(time))
  ) %>%
  arrange(time)

# Define 48 half-hour slots per day: 00:00..23:30
slots48 <- sprintf("%02d:%02d", rep(0:23, each = 2), rep(c(0, 30), times = 24))

# Wide format: rows = days, columns = 48 slots, values = east temperature
wide <- df %>%
  filter(slot %in% slots48) %>%
  mutate(slot = factor(slot, levels = slots48)) %>%
  select(date, slot, temp = temp_east) %>%
  distinct() %>%
  pivot_wider(names_from = slot, values_from = temp) %>%
  arrange(date)

# Temperature matrix
X <- wide %>% select(all_of(slots48)) %>% as.matrix()

# Check data coverage
obs_per_day <- rowMeans(!is.na(X))
cat("Total days:", nrow(X),
    " | Median observed per day:", signif(median(obs_per_day), 3), "\n")

slot_to_hours <- function(s) {
  hh <- as.numeric(substr(s, 1, 2))
  mm <- as.numeric(substr(s, 4, 5))
  hh + mm / 60
}

df_long <- wide %>%
  pivot_longer(cols = all_of(slots48), names_to = "slot", values_to = "temp") %>%
  mutate(hour_num = slot_to_hours(slot))

fd <- tf::tfd(df_long, arg = "hour_num", id = "date", value = "temp")

# ===============================================================
# Prepare data for plotting
# ===============================================================

# Assign groups for plotting
wide$group <- ifelse(group_A, "Complete", "Incomplete")

# Prepare long-format data for left panel
df_long_plot <- wide %>%
  tidyr::pivot_longer(cols = all_of(slots48), names_to = "slot", values_to = "temp") %>%
  mutate(slot = factor(slot, levels = slots48))

# --- Prepare data frame for plotting mean difference and bootstrap bands ---
grid_sub  <- res_sup$bands$grid
grid_full <- slot_to_hours(slots48)

# Identify positions of subgrid points in the full grid
idx_out <- vapply(grid_sub, function(g) which.min(abs(grid_full - g)), integer(1))

# Sort and define continuity blocks (for plotting gaps)
o          <- order(idx_out)
idx_sorted <- idx_out[o]
block      <- c(0, cumsum(diff(idx_sorted) > 1))

res_sup <- boot_mean_test(
  fd = fd, stat = "D", compute_bands = TRUE)

# Extract mean difference and bands
bands <- data.frame(
  time  = factor(slots48[idx_sorted], levels = slots48),
  diff  = as.numeric(as.matrix(res_sup$estimate$mean_diff))[o],
  lower = as.numeric(res_sup$bands$lower[o]),
  upper = as.numeric(res_sup$bands$upper[o]),
  block = block
)

# ===============================================================
# Visualization 
# ===============================================================

# --- Left panel: daily temperature curves ---
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
  labs(x = "time", y = expression("temperature (" * degree * C * ")")) +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    legend.position   = "none",
    text              = element_text(family = "Times New Roman"),
    axis.title        = element_text(size = 20),
    axis.text         = element_text(size = 16),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black", linewidth = 0.5),
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# --- Right panel: mean difference and bootstrap confidence bands ---
p_right <- ggplot(bands, aes(x = time)) +
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
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(x = "time", y = expression("difference in means (" * degree * C * ")")) +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    text              = element_text(family = "Times New Roman"),
    axis.title        = element_text(size = 18),
    axis.text         = element_text(size = 16),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "black", linewidth = 0.5),
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# Combine both panels
temperature_plot <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# # Save plot
# ggsave(
#   filename = "Plots/temperature_plot.png",
#   plot     = temperature_plot,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )

# ============================================================================
# Plot: Rejection probabilities vs. grid size
# Tests the power of T_{μ,L²} and T_{μ,D} under increasing number of gridpoints
# ===========================================================================
set.seed(2025)
m_grid <- c(24, 26, 28, 31, 33, 36, 38, 40, 43, 45)

# Groups fixed on the full grid to ensure consistent assignment
groups_fixed <- rowMeans(!is.na(X)) == 1

# Rank slots by overall availability
avail    <- colSums(!is.na(X))
slot_idx <- seq_len(ncol(X))

# --- Compute bootstrap p-values for each subgrid size ---
res_list <- pblapply(m_grid, function(m) {
  idx_top <- order(-avail, slot_idx)[1:m] |> sort()
  X_sub   <- X[, idx_top, drop = FALSE]
  
  rb_L2 <- boot_mean_test(X = X_sub, groups = groups_fixed,
                          stat = "L2", compute_bands = FALSE)
  
  rb_D  <- boot_mean_test(X = X_sub,  groups = groups_fixed,
                          stat = "D", compute_bands = FALSE)
  
  data.frame(
    m    = m,
    p_L2 = as.numeric(rb_L2$p.value),
    p_D  = as.numeric(rb_D$p.value)
  )
})

# Combine and reshape results
df_p <- dplyr::bind_rows(res_list) |>
  tidyr::pivot_longer(cols = c(p_L2, p_D), names_to = "test", values_to = "p")

df_table <- df_p |>
  tidyr::pivot_wider(names_from = test, values_from = p) |>
  dplyr::rename(`T[mu,L2]` = p_L2, `T[mu,D]` = p_D)

print(df_table, n = Inf)

temperature_rej_probs_plot <- ggplot(df_p, aes(x = m, y = p, shape = test, color = test)) +
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
    limits = c(min(m_grid) - 1, max(m_grid) + 1),
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
  guides(shape = guide_legend(override.aes = list(size = 4)))

print(temperature_rej_probs_plot)

# # Save figure
# ggsave(
#   filename = "Plots/temperature_rej_probs_plot.png",
#   plot     = temperature_rej_probs_plot,
#   width    = 10,
#   height   = 5,
#   dpi      = 300
# )
