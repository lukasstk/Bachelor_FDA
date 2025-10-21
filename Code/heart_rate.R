# ===============================================================
# Heart Rate Data Example
# Apply asymptotic L2 and Supremum tests (with simultaneous bands)
# ===============================================================

library(extrafont)

# ===============================================================
# Load and prepare data
# ===============================================================
load("Data/heart_rate.RData")

rate_dt <- rate_dt %>% mutate(ti = as.numeric(ti))
fd_hr <- tf::tfd(rate_dt, arg = "ti", id = "id", value = "y0")

# Convert to wide format (rows = id, columns = time points)
wide_hr <- rate_dt |>
  tidyr::pivot_wider(
    id_cols     = id,
    names_from  = ti,
    values_from = y0
  )

# Extract numeric matrix (without id)
X <- wide_hr |>
  dplyr::select(-id) |>
  as.matrix()

row.names(X) <- as.character(wide_hr$id)

# ===============================================================
# Prepare grouping and plotting data (left panel)
# ===============================================================
# Grouping: complete (all observed) vs incomplete
group_A <- rowMeans(!is.na(X)) == 1
grp_df  <- tibble(
  id    = row.names(X),
  group = ifelse(group_A, "Complete", "Incomplete")
)

# Ensure consistent ID types before join
plot_df <- rate_dt %>%
  select(id, ti, y0) %>%
  mutate(id = as.character(id)) %>%
  inner_join(grp_df, by = "id")

# Time axis breaks and labels
hr_breaks <- c(20, 22, 24, 26)
hr_labels <- c("20:00", "22:00", "24:00", "02:00")

# ===============================================================
# Prepare data frame for plotting mean difference and bands
# ===============================================================
res_sup <- asym_mean_sup_test(fd = fd_hr, bands_only = TRUE)

bands <- data.frame(
  ti    = res_sup$bands$grid,                                                    
  diff  = tf::tf_evaluate(res_sup$estimate$mean_diff, arg = res_sup$bands$grid)[[1]],  
  lower = res_sup$bands$lower,                                                   
  upper = res_sup$bands$upper,                                                   
  row.names = NULL
) %>%
  arrange(ti)

# Identify gaps in the subdomain to draw separate line segments
unique_ti <- sort(unique(bands$ti))
base_step <- min(diff(unique_ti))
block <- c(0, cumsum(diff(bands$ti) > 1.5 * base_step))
bands$block <- block

# ===============================================================
# Visualization 
# ===============================================================

# --- Left panel: heart rate trajectories (complete vs incomplete) ---
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
    expand = c(0, 0)
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
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# --- Right panel: mean difference with simultaneous confidence bands ---
p_right <- ggplot(bands, aes(x = ti)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, fill = "grey50") +
  geom_line(aes(y = diff), linewidth = 1, color = "black") +
  geom_line(aes(y = lower), linetype = 2, color = "black") +
  geom_line(aes(y = upper), linetype = 2, color = "black") +
  scale_x_continuous(
    breaks = hr_breaks,
    labels = hr_labels,
    expand = expansion(mult = c(0, 0.02))
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
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# --- Combine both panels with spacer ---
heart_rate_plot <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# Save figure 
# ggsave(
#   filename = "Plots/heart_rate_plot.png",
#   plot     = heart_rate_plot,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )
