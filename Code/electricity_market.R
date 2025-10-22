# ===============================================================
# Functional Electricity Market Data
# Apply asymptotic L2 and Supremum tests on log-bid curves
# ===============================================================
suppressPackageStartupMessages({
  library(extrafont)
  library(tidyverse)
  library(patchwork)
})

# ===============================================================
# Load and prepare data
# ===============================================================
data <- as.matrix(read.csv("Data/logbidcurves.csv", header = TRUE, check.names = FALSE))
grid_full <- seq(1700, 2500, length.out = 161)
X_full <- if (ncol(data) == length(grid_full)) data else t(data)

# ===============================================================
# Run asymptotic Supremum test (with bands)
# ===============================================================
set.seed(2025)
fd_full <- tf::tfd(X_full, arg = grid_full)
observed_ratio <- 0.43

res_sup <- asym_mean_sup_test(
  fd = fd_full,
  observed_ratio = observed_ratio,   
  compute_bands = TRUE,
  seed = 42
)

# ===============================================================
# Prepare data for plotting
# ===============================================================
lab_A <- sprintf("obs \u2265 %.0f%%", observed_ratio * 100)  # "obs ≥ 43%"
lab_B <- sprintf("obs < %.0f%%",  observed_ratio * 100)      # "obs < 43%"

group_A    <- rowMeans(!is.na(X_full)) >= observed_ratio
group_labs <- ifelse(group_A, lab_A, lab_B)

df_long <- as.data.frame(X_full)
colnames(df_long) <- as.character(grid_full)
df_long$id <- sprintf("curve_%05d", seq_len(nrow(df_long)))
df_long$group <- factor(group_labs, levels = c(lab_A, lab_B))    

df_long <- tidyr::pivot_longer(
  df_long, cols = -c(id, group),
  names_to = "demand", values_to = "logprice"
) |> 
  dplyr::mutate(demand = as.numeric(demand))

bands <- data.frame(
  t     = res_sup$bands$grid,
  diff  = tf::tf_evaluate(res_sup$estimate$mean_diff, arg = res_sup$bands$grid)[[1]],
  lower = res_sup$bands$lower,
  upper = res_sup$bands$upper
)
rownames(bands) <- NULL

# ===============================================================
# Visualization
# ===============================================================
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
    expand  = c(0, 0)
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
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

# --- Right panel: mean difference with confidence bands ---
p_right <- ggplot(bands, aes(x = t)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, fill = "grey50") +
  geom_line(aes(y = diff), linewidth = 1, color = "black") +
  geom_line(aes(y = lower), linetype = 2, color = "black") +
  geom_line(aes(y = upper), linetype = 2, color = "black") +
  scale_x_continuous(
    breaks = c(1800, 2000, 2200, 2400),
    expand = c(0, 0)
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
    axis.title.x      = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.ticks.length = unit(3, "pt")
  )

electricity_market_plot <- p_left + plot_spacer() + p_right +
  plot_layout(widths = c(1, 0.1, 1))

# # Save figure
# ggsave(
#   filename = "Plots/electricity_market_plot.png",
#   plot     = electricity_market_plot,
#   width    = 10,
#   height   = 4,
#   dpi      = 300
# )
