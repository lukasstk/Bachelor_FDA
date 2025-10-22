# ===============================================================
# Example: Functional observation X(t) with missing region and corresponding O(t)
# ===============================================================
suppressPackageStartupMessages({
  library(extrafont)
  library(tidyverse)
  library(patchwork)
})

# --- Data setup ---
t <- seq(0, 1, length.out = 200)
X <- 0.5 + sin(2 * pi * t) * 0.3                          # underlying function
O <- ifelse(t <= 0.55 | (t >= 0.9 & t <= 1), 1, 0)        # observation indicator
XO <- ifelse(O == 1, X, NA)                               # observed part only
df <- data.frame(t, X, O, XO)

# Compute midpoint of the missing region for annotation ("?")
x_missing <- mean(c(0.55, 0.9))
y_missing <- mean(range(X, na.rm = TRUE))

# --- Plot 1: Observed curve X(t) with missing part ---
p1 <- ggplot(df, aes(x = t, y = XO)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  annotate("text", x = x_missing, y = y_missing, label = "?", 
           size = 22, color = "grey40", family = "Times New Roman") +
  labs(y = "X(t)", x = NULL) +
  theme_minimal(base_size = 20, base_family = "Times New Roman") +
  theme(
    axis.text.y  = element_blank(),      
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),      
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 25, family = "Times New Roman",
                                margin = margin(r = 15))  
  )

# --- Plot 2: Observation indicator O(t) ---
p2 <- ggplot(df, aes(x = t, y = O)) +
  geom_step(linewidth = 1.2) +
  scale_y_continuous(breaks = c(0, 1), limits = c(-0.1, 1.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(y = "O(t)", x = "t") +
  theme_minimal(base_size = 20, base_family = "Times New Roman") +
  theme(
    axis.text       = element_text(size = 18, family = "Times New Roman"),
    axis.title.x    = element_text(size = 22, family = "Times New Roman",
                                   margin = margin(t = 15)),  
    axis.title.y    = element_text(size = 25, family = "Times New Roman",
                                   margin = margin(r = 15)),  
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    axis.ticks      = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(6, "pt")                         
  )

# --- Combine both panels vertically ---
p_final <- p1 / p2  

# --- Save figure (optional) ---
# ggsave("Plots/x_o_plot.png", p_final,
#        width = 10, height = 4, dpi = 120)
