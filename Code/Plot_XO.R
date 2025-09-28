library(ggplot2)
library(dplyr)
library(extrafont)
library(patchwork)

# Fonts laden
loadfonts(device = "win")

# Daten
t <- seq(0, 1, length.out = 200)
X <- 0.5 + sin(2 * pi * t) * 0.3
O <- ifelse(t <= 0.55 | (t >= 0.9 & t <= 1), 1, 0)
XO <- ifelse(O == 1, X, NA)

df <- data.frame(t, X, O, XO)

# Mittelpunkt der Lücke für "?"
x_missing <- mean(c(0.55, 0.9))
y_missing <- mean(range(X, na.rm = TRUE))

# --- Plot 1: X(t) ---
p1 <- ggplot(df, aes(x = t, y = XO)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  annotate("text", x = x_missing, y = y_missing, label = "?", 
           size = 22, color = "grey40", family = "Times New Roman") +
  labs(y = "X(t)", x = NULL) +
  theme_minimal(base_size = 20, base_family = "Times New Roman") +
  theme(
    axis.text.y = element_blank(),   # keine Y-Werte
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),   # keine X-Werte oben
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 25, family = "Times New Roman",
                                margin = margin(r = 15))  # mehr Platz rechts
  )

# --- Plot 2: O(t) ---
p2 <- ggplot(df, aes(x = t, y = O)) +
  geom_step(linewidth = 1.2) +
  scale_y_continuous(breaks = c(0, 1), limits = c(-0.1, 1.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(y = "O(t)", x = "t") +
  theme_minimal(base_size = 20, base_family = "Times New Roman") +
  theme(
    axis.text = element_text(size = 18, family = "Times New Roman"),
    axis.title.x = element_text(size = 22, family = "Times New Roman",
                                margin = margin(t = 15)),  # mehr Platz unten
    axis.title.y = element_text(size = 25, family = "Times New Roman",
                                margin = margin(r = 15)),  # mehr Platz rechts
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    # Achsenticks einschalten
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(6, "pt")
  )

# --- Kombinieren ---
p_final <- p1 / p2  

# Speichern
ggsave("Plots/XO_plot.png", p_final,
       width = 10, height = 4, dpi = 120)
