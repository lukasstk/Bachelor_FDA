# --- Setup ---------------------------------------------------------------
library(ggplot2)
library(extrafont)
loadfonts(device = "win")

set.seed(42)
t <- seq(0, 1, length.out = 400)

# Glatte Grundfunktion X(t)
shape <- function(t, phase = 0, amp = 1, offset = 0) {
  amp * (0.6 * exp(-((t - 0.6)^2) / (2 * 0.07^2)) + 0.15 * sin(2 * pi * (t + phase))) + offset
}

# Zufällige Missing-Intervals in ~0.2-Stücken erzeugen
mask_random_chunks <- function(t, observed_prop, chunk_prop = 0.2) {
  n <- length(t)
  miss <- max(0, 1 - observed_prop)
  if (miss <= 0) return(rep(TRUE, n))
  
  K <- floor(miss / chunk_prop)
  leftover <- miss - K * chunk_prop
  lengths <- c(rep(chunk_prop, K), if (leftover > 1e-6) leftover else numeric(0))
  
  observed <- rep(TRUE, n)
  taken <- matrix(0, nrow = 0, ncol = 2)
  
  for (L in lengths) {
    len_idx <- max(1, round(L * (n - 1)))
    tries <- 0
    placed <- FALSE
    while (!placed && tries < 2000) {
      s_idx <- sample.int(n - len_idx, 1)
      e_idx <- s_idx + len_idx
      # Check Overlap
      if (nrow(taken) == 0 ||
          all(e_idx < taken[,1] | s_idx > taken[,2])) {
        observed[s_idx:e_idx] <- FALSE
        taken <- rbind(taken, c(s_idx, e_idx))
        placed <- TRUE
      }
      tries <- tries + 1
    }
  }
  observed
}

gen_curve <- function(coverage, id, panel, phase = 0, amp = 1, offset = 0) {
  x <- shape(t, phase, amp, offset)
  obs_mask <- mask_random_chunks(t, observed_prop = coverage, chunk_prop = 0.2)
  x[!obs_mask] <- NA
  data.frame(t = t, x = x, id = id, panel = panel)
}

# Kleine Variation zwischen Kurven
jit <- function(v, s = 0.02) v + runif(1, -s, s)

# Gemeinsames Theme: größere Labels, Kasten pro Panel, keine Tick-Werte
base_theme <- theme_minimal(base_size = 20, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    strip.text = element_text(size = 30, family = "Times New Roman",
                              margin = margin(t = 10)),   # Abstand der Facettenlabels (A,B) nach unten
    strip.placement = "outside",                          # Labels nach außen legen
    panel.spacing.x = unit(12, "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 30, family = "Times New Roman",
                                margin = margin(r = 20))   # weiter weg von der Achse
  )



# --- Example 1: complete (A) vs incomplete (B) --------------------------
A_cov1 <- c(1.00, 1.00, 1.00)
B_cov1 <- c(0.90, 0.85, 0.80)

df_ex1 <- rbind(
  gen_curve(A_cov1[1], "A1", "A", phase = jit(0.00), offset = -0.05),
  gen_curve(A_cov1[2], "A2", "A", phase = jit(0.02), offset =  0.00),
  gen_curve(A_cov1[3], "A3", "A", phase = jit(-0.02), offset =  0.05),
  gen_curve(B_cov1[1], "B1", "B", phase = jit(0.00), offset = -0.05),
  gen_curve(B_cov1[2], "B2", "B", phase = jit(0.02), offset =  0.00),
  gen_curve(B_cov1[3], "B3", "B", phase = jit(-0.02), offset =  0.05)
)

p_ex1 <- ggplot(df_ex1, aes(t, x, group = id, color = panel)) +
  geom_line(linewidth = 1.2, lineend = "round", na.rm = TRUE) +
  scale_color_manual(values = c(A = "#bdbdbd", B = "#000000"), guide = "none") +
  facet_wrap(~panel, nrow = 1, strip.position = "bottom") +
  labs(y = expression(X(t))) +   
  base_theme

print(p_ex1)


# --- Example 2: threshold delta = 0.9 -----------------------------------
delta <- 0.9
A_cov2 <- c(1.00, 0.95, 0.90)  # >= delta
B_cov2 <- c(0.75, 0.50, 0.40)  # < delta

df_ex2 <- rbind(
  gen_curve(A_cov2[1], "A1", "A", phase = jit(0.00), offset = -0.05),
  gen_curve(A_cov2[2], "A2", "A", phase = jit(0.02), offset =  0.00),
  gen_curve(A_cov2[3], "A3", "A", phase = jit(-0.02), offset =  0.05),
  gen_curve(B_cov2[1], "B1", "B", phase = jit(0.00), offset = -0.05),
  gen_curve(B_cov2[2], "B2", "B", phase = jit(0.02), offset =  0.00),
  gen_curve(B_cov2[3], "B3", "B", phase = jit(-0.02), offset =  0.05)
)

p_ex2 <- ggplot(df_ex2, aes(t, x, group = id, color = panel)) +
  geom_line(linewidth = 1.2, lineend = "round", na.rm = TRUE) +
  scale_color_manual(values = c(A = "#bdbdbd", B = "#000000"), guide = "none") +
  facet_wrap(~panel, nrow = 1, strip.position = "bottom") +
  labs(y = expression(X(t))) +   
  base_theme

print(p_ex2)


# Plot 1 speichern
ggsave("Plots/grouping_partition_1.png", plot = p_ex1, width = 10, height = 4, dpi = 300)

# Plot 2 speichern
ggsave("Plots/grouping_partition_2.png", plot = p_ex2, width = 10, height = 4, dpi = 300)