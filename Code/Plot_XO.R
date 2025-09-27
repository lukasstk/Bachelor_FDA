# Save plot as PNG
png("XO_plot.png", width = 800, height = 500, res = 120)

t <- seq(0, 1, length.out = 200)
X <- 0.5 + sin(2 * pi * t) * 0.3
O <- rep(0, length(t))
O[t <= 0.55] <- 1
O[t >= 0.9 & t <= 1] <- 1
XO <- ifelse(O == 1, X, NA)

par(mfrow = c(2,1), mar = c(3,5,1,1), cex.lab = 1.5, cex.axis = 1.3)

# Oberes Panel: X(t)
plot(t, X, type = "n", ylim = range(X, na.rm=TRUE),
     ylab = expression(X(t)), xlab = "", xaxt = "n", yaxt = "n")
lines(t, XO, lwd = 2)
text(0.74, mean(range(X, na.rm = TRUE)), "?", col = "grey40", cex = 4)

# Unteres Panel: O(t)
plot(t, O, type = "n", ylim = c(-0.1, 1.1),
     ylab = expression(O(t)), xlab = "t", yaxt = "n")
axis(2, at = c(0,1), labels = c("0","1"), las = 1, cex.axis = 1.3)

segments(0, 1, 0.55, 1, lwd = 5)
segments(0.9, 1, 1, 1, lwd = 5)
segments(0.55, 0, 0.9, 0, lwd = 5)

dev.off()
