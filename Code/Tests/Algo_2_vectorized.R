library(tidyfun)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Parameter definieren ----
n     <- 100    # Anzahl Pfade
grid_spacing <- 100    # Gitterpunkte
grid  <- seq(0, 1, length.out = grid_spacing)
grid_plot <- grid
set.seed(42)

# ---- 1) IDs und Gruppenverteilung (hälfte vollständig, hälfte mit Lücken) ----
ids <- paste0("ID", seq_len(n))
drop_ids <- sample(ids, size = n/2)          # mit Lücken
full_ids <- setdiff(ids, drop_ids)           # vollständig

# ---- 2) Brownian Motion Funktion ----
simulate_bm <- function(grid, mean_shift = 0) {
  dt  <- diff(grid)[1]
  inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
  c(0, cumsum(inc)) + mean_shift
}

# ---- 3) Simulation: BM-Pfade ----
bm_mat <- t(replicate(n, simulate_bm(grid)))
rownames(bm_mat) <- ids

make_O_mcar <- function(grid) {
  if (runif(1) < 0.5) {
    rep(1, length(grid))  # vollständig
  } else {
    U1 <- runif(1); U2 <- runif(1)
    L  <- min(U1, U2); U <- max(U1, U2)
    as.numeric(grid >= L & grid < U)  # nur im Intervall beobachtet
  }
}

# Hinweis aus deinem Code:
make_O_mnar <- function(x_row) {
  if (all(x_row > -1 & x_row < 2)) {
    return(rep(1, length(x_row)))  # vollständig beobachtet
  } else {
    return(as.numeric(x_row > -1 & x_row < 2))  # sonst nur in Bereich beobachtet
  }
}

create_O_mat <- function(bm_mat, grid, ids, mechanism = c("MCAR", "MNAR")) {
  mechanism <- match.arg(mechanism)
  if (mechanism == "MCAR") {
    O_mat <- t(sapply(ids, function(id) make_O_mcar(grid)))
  } else if (mechanism == "MNAR") {
    O_mat <- t(apply(bm_mat, 1, make_O_mnar))
  }
  rownames(O_mat) <- ids
  colnames(O_mat) <- round(grid, 3)
  O_mat
}

O_mat <- create_O_mat(bm_mat, grid, ids, mechanism = "MCAR")

# ---- 3) Indikatoren A und B ----
# A = komplett, B = unvollständig
I_A <- as.numeric(rowSums(O_mat) == grid_spacing)
I_B <- 1 - I_A

# ---- 4) Beobachtete Werte X * O ----
X_obs <- bm_mat * O_mat
X_obs[O_mat == 0] <- NA

# ---- Subdomain I: t mit min. 10% beob. Werten in beiden Gruppen ----
valid_t_idx <- which(
  colSums(O_mat * I_A) > 0.1 * n &
    colSums(O_mat * I_B) > 0.1 * n
)

# Grid/Daten auf gültige t einschränken
grid   <- grid[valid_t_idx]
bm_mat <- bm_mat[, valid_t_idx, drop = FALSE]
O_mat  <- O_mat[,  valid_t_idx, drop = FALSE]
X_obs  <- X_obs[,  valid_t_idx, drop = FALSE]

# ---- 5) p_hat & Gruppenmittel (vektorisiert, identische Logik) ----
IA <- as.numeric(I_A); IB <- 1 - IA
nA <- sum(IA); nB <- sum(IB)
m  <- length(grid)
dt <- diff(grid)[1]

pA_hat <- colSums(O_mat * IA) / n
pB_hat <- colSums(O_mat * IB) / n

muA_hat <- colSums(X_obs * IA, na.rm = TRUE) / (n * pA_hat)
muB_hat <- colSums(X_obs * IB, na.rm = TRUE) / (n * pB_hat)

# ---- 6) Teststatistik T_mu_D ----
T_mu_D <- sqrt(n) * max(abs(muA_hat - muB_hat))
cat("T_mu_D =", round(T_mu_D, 6), "\n")

# ---- 7) Kovarianzschätzer k(s,t) (voll vektorisiert, gleiche Formel) ----
# Zentrieren je Gruppe, unobserved -> 0, spaltenweise durch p_hat teilen
Xc <- X_obs
A_idx <- which(IA == 1L); B_idx <- which(IB == 1L)
if (nA) Xc[A_idx, ] <- sweep(X_obs[A_idx, , drop = FALSE], 2, muA_hat, `-`)
if (nB) Xc[B_idx, ] <- sweep(X_obs[B_idx, , drop = FALSE], 2, muB_hat, `-`)
Xc[is.na(Xc)] <- 0
Xc <- Xc * O_mat

XA <- if (nA) Xc[A_idx, , drop = FALSE] else matrix(0, 0, m)
XB <- if (nB) Xc[B_idx, , drop = FALSE] else matrix(0, 0, m)

scaleA <- ifelse(pA_hat > 0, 1 / pA_hat, 0)
scaleB <- ifelse(pB_hat > 0, 1 / pB_hat, 0)
if (nA) XA <- sweep(XA, 2, scaleA, "*")
if (nB) XB <- sweep(XB, 2, scaleB, "*")

Sum_A <- if (nA) crossprod(XA) else matrix(0, m, m)
Sum_B <- if (nB) crossprod(XB) else matrix(0, m, m)
corrected_K <- (Sum_A + Sum_B) / n

# ---- KL wie bei dir: Eigen(K), dann lam <- ev$values * dt; phi <- ev$vectors / sqrt(dt) ----
ev  <- eigen(corrected_K, symmetric = TRUE)
phi <- ev$vectors / sqrt(dt)
lam <- ev$values * dt
explained <- cumsum(lam) / sum(lam)
q <- which(explained >= 0.95)[1]
lam_q <- lam[1:q]
phi_q <- phi[, 1:q, drop = FALSE]

# ---- GP-Simulation (Algorithmus 2 / Bänder) vektorisiert ----
Bstar <- 5000
Z     <- matrix(rnorm(q * Bstar), nrow = q)         # q x B
lam_phi <- sweep(phi_q, 2, sqrt(lam_q), "*")        # m x q
gp_vals <- lam_phi %*% Z                            # m x B
W <- apply(abs(gp_vals), 2, max)

# Approximierter p-Wert (Supremum)
p_value <- (sum(W > T_mu_D) + 1) / (Bstar + 1)
cat("Approximate p-value =", round(p_value, 4), "\n")

# ---- Simultane Konfidenzbänder (Algorithmus 3) ----
alpha <- 0.05
q_alpha <- as.numeric(stats::quantile(W, probs = 1 - alpha, names = FALSE))

diff_mu <- muA_hat - muB_hat
ci_upper <- diff_mu + q_alpha / sqrt(n)
ci_lower <- diff_mu - q_alpha / sqrt(n)

# ---- Datenrahmen & Plot (unverändert) ----
mu_df <- data.frame(
  grid = grid,
  muA = muA_hat,
  muB = muB_hat,
  diff = diff_mu,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

ggplot(mu_df, aes(x = grid)) +
  geom_line(aes(y = muA, color = "muA")) +
  geom_line(aes(y = muB, color = "muB")) +
  geom_line(aes(y = diff, color = "muA - muB"), linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "black", alpha = 0.2) +
  labs(title = "Mittelkurven, Differenz und simultane Konfidenzbänder",
       x = "t", y = "Wert", color = "Kurve") +
  theme_minimal() +
  scale_color_manual(values = c("muA" = "blue", "muB" = "orange", "muA - muB" = "black"))
