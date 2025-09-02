set.seed(12345)
S <- 100        # Anzahl Wiederholungen
alpha <- 0.05     # Signifikanzniveau
reject <- logical(S)

for (sim in seq_len(S)) {
  cat("Simulationslauf", sim, "\n")
  
  library(tidyfun)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # ---- Parameter definieren ----
  n     <- 500    # Anzahl Pfade
  grid_spacing     <- 100    # Gitterpunkte
  grid  <- seq(0, 1, length.out = grid_spacing)
  grid_plot <- grid
  
  # ---- 1) IDs und Gruppenverteilung (3 vollständig, 3 mit Lücken) ----
  ids <- paste0("ID", seq_len(n))
  drop_ids <- sample(ids, size = n/2)         # 3 mit Lücken
  full_ids <- setdiff(ids, drop_ids)        # 3 vollständig
  
  # ---- 2) Brownian Motion Funktion mit Mittelwertverschiebung ----
  simulate_bm <- function(grid, mean_shift = 0) {
    dt  <- diff(grid)[1]
    inc <- rnorm(length(grid) - 1, sd = sqrt(dt))
    c(0, cumsum(inc)) + mean_shift
  }
  
  # ---- 3) Simulation: 3 um +1 (vollständig), 3 um -1 (mit Lücken) ----
  bm_mat <- t(replicate(n, simulate_bm(grid)))
  rownames(bm_mat) <- ids
  
  make_O_mcar <- function(grid) {
    if (runif(1) < 0.5) {
      rep(1, length(grid))  # vollständig
    } else {
      U1 <- runif(1)
      U2 <- runif(1)
      L  <- min(U1, U2)
      U  <- max(U1, U2)
      as.numeric(grid >= L & grid < U)  # nur im Intervall beobachtet
    }
  }
  
  # bekomm pval < 0.05 nur bei hohem n, mit n = 100 wie im paper funktionierts bei mir nicht
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
    return(O_mat)
  }
  
  
  O_mat <- create_O_mat(bm_mat, grid, ids, mechanism = "MCAR")
  
  # ---- 3) Indikatoren A und B ----
  # A = komplett, B = unvollständig
  I_A <- as.numeric(rowSums(O_mat) == grid_spacing)
  I_B <- 1 - I_A
  
  # ---- 4) Beobachtete Werte X * O ----
  #    bm_mat enthält die ungestörten Pfade
  X_obs <- bm_mat * O_mat
  X_obs[O_mat == 0] <- NA
  
  # ---- Subdomain I aus dem Paper: Nur t mit min. 10% beobachteten Werten in beiden Gruppen ----
  valid_t_idx <- which(
    colSums(O_mat * I_A) > 0.1 * n &
      colSums(O_mat * I_B) > 0.1 * n
  )
  
  # Grid, Daten und Matrizen auf gültige t einschränken
  grid <- grid[valid_t_idx]
  bm_mat <- bm_mat[, valid_t_idx]
  O_mat <- O_mat[, valid_t_idx]
  X_obs <- X_obs[, valid_t_idx]
  
  # ---- 5) Schätzer für p_A, mu_A und analog für B ----
  
  # Anteil beobachteter Werte (pA_hat, pB_hat)
  pA_hat <- colSums(O_mat * I_A) / n
  pB_hat <- colSums(O_mat * I_B) / n
  
  # Gruppenspezifische Mittelwerte, NA ignorieren
  muA_hat <- colSums(X_obs * I_A, na.rm = TRUE) / (n * pA_hat) 
  
  # z.B. bei Punkt 0.25: soll das ca. 1.9 sein oder 0.6580639 -> sum(X_obs["ID1",25])/(n * pA_hat[25])
  muB_hat <- colSums(X_obs * I_B, na.rm = TRUE) / (n * pB_hat)
  
  # ---- 6) Teststatistik T_mu_D ----
  T_mu_D <- sqrt(n) * max(abs(muA_hat - muB_hat))
  cat("T_mu_D =", round(T_mu_D, 6), "\n")
  
  # ---- 7) Kovarianzschätzer k(s,t) nach Paper (Formel 2) ----
  
  X_tilde <- matrix(NA, nrow = n, ncol = length(grid))
  for (i in seq_len(n)) {
    if (I_A[i] == 1) {
      X_tilde[i, ] <- X_obs[i, ] - muA_hat
    } else {
      X_tilde[i, ] <- X_obs[i, ] - muB_hat
    }
  }
  
  # Initialisiere Kovarianzmatrix
  corrected_K <- matrix(0, nrow = length(grid), ncol = length(grid))
  
  # k(s, t) schätzen nach Paper
  for (s in seq_len(length(grid))) {
    for (t in seq_len(length(grid))) {
      # Zählerterme vorbereiten
      term_A <- (X_tilde[, s] * X_tilde[, t] * O_mat[, s] * O_mat[, t] * I_A)
      term_B <- (X_tilde[, s] * X_tilde[, t] * O_mat[, s] * O_mat[, t] * I_B)
      
      # Nenner vorbereiten
      denom_A <- pA_hat[s] * pA_hat[t]
      denom_B <- pB_hat[s] * pB_hat[t]
      
      # Kovarianzanteile, falls Nenner > 0
      part_A <- if (denom_A > 0) sum(term_A / denom_A, na.rm = TRUE) else 0
      part_B <- if (denom_B > 0) sum(term_B / denom_B, na.rm = TRUE) else 0
      
      # Finaler Eintrag
      corrected_K[s, t] <- (part_A + part_B) / n
    }
  }
  corrected_K <- corrected_K * diff(grid)[1]
  ev <- eigen(corrected_K, symmetric = TRUE)
  phi <- ev$vectors / sqrt(diff(grid)[1])
  lam <- ev$values 
  explained <- cumsum(lam) / sum(lam)
  q <- which(explained >= 0.95)[1]
  lam_q <- lam[1:q]
  phi_q <- phi[,1:q]
  Bstar <- 5000
  Z_all <- matrix(NA, nrow = Bstar, ncol = q)
  W <- numeric(Bstar)
  
  for(i in seq(Bstar)){
    Z_all[i, ] <- rnorm(length(lam_q))
    weighted_sum <- phi_q %*% (sqrt(lam_q) * Z_all[i, ])
    W[i] <- max(abs(weighted_sum))
  }
  
  # # or: macht genau das selbe
  #   for(i in seq(Bstar)){
  #     Z_all[i, ] <- rnorm(length(lam_q))
  #     lam_phi <- sweep(phi_q, MARGIN = 2, sqrt(lam_q), FUN = "*")
  #     weighted_sum <- rowSums(sweep(lam_phi, MARGIN = 2, Z_all[i, ], FUN = "*"))
  #     W[i] <- max(abs(weighted_sum))
  #   }
  
  # # or:
  # for (i in seq(Bstar)) {
  #   Z <- rnorm(q)
  #   Z_all[i, ] <- Z
  # 
  #   # Initialisiere Vektor für gewichtete Summe über alle j
  #   weighted_sum <- rep(0, length(grid))  # Länge = Anzahl Gitterpunkte
  # 
  #   for (j in seq_len(q)) {
  #     weighted_sum <- weighted_sum + Z[j] * sqrt(lam_q[j]) * phi_q[, j]
  #   }
  # 
  #   W[i] <- max(abs(weighted_sum))
  # }
  
  # Approximierter p-Wert
  p_value <- (sum(W > T_mu_D) + 1) / (Bstar + 1)
  
  # Ausgabe
  cat("Approximate p-value =", round(p_value, 4), "\n")
  
  
  reject[sim] <- (p_value < alpha)
}

type1_error <- mean(reject)
cat("Geschätzter Typ-I-Fehler (alpha =", alpha, ") =", round(type1_error, 3), "\n")
