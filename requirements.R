# ===============================================================
# Setup: Install, load, and configure required R packages and fonts
# ===============================================================

# --- Install 'pak' if not already installed ---
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# --- Define required packages ---
required_pkgs <- c(
  "dplyr", "extrafont", "patchwork", "pbapply", "tf",
  "tidyverse", "tidyr", "lubridate", "checkmate",
  "foreach", "doParallel", "doRNG", "testthat"
)

# --- Install all missing or outdated packages (pak skips already up-to-date ones) ---
pak::pak(required_pkgs)

# --- Load all packages quietly (no startup messages or warnings) ---
invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# ===============================================================
# Font setup: Import and register Times New Roman for plotting
# ===============================================================

# --- Import only Times New Roman font files from the Windows font directory ---
extrafont::font_import(
  path = "C:/Windows/Fonts",
  pattern = "times.ttf",
  prompt = FALSE
)

# --- Register fonts for Windows devices (e.g., PDF, plots) ---
extrafont::loadfonts(device = "win")
extrafont::loadfonts(device = "pdf")

# ===============================================================
# Done: All required packages and fonts are ready for use
# ===============================================================
