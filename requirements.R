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

# --- Clean up temporary variable ---
rm(required_pkgs)
# ===============================================================
# Font setup: Import and register Times New Roman for plotting
# ===============================================================

# --- Detect system font directory automatically (cross-platform) ---
font_dir <- switch(
  Sys.info()[["sysname"]],
  "Windows" = "C:/Windows/Fonts",
  "Darwin"  = "/System/Library/Fonts",  # macOS
  "Linux"   = "/usr/share/fonts",        # Linux
  NULL
)

# --- Import Times New Roman font if not already registered ---
if (!"Times New Roman" %in% extrafont::fonts()) {
  message("Importing Times New Roman font...")
  extrafont::font_import(
    path = font_dir,
    pattern = "times.ttf",
    prompt = FALSE
  )
}

# --- Register fonts for plotting devices (PDF and system default) ---
extrafont::loadfonts(device = "pdf")
extrafont::loadfonts(device = ifelse(.Platform$OS.type == "windows", "win", "postscript"))

rm(font_dir)
# ===============================================================
# Load all package functions manually (development version)
# ===============================================================

# If the package is not loaded with `library(mcartest)`,
# all functions can be made available manually using:
devtools::load_all("cran_package/mcartest")

# ===============================================================
# Done: All required packages and fonts are ready for use
# ===============================================================
