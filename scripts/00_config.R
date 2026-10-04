# ------------------------------------------------------------------------------
# 00_config.R -- shared setup, sourced at the top of every analysis script.
#
# Run all scripts from the REPOSITORY ROOT (the folder containing data/ and
# scripts/), e.g.   Rscript scripts/01_statistical_tests_and_exploratory_plots.R
# ------------------------------------------------------------------------------

if (!dir.exists("data") || !dir.exists("scripts")) {
  stop("Run this script from the repository root (the folder that contains data/ and scripts/).")
}

DATA_DIR <- "data"
FIG_DIR  <- "figures"

# Random seed. The only stochastic element in these scripts is geom_jitter()
# (point placement in violin plots); statistics are deterministic.
SEED <- 42
set.seed(SEED)

# Number of distinct NITROPLASTS per time point used to normalize counts in
# 03_combined_vesicle_hat_figures.R (previously hardcoded in several places).
# These are nitroplast counts, not tomogram counts: several tomograms can cover the
# same nitroplast, so they differ from the tomogram counts in data/hat_counts_by_tomogram.csv.
N_NITRO_DAY     <- 8
N_NITRO_NIGHT   <- 33
N_NITRO_MORNING <- 28

# Read a CSV from data/. fileEncoding = "UTF-8-BOM" strips the byte-order mark
# that Excel adds; without it the first column is misnamed (e.g. "X...Thick")
# and df$Thick silently returns NULL.
read_data <- function(file) {
  read.csv(file.path(DATA_DIR, file), fileEncoding = "UTF-8-BOM")
}

# When run non-interactively (Rscript), write every plot of a script, in order,
# to figures/<name>.pdf. In RStudio / interactive sessions plots appear as usual.
open_fig_device <- function(name, width = 7, height = 5) {
  if (!interactive()) {
    dir.create(FIG_DIR, showWarnings = FALSE)
    pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)
  }
}
