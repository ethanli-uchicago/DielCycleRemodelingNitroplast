# Diel remodeling and cellular integration of the nitroplast — analysis scripts

R scripts for the statistics and figures in the paper *Diel remodeling and cellular integration of the nitroplast*
(**TODO: authors, journal, DOI**).

The analyses compare nitroplast features (internal lamellae, globules, spiky vesicles, "hat" structures,
thylakoid widths, membrane contact sites) across four time points:
**Day, Evening, Night, Morning**.

## Repository layout

```
├── data/                      input CSVs (see data/README.md)
├── scripts/
│   ├── 00_config.R            shared settings: seed, paths, helper functions
│   ├── 01_statistical_tests_and_exploratory_plots.R
│   ├── 02_vesicle_hat_figures_4state.R
│   ├── 03_thylakoid_width.R
│   ├── 04_mcs_table.R
│   └── run_all.R              runs everything in order
├── archive/                   superseded code, kept for provenance (not run by run_all.R)
│   └── vesicle_hat_figures_superseded.R
│   ├── lamellae_figures.R
│   ├── combined_vesicle_hat_figures.R (overwritten by 02)
```

## What each script produces

| Script | Produces | Paper figure / table |
|---|---|---|
| `01_statistical_tests_and_exploratory_plots` | Day vs Night Welch t-tests and correlations (hats, lamellae, globules, vesicles); exploratory boxplots, violins, bar charts | TODO |
| `02_vesicle_hat_figures_4state` | **Four-time-point** Panels A–C: vesicle density, vesicle diameters, hat density; negative binomial model with Holm-adjusted pairwise contrasts | TODO |
| `03_thylakoid_width` | Thylakoid-width histograms by time point; logistic regression and t-test, Day vs Night | TODO |
| `04_mcs_table` | Membrane contact site percentages table | TODO |

## Running

Run everything **from the repository root** (the folder containing `data/` and `scripts/`):

Non-interactive runs write plots to `figures/<script>.pdf`; in RStudio, plots appear in the Plots pane
as usual (open the repo folder as the working directory).

## Dependencies

R ≥ 4.3 (developed/tested with R 4.3.3) and these packages:

```r
install.packages(c("dplyr", "ggplot2", "ggpubr", "patchwork", "MASS", "emmeans", "gt"))
```

`MASS` ships with R. `gt` is only needed for `06_mcs_table.R`. To lock exact versions for the paper, run
`renv::init()` then `renv::snapshot()` in the environment you used for the published figures.

## Reproducibility notes

- The only random element is `geom_jitter()` point placement (violin plots); the seed is set in `scripts/00_config.R`.
  All statistics are deterministic.
- Some inputs are hardcoded in the scripts rather than read from `data/`
  (Evening/Morning vesicle sizes, vesicle volumes, thylakoid widths, MCS percentages).
  See `data/README.md`.

## Data availability

TODO — data hosting is not yet decided (repo, Zenodo, institutional storage, or restricted). The three CSVs in
`data/` are included for convenience; remove them and link to the archive here if the data are hosted elsewhere.

## License and citation

TODO — add a `LICENSE` (MIT or BSD-3 are common for analysis code; CC-BY for data) and fill in `CITATION.cff`.
