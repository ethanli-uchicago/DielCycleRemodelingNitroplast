# ==============================================================================
# 06_mcs_table.R   (originally: MCS.R)
#
# Table of nitroplast membrane contact site (MCS) percentages, as % of nitroplast
# surface area, for each partner organelle.
#
# This is the FINAL version (cell 4 removed; cells renumbered 1-2). The earlier
# 3-cell table that preceded the "## cell 4 removed ##" marker was dropped.
#
# Inputs : none -- values are hardcoded from the author's measurements
# Output : gt table (RStudio Viewer).
# ==============================================================================
library(gt)
library(dplyr)

df <- data.frame(
  Cell = c("Cell 1", "Cell 2"),
  `Nitroplast-Mitochondria` = c("4.67 ± 0.06 %", "15.72 ± 0.02 %"),
  `Nitroplast-Chloroplast` = c("0.46 ± 0.05 %", "3.31 ± 0.06 %"),
  `Nitroplast-Lipid Vacuole` = c("None within threshold", "30.00 ± 0.03 %"),
  check.names = FALSE
)

df |>
  gt(rowname_col = "Cell") |>
  tab_header(
    title = "Nitroplast MCS Percentages (% of surface area)"
  ) |>
  cols_align(align = "center", columns = everything()) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_stub()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )
