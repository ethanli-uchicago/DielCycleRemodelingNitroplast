# ==============================================================================
# ARCHIVED -- superseded; not used in the paper.
# 02_vesicle_hat_figures_superseded.R   (originally: the vesicle/hat part of featurevisualization.R)
#
# Earlier 2-3 time point (Day / Night / Morning) versions of the vesicle and
# hat-structure figures: total vesicles, vesicles per nitroplast (boxplot),
# vesicle diameter (violin and histogram) and hat totals, each labeled + BLANK.
# The paper's versions come from scripts/03_combined_vesicle_hat_figures.R and
# scripts/04_vesicle_hat_figures_4state.R.
#
# Kept for provenance; still runs from the repository root:
#     Rscript archive/02_vesicle_hat_figures_superseded.R
# (Note: the per-nitroplast boxplot and vesicle-diameter violin exist ONLY here
#  and in 01; they have no counterpart in 03/04.)
# ==============================================================================
source("scripts/00_config.R")
library(dplyr)
library(ggplot2)
open_fig_device("archive_02_vesicle_hat_figures_superseded")

day   <- read_data("Day_cleaned.csv")
night <- read_data("Night_cleaned.csv")
hat_tbl <- read_data("hat_counts_by_tomogram.csv")   # per-tomogram hat counts: Day 1 / Evening 24 / Night 49 / Morning 31

# ---- (RAW) VESICLE AMOUNT ----
day_ves_counts <- day$Spiky.Vesicles.Observed[!is.na(day$Spiky.Vesicles.Observed)]
night_ves_counts <- night$Spiky.Vesicles.Observed[!is.na(night$Spiky.Vesicles.Observed)]
morning_ves_counts <- c(0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

total_day <- sum(day_ves_counts)
total_night <- sum(night_ves_counts)
total_morning <- 3 #I counted as of 7.25.25. There's only three.

#barchart_data <- data.frame(
#  Period = c("Day", "Night", "Morning"),
#  TotalCounts = c(total_day, total_night, total_morning)
#)

barchart_data <- data.frame(
  Period = c("Night", "Morning", "Day"),
  TotalCounts = c(total_night, total_morning, total_day)
)
barchart_data$Period <- factor(barchart_data$Period, levels = c("Night", "Morning", "Day"))

ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  geom_text(aes(label = TotalCounts), vjust = 1.5, size = 3, color = "white") +
  labs(
    title = "Total Number of Vesicles Observed",
    #subtitle = "Sum of all vesicles found during the day versus the night",
    x = "Metabolic State",
    y = "Total Count of Vesicles"
  ) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue", "Morning" = "lightblue")) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- VESICLES PER NITROPLAST ----
#all_ves_counts <- c(day_ves_counts, night_ves_counts, morning_ves_counts)
all_ves_counts <- c(night_ves_counts, morning_ves_counts, day_ves_counts)
ves_day_labels <- rep("Day", length(day_ves_counts))
ves_night_labels <- rep("Night", length(night_ves_counts))
ves_trans_labels <- rep("Morning", length(morning_ves_counts))
#ves_time_labels <- c(ves_day_labels, ves_night_labels, ves_trans_labels)
ves_time_labels <- c(ves_night_labels, ves_trans_labels, ves_day_labels)
plot_data <- data.frame(
  Counts = all_ves_counts,
  Period = factor(ves_time_labels, levels = c("Night","Morning","Day")) # factors should be in the correct order
)
#plot_data$Counts <- factor(plot_data$Counts, levels = c("Night", "Morning", "Day"))


ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  theme_classic() +
  labs(
    title = "Comparison of Vesicle Amounts: Day vs. Night",
    subtitle = "Each box shows the distribution of vesicle counts per nitroplast",
    x = "Metabolic States",
    y = "Number of Vesicles Observed"
  ) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue", "Morning" = "lightblue")) +
  annotate(
    geom = "text",
    x = c("Day", "Night", "Morning"),
    y = -Inf,       # Sets the y-position to the very bottom
    vjust = -1.5,   # Adjusts the labels to sit just below the plot
    label = paste0("n=", c(length(day_ves_counts), length(night_ves_counts), length(morning_ves_counts))),   # was hardcoded c("n=6", "n=17", "n=20")
    size = 4.5
  )

# ---- VESICLE SIZE VIOLIN PLOT ----
day_ves_size <- day$Vesicle.Diameter..widest.[!is.na(day$Vesicle.Diameter..widest.)]
night_ves_size <- night$Vesicle.Diameter..widest.[!is.na(night$Vesicle.Diameter..widest.)]
morning_ves_size <- c(54,50,50)
all_ves_size <- c(day_ves_size, night_ves_size, morning_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_morning_labels <- rep("Morning", 3)
ves_time_labels <- c(ves_day_labels, ves_night_labels, ves_morning_labels)
# Violin for Vesicle Size: Day vs. Night
#all_ves_size <- c(night_ves_size, day_ves_size)
#ves_time_labels <- c(ves_night_labels, ves_day_labels)
#plot_data <- data.frame(
#  Counts = all_ves_size,
#  Period = factor(ves_time_labels, levels = c("Day","Night"))
#)
plot_data <- data.frame(
  Counts = all_ves_size,
  Period = factor(ves_time_labels, levels = c("Day","Night","Morning")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  theme_classic() +
  geom_jitter(width = 0.1, alpha = 0.4, color = "black") +
  labs(
    title = "Distribution of Vesicle Diameters: Day vs. Night",
    subtitle = "(Dots are individual nitroplasts)",
    x = "Metabolic State",
    y = "Size of Vesicles Diameters (nm)"
  ) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue", "Morning" = "lightblue")) +
  guides(fill = "none")

# ---- VESICLE SIZE DENSITY PLOT ----
# Assuming 'plot_data' is already loaded
ggplot(plot_data, aes(x = Counts, fill = Period)) +
  geom_histogram(binwidth = 10, alpha = 0.7, color="black", linetype="dashed") +
  #geom_histogram(binwidth = 10, alpha = 0.7) +
  geom_density(aes(y=after_stat(count)*10),alpha = 0.4) +
  #geom_histogram() +
  theme_classic() +
  labs(
    title = "Distribution of Vesicle Diameters",
    x = "Size of Vesicles Diameters (nm)",
    y = "Count"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Night" = "slateblue", "Day" = "orange", "Morning" = "lightblue")) +
  coord_cartesian(xlim = c(NA, 140))

# ---- HAT STRUCTURE BAR GRAPH ----
day_hat_number   <- hat_tbl$hats[hat_tbl$state == "Day"]     # 8 tomograms, 1 hat
night_hat_number <- hat_tbl$hats[hat_tbl$state == "Night"]   # 27 tomograms, 49 hats
total_day <- sum(day_hat_number)
total_night <- sum(night_hat_number)
total_morning <- sum(hat_tbl$hats[hat_tbl$state == "Morning"])   # 31 (was hardcoded 25)
barchart_data <- data.frame(
  Period = factor(c("Night", "Morning", "Day"), 
                  levels = c("Day", "Night", "Morning")),
  TotalCounts = c(total_night, total_morning, total_day)
)

ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(
    aes(y = TotalCounts + 0.3),
    stat = "identity",
    width = 0.4,
    color = "black"
  ) +
  theme_classic() +
  geom_text(aes(label = TotalCounts), y = 1.2, size = 5) +
  labs(
    title = "Total Hat-Structures Observed",
    x = "Metabolic State",
    y = "Total Hat-Structures Observed"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Night" = "slateblue", "Morning" = "lightblue", "Day" = "orange")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none")


# ---- BLANKS (no titles/labels, for figure layout) ----
# RAW VESICLE AMOUNT
total_day <- sum(day_ves_counts)
total_night <- sum(night_ves_counts)
total_morning <- 3
barchart_data <- data.frame(
  Period = factor(c("Day", "Night", "Morning"), 
                  levels = c("Day", "Night", "Morning")),
  TotalCounts = c(total_day, total_night, total_morning)
)
ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue", "Morning" = "lightblue")) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# VESICLES PER NITROPLAST
all_ves_counts <- c(day_ves_counts, night_ves_counts)
ves_day_labels <- rep("Day", length(day_ves_counts))
ves_night_labels <- rep("Night", length(night_ves_counts))
ves_time_labels <- c(ves_day_labels, ves_night_labels)
plot_data <- data.frame(
  Counts = all_ves_counts,
  Period = factor(ves_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5), 
    axis.text.x = element_blank(),
    legend.text = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL)

# VESICLE SIZE VIOLINPLOT
day_ves_size <- day$Vesicle.Diameter..widest.[!is.na(day$Vesicle.Diameter..widest.)]
night_ves_size <- night$Vesicle.Diameter..widest.[!is.na(night$Vesicle.Diameter..widest.)]
all_ves_size <- c(day_ves_size, night_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_time_labels <- c(ves_day_labels, ves_night_labels)
plot_data <- data.frame(
  Counts = all_ves_size,
  Period = factor(ves_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  theme_classic() +
  geom_jitter(width = 0.1, alpha = 0.4, color = "black") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none") + 
  labs(title = NULL, subtitle = NULL, x = NULL, y = "Vesicle Diameter (nm)")

# VESICLE SIZE DENSITY CURVE
ggplot(plot_data, aes(x = Counts, fill = Period)) +
  #geom_histogram(binwidth = 10, alpha = 0.7, color="black", linetype="dashed") +
  geom_histogram(binwidth = 10, alpha = 0.7) +
  geom_density(aes(y=after_stat(count)*10),alpha = 0.4) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_blank(),
        legend.title = element_blank()
        ) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  coord_cartesian(xlim = c(NA, 140))


# HAT STRUCTURE BAR GRAPH
day_hat_number   <- hat_tbl$hats[hat_tbl$state == "Day"]     # 8 tomograms, 1 hat
night_hat_number <- hat_tbl$hats[hat_tbl$state == "Night"]   # 27 tomograms, 49 hats
total_day <- sum(day_hat_number)
total_night <- sum(night_hat_number)
barchart_data <- data.frame(
  Period = c("Day", "Night"),
  TotalCounts = c(total_day, total_night)
)
ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(
    aes(y = TotalCounts + 0.3),
    stat = "identity",
    width = 0.4,
    color = "black"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    ) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL)


