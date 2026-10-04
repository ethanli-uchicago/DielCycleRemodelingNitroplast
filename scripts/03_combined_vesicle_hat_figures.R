# ==============================================================================
# 03_combined_vesicle_hat_figures.R   (originally: combinedplots.R)
#
# Vesicle totals, vesicle-diameter histograms and hat-structure totals, raw and
# normalized per nitroplast, combined side-by-side with patchwork.
# Time points here: Day / Night / Morning ("Morning" in the original; renamed).
#
# Inputs : data/Day_cleaned.csv, data/Night_cleaned.csv
# Output : figures/03_combined_vesicle_hat_figures.pdf
#
# Used in the paper (together with 04_vesicle_hat_figures_4state.R, which adds the
# Evening time point). The 2-3 time point versions of these plots in 01 and in
# archive/02_vesicle_hat_figures_superseded.R are superseded by this script.
# ==============================================================================
source("scripts/00_config.R")
library(dplyr)
library(ggplot2)
library(patchwork)
open_fig_device("03_combined_vesicle_hat_figures", width = 15, height = 5)

day   <- read_data("Day_cleaned.csv")
night <- read_data("Night_cleaned.csv")
hat_tbl <- read_data("hat_counts_by_tomogram.csv")   # per-tomogram hat counts: Day 1 / Evening 24 / Night 49 / Morning 31

day_ves_counts <- day$Spiky.Vesicles.Observed[!is.na(day$Spiky.Vesicles.Observed)]
night_ves_counts <- night$Spiky.Vesicles.Observed[!is.na(night$Spiky.Vesicles.Observed)]
morning_ves_counts <- c(0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

## SPIKY VESICLE AMOUNTS ##
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


# Normalized (per nitroplast count)
avg_day <- sum(day_ves_counts) / N_NITRO_DAY
avg_night <- sum(night_ves_counts) / N_NITRO_NIGHT
avg_morning <- 3 / N_NITRO_MORNING
barchart_data <- data.frame(
  Period = factor(c("Day", "Night", "Morning"), 
                  levels = c("Day", "Night", "Morning")),
  NormCounts = c(avg_day, avg_night, avg_morning)
)
p1_2 <- ggplot(barchart_data, aes(x = Period, y = NormCounts, fill = Period)) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  scale_y_continuous(
    limits = c(0, 12),           # Sets the total range
    breaks = seq(0, 12, by = 2), # Optional: Sets specific tick intervals (e.g., every 2 units)
    expand = c(0, 0)             # Prevents the bars from "floating" above the x-axis
  ) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue", "Morning" = "lightblue")) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    #axis.line.x = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## SPIKY VESICLE DIAMETERS ##
day_ves_size <- day$Vesicle.Diameter..widest.[!is.na(day$Vesicle.Diameter..widest.)]
night_ves_size <- night$Vesicle.Diameter..widest.[!is.na(night$Vesicle.Diameter..widest.)]
morning_ves_size <- c(54,50,50)
all_ves_size <- c(day_ves_size, night_ves_size, morning_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_morning_labels <- rep("Morning", 3)
ves_time_labels <- c(ves_day_labels, ves_night_labels, ves_morning_labels)
plot_data <- data.frame(
  Counts = all_ves_size,
  Period = factor(ves_time_labels, levels = c("Day","Night","Morning"))
)
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

# Normalized (per nitroplast count)
# 1. Define your nitroplast counts (Replace these with your actual data)
num_nitro_day <- N_NITRO_DAY
num_nitro_night <- N_NITRO_NIGHT
num_nitro_morning <- N_NITRO_MORNING

# 2. Create a lookup vector or add weights to your dataframe
# This calculates (1 / N) for each row, so the sum of counts = Frequency per Nitroplast
weights_map <- c("Day" = 1/num_nitro_day, "Night" = 1/num_nitro_night, "Morning" = 1/num_nitro_morning)
plot_data$weight <- weights_map[as.character(plot_data$Period)]

# 3. Plot using the weight aesthetic
p2_2 <- ggplot(plot_data, aes(x = Counts, fill = Period, weight = weight)) +
  # Using weight here makes the y-axis (sum of weights) equal to counts/N
  geom_histogram(binwidth = 10, alpha = 0.7, color="black", linetype="dashed", position = "identity") +
  
  # For density, we scale it by the binwidth (10) to match the histogram height
  geom_density(aes(y = after_stat(density) * 10), alpha = 0.4) +
  
  theme_classic() +
  labs(x = "", y = "") +
  #labs(
  #  title = "Distribution of Vesicle Diameters per Nitroplast",
  #  x = "Size of Vesicles Diameters (nm)",
  #  y = "Frequency per Nitroplast"
  #) +
  theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
  #scale_fill_manual(values = c("Night" = "slateblue", "Day" = "orange", "Morning" = "lightblue")) +
  scale_y_continuous(
    expand = c(0, 0)             # Prevents the bars from "floating" above the x-axis
  ) +
  scale_fill_manual(values = c("Night" = "slateblue", "Day" = "orange", "Morning" = "lightblue"),labels = c("", "", "")) +
  coord_cartesian(xlim = c(NA, 140))

## HAT LIKE STRUCTURE AMOUNTS ##
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

# Normalized
# 1. Define your nitroplast counts (ensure these match your previous plot)
num_nitro_day <- N_NITRO_DAY
num_nitro_night <- N_NITRO_NIGHT
num_nitro_morning <- N_NITRO_MORNING

# 2. Add the nitroplast counts to your dataframe
barchart_data$NitroCount <- c(num_nitro_night, num_nitro_morning, num_nitro_day)

# 3. Calculate the normalized frequency
barchart_data$HatsPerNitro <- barchart_data$TotalCounts / barchart_data$NitroCount

# 4. Plot
p3_2 <- ggplot(barchart_data, aes(x = Period, y = HatsPerNitro, fill = Period)) +
  geom_bar(
    stat = "identity",
    width = 0.4,
    color = "black"
  ) +
  theme_classic() +
  # We round the label for readability
  geom_text(aes(label = round(HatsPerNitro, 2)), vjust = -0.5, size = 5) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()) +
  scale_fill_manual(values = c("Night" = "slateblue", "Morning" = "lightblue", "Day" = "orange")) +
  # Ensure the y-axis starts exactly at 0
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  guides(fill = "none")



## --- Your first plot ---
total_day <- sum(day_ves_counts)
total_night <- sum(night_ves_counts)
total_morning <- 3
barchart_data <- data.frame(
  Period = factor(c("Day", "Night", "Morning"), 
                  levels = c("Day", "Night", "Morning")),
  TotalCounts = c(total_day, total_night, total_morning)
)
p1 <- ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(
    aes(y = TotalCounts + 0.3),
    stat = "identity",
    width = 0.4,
    color = "black"
  ) +
  labs(
    title = "Total Vesicles Observed",
    x = NULL,
    y = "Number of Vesicles"
  ) +
  scale_fill_manual(values = c("Night" = "slateblue", "Morning" = "lightblue", "Day" = "orange")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  guides(fill = "none") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),   # remove built-in x-axis
    axis.ticks.x = element_blank()
  ) +
  geom_hline(yintercept = 0, color = "black")   # standardized baseline
p1


## --- Your second plot ---
p2 <- ggplot(plot_data, aes(x = Counts, fill = Period)) +
  geom_histogram(binwidth = 10, alpha = 0.7, color="black", linetype="dashed") +
  geom_density(aes(y=after_stat(count)*10), alpha = 0.4) +
  labs(
    title = "Distribution of Vesicle Diameters",
    x = "Vesicles diameter (nm)",
    y = "Number of vesicles"
  ) +
  coord_cartesian(xlim = c(NA, 140), ylim = c(-1.5,NA)) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) +
  scale_fill_manual(
    name = "Sample Type",
    values = c("Night" = "slateblue", "Day" = "orange", "Morning" = "lightblue"),
    labels = c(
      "Day" = paste0("Day, ", N_NITRO_DAY, " tomograms"), 
      "Night" = paste0("Night, ", N_NITRO_NIGHT, " tomograms"), 
      "Morning" = paste0("Morning, ", N_NITRO_MORNING, " tomograms")
    )
  ) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "black") + 
  geom_segment(data = data.frame(x = seq(0, 120, by = 40)), 
               aes(x = x, xend = x, y = 0, yend = -1), inherit.aes = FALSE) +
  annotate("text", x = seq(0, 120, by = 40), y = -5, 
           label = seq(0, 120, by = 40), vjust = -2) +
  theme(
    # Center the plot title
    plot.title = element_text(hjust = 0.5),
    
    # --- Legend Customization ---
    legend.position = "bottom",          # Position the legend at the top
    legend.direction = "horizontal",  # Arrange legend items horizontally
    
    # --- Custom X-axis ---
    axis.line.x = element_blank(),    # Remove built-in x-axis line
    axis.text.x = element_blank(),    # Remove built-in x-axis text
    axis.ticks.x = element_blank(),    # Remove built-in x-axis ticks
    
    axis.title.x = element_text(margin = margin(t = -10, unit = "pt"))
  )

# Print the plot
p2



## --- Your third plot ---
day_hat_number   <- hat_tbl$hats[hat_tbl$state == "Day"]     # 8 tomograms, 1 hat
night_hat_number <- hat_tbl$hats[hat_tbl$state == "Night"]   # 27 tomograms, 49 hats
total_day <- sum(day_hat_number)
total_night <- sum(night_hat_number)
total_morning <- sum(hat_tbl$hats[hat_tbl$state == "Morning"])   # 31 (was hardcoded 25)
barchart_data2 <- data.frame(
  Period = factor(c("Night", "Morning", "Day"), 
                  levels = c("Day", "Night", "Morning")),
  TotalCounts = c(total_night, total_morning, total_day)
)
p3 <- ggplot(barchart_data2, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(
    aes(y = TotalCounts + 0.3),
    stat = "identity",
    width = 0.4,
    color = "black"
  ) +
  labs(
    title = "Hat-like structures observed",
    x = NULL,
    y = "Number of structures"
  ) +
  scale_fill_manual(values = c("Night" = "slateblue", "Morning" = "lightblue", "Day" = "orange")) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),   # remove built-in x-axis
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) +
  geom_hline(yintercept = 0, color = "black")   # standardized baseline

## --- Combine them side by side ---
p1 | p2 | p3

p1_2 | p2_2 | p3_2
