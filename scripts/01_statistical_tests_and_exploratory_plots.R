# ==============================================================================
# 01_statistical_tests_and_exploratory_plots.R   (originally: analysis.R)
#
# Day vs Night comparisons of hat structures, internal lamellae type, globules
# and spiky vesicles (Welch t-tests and Pearson correlations), followed by
# exploratory plots (boxplots, violins, bar charts).
#
# Inputs : data/Day_cleaned.csv, data/Night_cleaned.csv
# Outputs: console statistics (see results/ when run via run_all.R) and
#          figures/01_statistical_tests_and_exploratory_plots.pdf
#
# NOTE: only two time points (Day, Night) appear here. The paper's vesicle and hat
# figures come from 03 and 04 (four time points), so the vesicle/hat plots in the
# "Visualizations" section at the end of this script are exploratory/superseded.
#
# (The globule plots exist only here)
# Recorded p-values in comments are the author's; where they no longer match
# the current data they are marked NOTE(review). See REVIEW_NOTES.md.
# ==============================================================================
source("scripts/00_config.R")
library(dplyr)
library(ggplot2)
open_fig_device("01_statistical_tests_and_exploratory_plots")

# --- Load data ---
day   <- read_data("Day_cleaned.csv")
night <- read_data("Night_cleaned.csv")
hat_tbl <- read_data("hat_counts_by_tomogram.csv")   # per-tomogram hat counts: Day 1 / Evening 24 / Night 49 / Morning 31


# --- Hat Analysis ---
# --- (measure is at least one hat structure observed) ---

#day_hat_counts <- day$Hat.Structure.Observed[!is.na(day$Hat.Structure.Observed)]
#night_hat_counts <- night$Hat.Structures.Observed.[!is.na(night$Hat.Structures.Observed.)]
day_hat_counts <- numeric()
for(i in 1:length(day$Hat.Structure.Observed)){
  #print(i)
  if(day$Hat.Structure.Observed[i] == "No"){
    day_hat_counts[i] <- 0
  }
  else if(day$Hat.Structure.Observed[i] == "Yes"){
    day_hat_counts[i] <- 1
  }
  else if(day$Hat.Structure.Observed[i] == "Maybe"){
    day_hat_counts[i] <- 0.5
  }
  else{
    day_hat_counts[i] <- NA
  }
}
day_hat_counts <- day_hat_counts[!is.na(day_hat_counts)]
#vector showing if there was a hat or not. 1 is yes, 0 is now, 0.5 is uncertain.
#print(day_hat_counts)

night_hat_counts <- numeric()
for(i in 1:length(night$Hat.Structures.Observed.)){
  #print(i)
  if(night$Hat.Structures.Observed.[i] == "No"){
    night_hat_counts[i] <- 0
  }
  else if(night$Hat.Structures.Observed.[i] == "Yes"){
    night_hat_counts[i] <- 1
  }
  else if(night$Hat.Structures.Observed.[i] == "Maybe"){
    night_hat_counts[i] <- 0.5
  }
  else{
    night_hat_counts[i] <- NA
  }
}
night_hat_counts <- night_hat_counts[!is.na(night_hat_counts)]
#print(night_hat_counts)

t.test(day_hat_counts, night_hat_counts)
# t = -3.661, df = 14.803, p-value = 0.002363. Statistically significant
# NOTE(review): recorded value is from an earlier data version (32 coded tomograms; now 30).
#   Current data (Day tomogram 1 = "Yes") give t = -3.1231, df = 11.499, p = 0.009227.

all_hat_counts <- c(day_hat_counts, night_hat_counts)
hat_day_labels <- rep("Day", length(day_hat_counts))
hat_night_labels <- rep("Night", length(night_hat_counts))
hat_time_labels <- c(hat_day_labels, hat_night_labels)
hat_combined_data <- data.frame(
  Counts = all_hat_counts,
  Time = hat_time_labels
)
# Pearson's correlation
hat_correlation_test_amount <- cor.test(hat_combined_data$Counts, as.numeric(factor(hat_combined_data$Time)))
print(hat_correlation_test_amount)
# current data gives t = 2.8654, df = 28, p = 0.007815. Statistically significant.


# --- Hat Number Analysis ---
# --- (Of the observable ones I found, there will be some error with this) ---
# Hat NUMBER per tomogram comes from data/hat_counts_by_tomogram.csv (the same table 04 uses):
# Day = 8 tomograms / 1 hat, Night = 27 tomograms / 49 hats. (Night_cleaned.csv's NumHat column is
# incomplete -- 23 tomograms / 38 hats -- so it is no longer used for this test; see REVIEW_NOTES.md.)
day_hat_number   <- hat_tbl$hats[hat_tbl$state == "Day"]
night_hat_number <- hat_tbl$hats[hat_tbl$state == "Night"]

t.test(day_hat_number, night_hat_number)
# NOTE(review): now computed from data/hat_counts_by_tomogram.csv (Day 8 tomograms / 1 hat; Night 27 / 49):
#   t = -4.9131, df = 31.782, p = 2.599e-05.

  
# --- Lamellae Type Analysis ---
# --- (I classified as thick, thin, thick/thin (or otherway around), and none) ---
day_lamellae <- numeric()
for(i in 1:length(day$Internal.Lamellae.Type)){
  #print(i)
  if(day$Internal.Lamellae.Type[i] == "Thin"){
    day_lamellae[i] <- 0
  }
  else if(day$Internal.Lamellae.Type[i] == "Thick"){
    day_lamellae[i] <- 1
  }
  else if(day$Internal.Lamellae.Type[i] == "None"){
    day_lamellae[i] <- -1
  }
  else if(grepl("and", day$Internal.Lamellae.Type[i])){
    day_lamellae[i] <- 0.5
  }
  else{
    day_lamellae[i] <- NA
  }
}
day_lamellae <- day_lamellae[!is.na(day_lamellae)] 
#vector of each nitroplast lamellae type: 1 is thick, 0.5 is containing both, 0 is thin, -1 is not observed

night_lamellae <- numeric()
for(i in 1:length(night$Internal.Lamellae.Type)){
  #print(i)
  if(night$Internal.Lamellae.Type[i] == "Thin"){
    night_lamellae[i] <- 0
  }
  else if(night$Internal.Lamellae.Type[i] == "Thick"){
    night_lamellae[i] <- 1
  }
  else if(night$Internal.Lamellae.Type[i] == "None"){
    night_lamellae[i] <- -1
  }
  else if(grepl("and", night$Internal.Lamellae.Type[i])){ ### FIX IS HERE ###
    night_lamellae[i] <- 0.5
  }
  else{
    night_lamellae[i] <- NA
  }
}
night_lamellae <- night_lamellae[!is.na(night_lamellae)]
#vector of each nitroplast lamellae type: 1 is thick, 0.5 is containing both, 0 is thin, -1 is not observed

t.test(day_lamellae,night_lamellae)
# t = -0.25324, df = 13.053, p-value = 0.804. Not statistically significant.

all_lamellae <- c(day_lamellae, night_lamellae)
lamellae_day_labels <- rep("Day", length(day_lamellae))
lamellae_night_labels <- rep("Night", length(night_lamellae))
lamellae_time_labels <- c(lamellae_day_labels, lamellae_night_labels)
lamellae_combined_data <- data.frame(
  Type = all_lamellae,
  Time = lamellae_time_labels
)
lamellae_correlation_test_size <- cor.test(lamellae_combined_data$Type, as.numeric(factor(lamellae_combined_data$Time)))
print(lamellae_correlation_test_size)
# t = 0.21797, df = 28, p-value = 0.829. Not statistically significant

lamellae_hat_ttest <- t.test(lamellae_combined_data$Type, hat_combined_data$Counts)
# NOTE(review): current data give t = -2.4039, df = 57.987, p = 0.01944. Statistically significant
# Exploratory and not used in the paper
















# --- Globule Analysis ---
# 1. AMOUNT: Day vs Night
day_globule_counts <- day$Globules.Observed[!is.na(day$Globules.Observed)]
night_globule_counts <- night$Globules.Observed[!is.na(night$Globules.Observed)]
#print(day_globule_counts) # vector showing amount of globules per vector per nitroplast
#print(night_globule_counts)

# Welch two-sample t-test,
t.test(day_globule_counts, night_globule_counts) 
# this test will tell if there is a statistically significant difference between the average number of globules per nitroplast during the day versus during the night
# t = 1.1186, df = 14.607, p-value = 0.2814. Not Statistically significant

# 2. SIZE: Day vs Night
day_globule_sizes <- day$Globule.Diameters..Widest.[!is.na(day$Globule.Diameters..Widest.)]
night_globule_sizes <- night$Globule.Diameters..Widest.[!is.na(night$Globule.Diameters..Widest.)]
#print day_globule_sizes
#print night_globule_sizes

t.test(day_globule_sizes,  night_globule_sizes)
# t = 1.4339, df = 33.53, p-value = 0.1609. Not statistically significant

# 3. CORRELATION
# 3.1 AMOUNT
all_globule_counts <- c(day_globule_counts, night_globule_counts)
glob_day_labels <- rep("Day", length(day_globule_counts))
glob_night_labels <- rep("Night", length(night_globule_counts))
glob_time_labels <- c(glob_day_labels, glob_night_labels)
glob_combined_data <- data.frame(
  Counts = all_globule_counts,
  Time = glob_time_labels
)
# Pearson's correlation
glob_correlation_test_amount <- cor.test(glob_combined_data$Counts, as.numeric(factor(glob_combined_data$Time)))
print(glob_correlation_test_amount)
# t = -0.90258, df = 29, p-value = 0.3742. Not statistically significant

# 3.2 SIZE
all_globule_size <- c(day_globule_sizes, night_globule_sizes)
glob_day_labels <- rep("Day", length(day_globule_sizes))
glob_night_labels <- rep("Night", length(night_globule_sizes))
glob_time_labels <- c(glob_day_labels, glob_night_labels)
glob_combined_data <- data.frame(
  Sizes = all_globule_size,
  Time = glob_time_labels
)
glob_correlation_test_size <- cor.test(glob_combined_data$Sizes, as.numeric(factor(glob_combined_data$Time)))
print(glob_correlation_test_size)
# t = -1.6334, df = 69, p-value = 0.1069. Not statistically significant



# --- Vesicle (spiky) Analysis ---
# 1. AMOUNT: Day vs Night
day_ves_counts <- day$Spiky.Vesicles.Observed[!is.na(day$Spiky.Vesicles.Observed)]
night_ves_counts <- night$Spiky.Vesicles.Observed[!is.na(night$Spiky.Vesicles.Observed)]
#print(day_globule_counts) # vector showing amount of vesicles per vector per nitroplast
#print(night_globule_counts)

t.test(day_ves_counts, night_ves_counts)
# t = 2.8613, df = 5.0198, p-value = 0.03519. Statistically significant

#2. SIZE: Day vs Night
day_ves_size <- day$Vesicle.Diameter..widest.[!is.na(day$Vesicle.Diameter..widest.)]
night_ves_size <- night$Vesicle.Diameter..widest.[!is.na(night$Vesicle.Diameter..widest.)]

day_ves_Q1 <- quantile(day_ves_size, 0.25)
day_ves_Q3 <- quantile(day_ves_size, 0.75)
day_ves_IQR_val <- IQR(day_ves_size)
day_ves_lower_bound <- day_ves_Q1 - 1.5 * day_ves_IQR_val
day_ves_upper_bound <- day_ves_Q3 + 1.5 * day_ves_IQR_val
day_ves_size_no_outliers <- day_ves_size[day_ves_size >= day_ves_lower_bound & day_ves_size <= day_ves_upper_bound]
print(day_ves_size_no_outliers)
print(mean(day_ves_size_no_outliers))

t.test(day_ves_size,night_ves_size)
# t = 4.7155, df = 43.432, p-value = 2.51e-05
t.test(day_ves_size_no_outliers,night_ves_size)
# t = 4.6605, df = 10.664, p-value = 0.0007529

# 3. CORRELATION
# 3.1 AMOUNT
all_ves_counts <- c(day_ves_counts, night_ves_counts)
ves_day_labels <- rep("Day", length(day_ves_counts))
ves_night_labels <- rep("Night", length(night_ves_counts))
ves_time_labels <- c(ves_day_labels, ves_night_labels)
ves_combined_data <- data.frame(
  Counts = all_ves_counts,
  Time = ves_time_labels
)
# Pearson's correlation
ves_correlation_test_amount <- cor.test(ves_combined_data$Counts, as.numeric(factor(ves_combined_data$Time)))
print(ves_correlation_test_amount)
# t = --5.9549, df = 28, p-value = 2.067e-06. Statistically significant

# 3.2 SIZE
all_ves_size <- c(day_ves_size, night_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_time_labels <- c(ves_day_labels, ves_night_labels)
ves_combined_data <- data.frame(
  Sizes = all_ves_size,
  Time = ves_time_labels
)
ves_correlation_test_size <- cor.test(ves_combined_data$Sizes, as.numeric(factor(ves_combined_data$Time)))
print(ves_correlation_test_size)
# t = -1.762, df = 94, p-value = 0.08132. Not statistically significant, but close



# --- Visualizations ---

glob_day_labels <- rep("Day", length(day_globule_counts))
glob_night_labels <- rep("Night", length(night_globule_counts))
glob_time_labels <- c(glob_day_labels, glob_night_labels)

# Boxplot for Globule Amount: Day vs. Night
plot_data <- data.frame(
  Counts = all_globule_counts,
  Period = factor(glob_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Comparison of Globule Amounts: Day vs. Night",
    subtitle = "Each box shows the distribution of globule counts per nitroplast",
    x = "Metabolic States",
    y = "Number of Globules Observed"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue"))

# Violinplot for Globule Amount
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  geom_jitter(width = 0.1, alpha = 0.6, color = "black") +
  labs(
    title = "Distribution of Globule Amounts: Day vs. Night",
    subtitle = "(Dots are individual nitroplasts)",
    x = "Metabolic State",
    y = "Number of Globules Observed"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")

# Bar Chart for raw globule amount
total_day <- sum(day_globule_counts)
total_night <- sum(night_globule_counts)
barchart_data <- data.frame(
  Period = c("Day", "Night"),
  TotalCounts = c(total_day, total_night)
)
ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(stat = "identity", position = "dodge") + # stat="identity" plots the raw values
  geom_text(aes(label = TotalCounts), vjust = -0.5, size = 5) + # Add count labels above bars
  labs(
    title = "Total Number of Globules Observed",
    subtitle = "Comparing the sum of all globules found during the day versus the night",
    x = "Metabolic State",
    y = "Total Raw Count of Globules"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")

# --- globule size ---
glob_day_labels <- rep("Day", length(day_globule_sizes))
glob_night_labels <- rep("Night", length(night_globule_sizes))
glob_time_labels <- c(glob_day_labels, glob_night_labels)
# Violin for Globule Size: Day vs. Night
plot_data <- data.frame(
  Counts = all_globule_size,
  Period = factor(glob_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  geom_jitter(width = 0.1, alpha = 0.6, color = "black") +
  labs(
    title = "Distribution of Globule Sizes: Day vs. Night",
    subtitle = "(Dots are individual nitroplasts)",
    x = "Metabolic State",
    y = "Size of Globules Observed (nm)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")







all_ves_counts <- c(day_ves_counts, night_ves_counts)
ves_day_labels <- rep("Day", length(day_ves_counts))
ves_night_labels <- rep("Night", length(night_ves_counts))
ves_time_labels <- c(ves_day_labels, ves_night_labels)


# Boxplot for Vesicle Amount: Day vs. Night
plot_data <- data.frame(
  Counts = all_ves_counts,
  Period = factor(ves_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Comparison of Vesicle Amounts: Day vs. Night",
    subtitle = "Each box shows the distribution of vesicle counts per nitroplast",
    x = "Metabolic States",
    y = "Number of Vesicles Observed"
  ) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  annotate(
    geom = "text",
    x = c("Day", "Night"),
    y = 23,
    vjust = 2,
    label = paste0("n=", c(length(day_ves_counts), length(night_ves_counts))),   # was hardcoded c("n=6", "n=17")
    size = 4.5
  )

# Violinplot for Vesicle Amount
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  geom_jitter(width = 0.1, alpha = 0.6, color = "black") +
  labs(
    title = "Distribution of Vesicle Amounts: Day vs. Night",
    subtitle = "(Dots are individual nitroplasts)",
    x = "Metabolic State",
    y = "Number of Vesicles Observed"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")

# Bar Chart for raw vesicle amount
total_day <- sum(day_ves_counts)
total_night <- sum(night_ves_counts)
barchart_data <- data.frame(
  Period = c("Day", "Night"),
  TotalCounts = c(total_day, total_night)
)
ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(stat = "identity", position = "dodge") + # stat="identity" plots the raw values
  geom_text(aes(label = TotalCounts), vjust = -0.5, size = 5) + # Add count labels above bars
  labs(
    title = "Total Number of Vesicles Observed",
    subtitle = "Comparing the sum of all vesicles found during the day versus the night",
    x = "Metabolic State",
    y = "Total Raw Count of Vesicles"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")




# --- vesicle size ---
all_ves_size <- c(day_ves_size, night_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_time_labels <- c(ves_day_labels, ves_night_labels)
# Violin for Vesicle Size: Day vs. Night
plot_data <- data.frame(
  Counts = all_ves_size,
  Period = factor(ves_time_labels, levels = c("Day", "Night")) # factors should be in the correct order
)
ggplot(plot_data, aes(x = Period, y = Counts, fill = Period)) +
  geom_violin() + # Creates the violin shape
  geom_jitter(width = 0.1, alpha = 0.6, color = "black") +
  labs(
    title = "Distribution of Vesicle Diameters: Day vs. Night",
    subtitle = "(Dots are individual nitroplasts)",
    x = "Metabolic State",
    y = "Size of Vesicles Diameters (nm)"
  ) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")



# --- bar chart for hats observed ---
total_day <- sum(day_hat_number)
total_night <- sum(night_hat_number)
barchart_data <- data.frame(
  Period = c("Day", "Night"),
  TotalCounts = c(total_day, total_night)
)
ggplot(barchart_data, aes(x = Period, y = TotalCounts, fill = Period)) +
  geom_bar(stat = "identity", position = "dodge") + # stat="identity" plots the raw values
  geom_text(aes(label = TotalCounts), vjust = -0.5, size = 5) + # Add count labels above bars
  labs(
    title = "Total Number of Hat-Structures Observed",
    #subtitle = "Comparing the sum of all vesicles found during the day versus the night",
    x = "Metabolic State",
    y = "Total Hat-Structures Observed"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Day" = "orange", "Night" = "slateblue")) +
  guides(fill = "none")
