# ==============================================================================
# 04_vesicle_hat_figures_4state.R   (originally: Vesicle_HatPlot.R)
#
# Used in the paper. FINAL four-time-point figures: Day / Evening / Night / Morning.
#   Panel A: vesicles per um^3 of matrix
#   Panel B: vesicle diameter histograms + density
#   Panel C: hat structures per um^3 of matrix, plus a negative binomial model
#            (hats ~ state + offset(log(volume))) with Holm-adjusted pairwise tests
# Each panel is drawn twice: labeled, then BLANK (for figure layout).
#
# Inputs : data/Day_cleaned.csv, data/Night_cleaned.csv
#          (Evening/Morning vesicle sizes, per-period counts/volumes and
#           per-tomogram hat counts are hardcoded below -- see data/README.md)
# Output : figures/04_vesicle_hat_figures_4state.pdf
# ==============================================================================
source("scripts/00_config.R")
library(ggplot2)
library(MASS)       # glm.nb()
library(emmeans)    # pairwise contrasts
open_fig_device("04_vesicle_hat_figures_4state")

day   <- read_data("Day_cleaned.csv")
night <- read_data("Night_cleaned.csv")

day_ves_counts <- day$Spiky.Vesicles.Observed[!is.na(day$Spiky.Vesicles.Observed)]
night_ves_counts <- night$Spiky.Vesicles.Observed[!is.na(night$Spiky.Vesicles.Observed)]
morning_ves_counts <- c(0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

day_ves_size <- day$Vesicle.Diameter..widest.[!is.na(day$Vesicle.Diameter..widest.)]
night_ves_size <- night$Vesicle.Diameter..widest.[!is.na(night$Vesicle.Diameter..widest.)]
morning_ves_size <- c(54,50,50)
all_ves_size <- c(day_ves_size, night_ves_size, morning_ves_size)
ves_day_labels <- rep("Day", length(day_ves_size))
ves_night_labels <- rep("Night", length(night_ves_size))
ves_morning_labels <- rep("Morning", 3)
ves_time_labels <- c(ves_day_labels, ves_night_labels, ves_morning_labels)

#### -- PANEL A -- ####

# Vesicle counts per period (Day 88 and Night 8 are checked against the CSVs just below; Evening 39 and
# Morning 3 are the lengths of the hardcoded diameter vectors in Panel B).
stopifnot(sum(day_ves_counts) == 88, sum(night_ves_counts) == 8)
ves_density <- data.frame(
  Period   = factor(c("Day", "Evening", "Night", "Morning"),
                    levels = c("Day", "Evening", "Night", "Morning")),
  Vesicles = c(88, 39, 8, 3),
  Volume   = c(1.601, 2.957, 4.979, 5.859)   # um^3
)
ves_density$Density <- ves_density$Vesicles / ves_density$Volume

# ---- LABELED ----
ggplot(ves_density, aes(x = Period, y = Density, fill = Period)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  labs(x = "Metabolic State",
       y = expression("Vesicles per " * mu * "m"^3 * " of matrix")) +
  scale_fill_manual(values = c("Day" = "orange", "Evening" = "coral",
                               "Night" = "slateblue", "Morning" = "lightblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# ---- BLANK (for figure layout) ----
ggplot(ves_density, aes(x = Period, y = Density, fill = Period)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("Day" = "orange", "Evening" = "coral",
                               "Night" = "slateblue", "Morning" = "lightblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank())


#### -- PANEL B -- ####
evening_ves_size <- c(62, 46, 49, 38, 61, 50, 55, 46, 39, 40,
                      62, 51, 57, 54, 57, 57, 50, 63, 46, 55,
                      48, 46, 60, 47, 53, 55, 48, 47, 41, 59,
                      54, 49, 53, 40, 46, 44, 43, 44, 45)
length(evening_ves_size)   # should be 39

morning_ves_size <- c(54, 50, 50)   # morning

all_ves_size <- c(day_ves_size, evening_ves_size, night_ves_size, morning_ves_size)
ves_labels <- c(rep("Day", length(day_ves_size)),
                rep("Evening", length(evening_ves_size)),
                rep("Night", length(night_ves_size)),
                rep("Morning", length(morning_ves_size)))

plot_data <- data.frame(
  Counts = all_ves_size,
  Period = factor(ves_labels, levels = c("Day", "Evening", "Night", "Morning"))
)

cols <- c("Day" = "orange", "Evening" = "coral",
          "Night" = "slateblue", "Morning" = "lightblue")

# ---- LABELED ----
ggplot(plot_data, aes(x = Counts, fill = Period)) +
  geom_histogram(binwidth = 10, alpha = 0.7, color = "black", linetype = "dashed") +
  geom_density(aes(y = after_stat(count) * 10), alpha = 0.4) +
  theme_classic() +
  labs(x = "Vesicle diameter (nm)", y = "Number of vesicles") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = cols) +
  coord_cartesian(xlim = c(NA, 140))

# ---- BLANK ----
ggplot(plot_data, aes(x = Counts, fill = Period)) +
  geom_histogram(binwidth = 10, alpha = 0.7) +
  geom_density(aes(y = after_stat(count) * 10), alpha = 0.4) +
  theme_classic() +
  theme(legend.text = element_blank(), legend.title = element_blank()) +
  scale_fill_manual(values = cols) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  coord_cartesian(xlim = c(NA, 140))


#### -- PANEL C -- ####
# Per-tomogram hat counts and matrix volumes (um^3): data/hat_counts_by_tomogram.csv
hat_data <- read_data("hat_counts_by_tomogram.csv")
hat_data$state <- factor(hat_data$state, levels = c("Day", "Evening", "Night", "Morning"))
hat_totals <- aggregate(cbind(hats, vol) ~ state, hat_data, sum)
hat_density <- data.frame(
  Period  = hat_totals$state,
  Density = round(hat_totals$hats / hat_totals$vol, 2)   # hats per um^3: 0.62, 8.12, 9.84, 5.29
)

cols <- c("Day" = "orange", "Evening" = "coral",
          "Night" = "slateblue", "Morning" = "lightblue")

# ---- LABELED ----
ggplot(hat_density, aes(x = Period, y = Density, fill = Period)) +
  geom_bar(stat = "identity", width = 0.4, color = "black") +
  geom_text(aes(label = Density), vjust = -0.5, size = 4) +
  labs(x = "Metabolic State",
       y = expression("Hat-structures per " * mu * "m"^3 * " of matrix")) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(fill = "none") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# ---- BLANK ----
ggplot(hat_density, aes(x = Period, y = Density, fill = Period)) +
  geom_bar(stat = "identity", width = 0.4, color = "black") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank())



# Check against your spreadsheet totals
aggregate(cbind(hats, vol) ~ state, hat_data, sum)
# Expect: Day 1 / 1.601, Evening 24 / 2.957, Night 49 / 4.979, Morning 31 / 5.859

m  <- glm.nb(hats ~ state + offset(log(vol)), data = hat_data)
m0 <- glm.nb(hats ~ 1 + offset(log(vol)), data = hat_data)
summary(m)
anova(m0, m)

pairs(emmeans(m, ~ state, type = "response", offset = 0), adjust = "holm")
#day is significant compared to everything else
