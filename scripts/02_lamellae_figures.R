# ==============================================================================
# 02_lamellae_figures.R
#
# Lamellae-width density plots (thick, thin, all) and a thick-vs-thin boxplot,
# from 322 lamellae measurements. Each figure is drawn twice: labeled, then
# BLANK (no titles/labels) for figure layout.
#
# The dashed reference lines mark the MEAN width of each class (thick ~34.3 nm,
# thin ~17.3 nm). They were previously drawn at the medians (31, 17) while the
# plots were labeled "average"; they now match the labels.
#
# Input : data/lamellae_thickthincombined.csv
# Output: figures/02_lamellae_figures.pdf
#
# The vesicle / hat-structure figures that used to share this file were
# superseded by 03 and 04 and moved to archive/02_vesicle_hat_figures_superseded.R.
#
# The figures here were rewritten by scripts/05_thylakoid_width.R
# ==============================================================================
source("scripts/00_config.R")
library(ggplot2)
library(ggpubr)   # ggdensity(); was missing in the original script
open_fig_device("02_lamellae_figures")

# Collected from 322 total lamellae measurements from 18 nitroplast tomograms
lamellae_widths <- read_data("lamellae_thickthincombined.csv")
lamellae_widths_Thin_cleaned <- na.omit(lamellae_widths$Thin)
lamellae_widths_Thick_cleaned <- na.omit(lamellae_widths$Thick)

# Class means, used for the dashed reference lines
thick_mean <- mean(lamellae_widths_Thick_cleaned)   # 34.3 nm
thin_mean  <- mean(lamellae_widths_Thin_cleaned)    # 17.3 nm

# Figures: THREE density plots (thick / thin / combined) and a boxplot of the two classes.

ggplot(lamellae_widths, aes(x = Thick)) +
  geom_density(aes(y = after_stat(count)), fill = "grey", color = "black") +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey") +
    labs(title = "Density Plot of Thick Lamellae",
       x = "Lamellae Width (nm)",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) +
  annotate("text",
           x = 41,
           y = max(density(na.omit(lamellae_widths$Thick))$y) * 10,
           label = paste0("n = ", length(lamellae_widths_Thick_cleaned)),
           hjust = 1)

# ---- DENSITY PLOTS ----
ggdensity(lamellae_widths$Thick, 
          main = "Density Plot of Thick Lamellae",
          xlab = "Lamellae Width (nm)",
          ylab = "Proportion",
          fill = "grey") +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text",
           x = 41,
           y = max(density(na.omit(lamellae_widths$Thick))$y) * 0.1,
           label = paste0("n = ", length(lamellae_widths_Thick_cleaned)),
           hjust = 1)

ggdensity(lamellae_widths$Thin,
          main = "Density Plot of Thin Lamellae",
          xlab = "Lamellae Width (nm)",
          ylab = "Proportion",
          fill = "grey") +
  geom_vline(xintercept = thin_mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text",
           x = 18,
           y = max(density(na.omit(lamellae_widths$Thin))$y) * 0.1,
           label = paste0("n = ", length(lamellae_widths_Thin_cleaned)),
           hjust = 1)

ggdensity(lamellae_widths$Combined, 
          main = "Density Plot of All Lamellae",
          xlab = "Lamellae Width (nm)",
          ylab = "Proportion",
          fill = "grey") +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_vline(xintercept = thin_mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text",
           x = 35.5,
           y = max(density(na.omit(lamellae_widths$Combined))$y) * 0.1,
           label = paste0("n = ", length(lamellae_widths$Combined)),
           hjust = 1) +
  annotate("text",
           x = 55,
           y = max(density(lamellae_widths$Combined)$y) * 0.9,
           label = "Thick-type Average",
           col = "red") + 
  annotate("text",
           x = 40,
           y = max(density(lamellae_widths$Combined)$y) * 1.05,
           label = "Thin-type Average",
           col = "blue")

# ---- LAMELLAE BOXPLOT ----
boxplot(lamellae_widths[,1:2], main = "Thick vs Thin Lamellae Types", ylab = "Lamellae Width (nm)", col=c("red","blue"))

# ---- BLANKS ----
# DENSITY PLOTS
ggdensity(lamellae_widths$Thick, 
          fill = "grey") +  
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey")

ggdensity(lamellae_widths$Thin,
          fill = "grey") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_vline(xintercept = thin_mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey")

ggdensity(lamellae_widths$Combined, 
          fill = "grey") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_vline(xintercept = thin_mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey")

#changed from density to counts
ggplot(lamellae_widths, aes(x = Combined)) +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  geom_density(aes(y = after_stat(count)), fill = "grey", color = "black") +
  geom_vline(xintercept = thick_mean, linetype = "dashed", color = "red") +
  geom_vline(xintercept = thin_mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

# LAMELLAE BOXPLOT
boxplot(lamellae_widths[,1:2], col=c("red","blue"), xlab ="", xaxt = "n")

