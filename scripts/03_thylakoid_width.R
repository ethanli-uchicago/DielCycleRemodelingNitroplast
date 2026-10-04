# ==============================================================================
# 03_thylakoid_width.R
#
# Thylakoid widths (nm) by time point (Day / Night / Morning / Evening):
# histograms of counts, a logistic regression and Welch t-test of Day vs Night,
# and a check of how many widths exceed 30 nm per group.
#
# Inputs : none -- width measurements are hardcoded below (see data/README.md)
# ==============================================================================
source("scripts/00_config.R")
library(ggplot2)
library(ggpubr)   # was loaded mid-script in the original
open_fig_device("05_thylakoid_width")

day <- c(28, 16, 38, 24, 24, 30, 20, 25, 30, 18, 30, 28, 12, 15, 12, 12, 15, 14, 14, 15,
         17, 17, 16, 19, 16, 23, 22, 16, 12, 20, 22, 21, 12, 15, 11, 14, 17, 16, 12, 18,
         15, 15, 13, 16, 11, 15, 15, 14, 12, 16, 15, 15, 17, 24, 17, 18, 19, 20, 22, 18,
         19, 19, 18, 25, 21, 19, 12, 16, 16, 14, 13, 11, 10, 20, 15, 18, 19, 14, 13, 15,
         14)

night <- c(22, 19, 15, 23, 16, 19, 19, 19, 18, 22,
           25, 29, 18, 28, 43, 28, 29, 26, 18, 27,
           21, 13, 14, 35, 17, 32, 15, 21, 32, 28,
           23, 22, 18, 17, 18, 20, 16, 17, 17, 16,
           18, 18, 16, 16, 14, 18, 16, 19, 19, 20,
           43, 40, 18, 23, 21, 18, 14, 15, 13, 23,
           14, 18, 17, 20, 16, 20, 15, 13, 19, 27,
           19, 18, 16, 20, 20, 20, 17, 15, 24, 18,
           18, 24, 17, 16, 13, 21, 24, 25, 17, 17,
           22, 17, 17, 16, 17, 16, 14, 16, 15, 14,
           20, 17, 24, 20, 17, 14, 17, 25, 19, 36,
           30, 15, 15)

morning <- c(14, 20, 33, 33, 14, 43, 24, 25, 34, 39,
             39, 34, 34, 42, 36, 32, 17, 32, 33, 44,
             24, 31, 19, 15, 20, 19, 19, 21, 18, 15,
             17, 16, 15, 15, 22, 19, 23, 23, 19, 24,
             28, 25, 36, 25, 34, 24, 24, 10, 11, 18,
             13, 17, 15, 12, 25, 25, 56, 40, 28, 26,
             43, 20, 54, 12, 20, 25, 32, 30, 36, 43,
             58, 61, 36, 40, 29, 37, 25, 46, 38, 23,
             28, 35, 24, 24, 32, 35, 26, 139, 29, 32,
             22, 18, 38, 25, 30, 23, 50, 34, 28, 27,
             42, 44, 89, 10, 16, 11, 22, 14, 24, 12,
             16, 17, 26, 16, 9, 15, 15, 13, 16, 15,
             21, 14, 20, 28, 22)

evening <- c(18, 14, 16, 15, 15, 14, 17, 20, 18, 23,
             17, 25, 23, 22, 32, 21, 24, 22, 29, 37,
             31, 33, 22, 47, 27, 43, 117, 40, 23, 27,
             26, 33, 27, 25, 23, 34, 32, 19, 21, 19,
             19, 30, 24, 25, 27, 23, 35, 30, 18, 29,
             21, 15, 14, 16, 13, 14, 23, 9, 22, 25,
             29, 27, 30, 16, 22, 14, 21, 22)

# Combined: all 387 measurements. (Originally retyped by hand; verified identical to this concatenation.)
all_thylakoid <- c(day, night, morning, evening)


# ---- COUNT VERSION (matches your "changed from density to counts" plot) ----
ggplot(data.frame(width = all_thylakoid), aes(x = width)) +
  geom_histogram(aes(y = after_stat(count)), fill = "grey", color = "black") +
  scale_x_continuous(limits = c(0, 150)) +
  #geom_vline(xintercept = 31, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 17, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

ggplot(data.frame(width = day), aes(x = width)) +
  geom_histogram(aes(y = after_stat(count)), fill = "orange", color = "black") +
  scale_x_continuous(limits = c(0, 150)) +
  #geom_vline(xintercept = 31, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 17, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

ggplot(data.frame(width = night), aes(x = width)) +
  geom_histogram(aes(y = after_stat(count)), fill = "slateblue", color = "black") +
  scale_x_continuous(limits = c(0, 150)) +
  #geom_vline(xintercept = 31, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 17, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

ggplot(data.frame(width = morning), aes(x = width)) +
  geom_histogram(aes(y = after_stat(count)), fill = "lightblue", color = "black") +
  scale_x_continuous(limits = c(0, 150)) +
  #geom_vline(xintercept = 31, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 17, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

ggplot(data.frame(width = evening), aes(x = width)) +
  geom_histogram(aes(y = after_stat(count)), fill = "coral", color = "black") +
  scale_x_continuous(limits = c(0, 150)) +
  #geom_vline(xintercept = 31, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 17, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())




# Statistical test to show if there is no relationship between thylakoid size and state
# Ultimately, the statistical significance (OR = 1.09 per nm, 95% CI 1.03-1.16, p = 0.0046) of this test
# is trivial because length measurement is inconsistent due to pixel sizes, binning, and noise
dat <- data.frame(
  width = c(day, night),
  state = factor(rep(c("Day", "Night"), c(length(day), length(night))),
                 levels = c("Day", "Night"))
)

m <- glm(state ~ width, data = dat, family = binomial)
summary(m)
exp(cbind(OR = coef(m), confint(m)))   # odds ratios

### how well does it discriminate? ###
t.test(log(day), log(night))
# There is a difference between night and day thylakoids, albeit it is so small that it is hard to make a statement


# Below is just checking thresholds

count <- 0
for(i in 1:length(day)){
  if(day[i] > 30){
    count <- count + 1
  } 
}
count

count <- 0
for(i in 1:length(morning)){
  if(morning[i] > 30){
    count <- count + 1
  } 
}
count

count <- 0
for(i in 1:length(night)){
  if(night[i] > 30){
    count <- count + 1
  } 
}
count

count <- 0
for(i in 1:length(evening)){
  if(evening[i] > 30){   # BUGFIX: original tested night[i] here
    count <- count + 1
  } 
}
count