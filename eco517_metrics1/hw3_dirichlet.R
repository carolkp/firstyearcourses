rm(list = ls())
library(dplyr)
library(ggplot2)
set.seed(0)

####### Part 1 - Load data and organize information

# Load AK data data
load("~/Documents/R/metrics1/data/akdata.RData")
# To download the data:
#load(url("http://sims.princeton.edu/yftp/emet1_2020/kmeans/akdata.RData"))

df <- akdataf

# Create vector with cuts to be made in yob to make cohorts
first_year = min(df$yob)
number_of_cuts = 6 # 5 cohorts; count the end of fifth cohort as a cut
groups <- rep(first_year - 1, number_of_cuts)
for (i in 2:number_of_cuts){
  groups[i] <- first_year - 1 + 2*(i-1)
}

# Define categorical variable for cohorts (and equivalent integer variable)
levels_for_cohort <- c("30-31","32-33","34-35","36-37", "38-39")
df$cohort = cut(df$yob, groups, labels = levels_for_cohort )
df$cohort_no <- match(df$cohort, levels_for_cohort)

# Group data by cohort and education level
educ_counts_per_cohort <- df %>%
  group_by(cohort, educ) %>%
  tally()
# Select only no education for all cohorts
no_education <- educ_counts_per_cohort[educ_counts_per_cohort$educ == 0, ]
colnames(no_education) <- c("cohort", "educ", "count_no_education_per_cohort")
# Drop educ column
no_education <- subset(no_education, select = -educ)

# Get number of people in each cohort for the whole dataset
cohort_counts <- df %>%
  group_by(cohort) %>%
  tally()
colnames(cohort_counts) <- c("cohort", "count_people_cohort")

# Joint cohort counts with no_education counts
no_education <- left_join(
  no_education, cohort_counts, by = "cohort") 

no_education = no_education%>% 
  mutate(proportion_no_education_per_cohort = count_no_education_per_cohort / count_people_cohort)

normalization_constant = sum(no_education$proportion_no_education_per_cohort)
no_education$normalized_proportion_per_cohort = no_education$proportion_no_education_per_cohort/normalization_constant

############ Part 2 - Draws from Dirichlet

rdirichlet <- function(n, alpha) {
  m <- length(alpha)
  outmat <- matrix(0, n, m)
  for (ic in 1:m) outmat[ , ic] <- rgamma(n, alpha[ic])
  apply(outmat, 1, function(x) x / sum(x))
}

# Make draws from Dirichlet:
draws <- rdirichlet(1000, no_education$normalized_proportion_per_cohort)

draws <- t(draws)
draws <- as.data.frame(draws)
colnames(draws) <- c("X1", "X2", "X3", "X4", "X5")
draws$condition <- with(draws, ifelse((X1 > X2) & (X2 > X3) & (X3 > X4) & (X4 > X5), 1, 0)) 
mean(draws$condition)

g1 <- ggplot(no_education, aes(cohort, count_no_education_per_cohort)) + 
  geom_point(color="#1B9E77", size=2) +
  ylab("No of Obs. with no schooling") +
  xlab("Cohort") + theme_gray()
g1

g2 <- ggplot(no_education, aes(cohort, proportion_no_education_per_cohort)) + 
  geom_point(color="#E7298A", size=2) +
  ylab("Proportion of Obs. with no schooling") +
  xlab("Cohort") + theme_gray()
g2
