rm(list = ls())
library(dplyr)
library(ggplot2)
library(MCMCpack)
set.seed(18022020)

f = function(x) {
  !any(diff(x) < 0)
} 

# Load AK data data
load("~/Documents/R/metrics1/data/akdata.RData")
# To download the data:
#load(url("http://sims.princeton.edu/yftp/emet1_2020/kmeans/akdata.RData"))
df <- akdataf

### First Approach
# Idea is to treat sigma_s as known for every S = years of education and 
# draw mu from the posterior distribution given the sigma_s, which is normal

# Estimate variance, mean and count observations for each group
groupedby_schooling <- df %>%
  group_by(educ) %>%
  summarise(std_deviation = sd(logwage), variance = var(logwage), sample_mean = mean(logwage), count = n()) 


# Given sigma_s, the posterior distribution of mu_s is normal with mean equal 
# to the sample mean and variance sigma_s / N_s, N_s is the number of obs 
# in group s
groupedby_schooling$posterior_var <- groupedby_schooling$variance / groupedby_schooling$count
groupedby_schooling$posterior_sd <- (groupedby_schooling$std_deviation)/sqrt(groupedby_schooling$count)

# Now we'll make the draws from all distributions
n <- 1000
m <- 21
first_approach <- matrix(0, n, m)

for (i in 1:m)
  first_approach[, i] <- rnorm(n, groupedby_schooling$sample_mean[i], groupedby_schooling$posterior_sd[i])


first_df <- data.frame(first_approach)
first_df$end_highschool <- first_df["X13"] - first_df["X12"]
first_df$notend_highschool <- first_df["X12"] - first_df["X11"]
first_df$all_college <- first_df["X17"] - first_df["X13"]
first_df$all_highschool <- first_df["X13"] - first_df["X9"]
first_approach_tests <- matrix(0, n, 5)

## Evalutate the probabilities:
# An additional year of schooling always corresponds to a higher expected log wage, whatever the year
first_approach_tests[, 1] <- apply(first_df, 1, f)
# Staying the last year of high school, and thereby graduating, gives a greater boost to expected income than staying the previous year
first_approach_tests[, 2] <- with(first_df, ifelse(end_highschool > notend_highschool, TRUE, FALSE)) 
# Every year of schooling from 5th grade through 4 year college graduation increases expected log wage
first_approach_tests[, 3] <- apply(first_df[,5:17], 1, f)
# Every year of schooling beyond college graduation adds to expected log wage.
first_approach_tests[, 4] <- apply(first_df[,17:21], 1, f)
# High school adds less to expected income than four years of college does
first_approach_tests[, 5] <- with(first_df, ifelse((all_college > all_highschool), TRUE, FALSE)) 

print(colMeans(first_approach_tests))

### Second Approach
# Use the fact that the joint probability of mu, sigma^2 is normal- inverse -gamma
# First, we draw 1000 values of d_s (estimate for sigma_s) from the Inverse Gamma
# The, for each of these values, we draw a mu from a normal distribution
second_approach_variance <- matrix(0, n, m)
second_approach_mean <- matrix(0, n, m)

for (i in 1:m){
  shape_parameter <- (groupedby_schooling$count[i] - 3)/2
  scale_parameter <- (groupedby_schooling$count[i]*groupedby_schooling$variance[i]/2)^(-1)
  #second_approach_variance[, i] <- rinvgamma(n, shape_parameter, scale = scale_parameter)
  second_approach_variance[, i] <- rinvgamma(n, shape_parameter, scale = scale_parameter)
  second_approach_variance[, i] <- sqrt(second_approach_variance[, i])
  # pass vector with standard deviations drawn from IG, but same mean
  second_approach_mean[, i] <- rnorm(n, groupedby_schooling$sample_mean[i], second_approach_variance[, i])}

second_df <- data.frame(second_approach_mean)
second_df$end_highschool <- second_df["X13"] - second_df["X12"]
second_df$notend_highschool <- second_df["X12"] - second_df["X11"]
second_df$all_college <- second_df["X17"] - second_df["X13"]
second_df$all_highschool <- second_df["X13"] - second_df["X9"]
second_approach_tests <- matrix(0, n, 5)

## Evalutate the probabilities:
# An additional year of schooling always corresponds to a higher expected log wage, whatever the year
second_approach_tests[, 1] <- apply(second_df, 1, f)
# Staying the last year of high school, and thereby graduating, gives a greater boost to expected income than staying the previous year
second_approach_tests[, 2] <- with(second_df, ifelse(end_highschool > notend_highschool, TRUE, FALSE)) 
# Every year of schooling from 5th grade through 4 year college graduation increases expected log wage
second_approach_tests[, 3] <- apply(second_df[,5:17], 1, f)
# Every year of schooling beyond college graduation adds to expected log wage.
second_approach_tests[, 4] <- apply(second_df[,17:21], 1, f)
# High school adds less to expected income than four years of college does
second_approach_tests[, 5] <- with(second_df, ifelse((all_college > all_highschool), TRUE, FALSE)) 

print(colMeans(second_approach_tests))
