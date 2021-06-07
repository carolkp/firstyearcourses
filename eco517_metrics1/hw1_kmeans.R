library(RColorBrewer)
library(ggplot2)
library(dplyr)

load("~/R/scripts/econometric_theory/hw1/akdata.RData")
# To download the data:
#load(url("http://sims.princeton.edu/yftp/emet1_2020/kmeans/akdata.RData"))

#Get first two columns
ak_data_matrix <-  data.matrix(akdataf[ , 1:2])
akv1 <- ak_data_matrix

# Calculate standard deviation
standard_deviation_akv1 <- apply(akv1, 2, sd)
## normalizes so variances are the same
akv1 <- akv1 %*% diag(1/standard_deviation_akv1)

## Run k-means algorithm
koutput2 <- kmeans(akv1, 2, iter.max=100)
koutput3 <- kmeans(akv1, 3, iter.max=100)
koutput4 <- kmeans(akv1, 4, iter.max=100)
koutput5 <- kmeans(akv1, 5, iter.max=100)

unnormalized_koutput3 <- kmeans(ak_data_matrix , 3, iter.max=100)

# Add column for clusters
akdataf$cluster2 <- koutput2$cluster
akdataf$cluster3 <- koutput3$cluster
akdataf$cluster4 <- koutput4$cluster
akdataf$cluster5 <- koutput5$cluster

cbbPalette <- c("#FF6666", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")

#Calculate average within clusters

groups2_average <- akdataf %>% 
  group_by(cluster2) %>% 
  summarise(educ = mean(educ),
            logwage  = mean(logwage))

groups3_average <- akdataf %>% 
  group_by(cluster3) %>% 
  summarise(educ = mean(educ),
            logwage  = mean(logwage))

groups4_average <- akdataf %>% 
  group_by(cluster4) %>% 
  summarise(educ = mean(educ),
            logwage  = mean(logwage))

groups5_average <- akdataf %>% 
  group_by(cluster5) %>% 
  summarise(educ = mean(educ),
            logwage  = mean(logwage))

# Scatter plots

## Unnormalized with 3 clusters
unnormalized_plot_3groups <- ggplot(akdataf, aes(educ, logwage, color = factor(unnormalized_koutput3$cluster))) + 
  geom_point(alpha = 0.4) + 
  xlab("years of education") + 
  ylab("log of wage") +
  theme_bw() + 
  scale_colour_manual(values=cbbPalette, name = "Groups") 

print(unnormalized_plot_3groups)


## Normalized for 2-5 clusters
plot_2groups <- ggplot(akdataf, aes(educ, logwage, color = factor(koutput2$cluster))) + 
  geom_point(alpha = 0.4) + 
  xlab("years of education") + 
  ylab("log of wage") +
  theme_bw() + 
  scale_colour_manual(values=cbbPalette, name = "Groups") 

print(plot_2groups)

plot_3groups <- ggplot(akdataf, aes(educ, logwage, color = factor(koutput3$cluster))) + 
  geom_point(alpha = 0.4) + 
  xlab("years of education") + 
  ylab("log of wage") +
  theme_bw() + 
  scale_colour_manual(values=cbbPalette, name = "Groups") 

print(plot_3groups)

plot_4groups <- ggplot(akdataf, aes(educ, logwage, color = factor(koutput4$cluster))) + 
  geom_point(alpha = 0.4) + 
  xlab("years of education") + 
  ylab("log of wage") +
  theme_bw() + 
  scale_colour_manual(values=cbbPalette, name = "Groups") 

print(plot_4groups)

plot_5groups <- ggplot(akdataf, aes(educ, logwage, color = factor(koutput5$cluster))) + 
  geom_point(alpha = 0.4) + 
  xlab("years of education") + 
  ylab("log of wage") +
  theme_bw() + 
  scale_colour_manual(values=cbbPalette, name = "Groups") 

print(plot_5groups)

## Histograms
ggplot(akdataf, aes(x=logwage)) + 
  geom_histogram(aes(y=..density..), fill="#FF6666",  color="#860303") + 
  theme_bw() +
  xlab("Log of wage") 

ggplot(akdataf, aes(x=educ)) + 
  geom_histogram(aes(y=..density..), fill="#8CFCA6",  color="#66A475") + 
  theme_bw() +
  xlab("Education") 