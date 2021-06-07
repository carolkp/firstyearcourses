library(dplyr)
library(ggplot2)

load("~/R/scripts/econometric_theory/hw1/akdata.RData")
# To download the data:
#load(url("http://sims.princeton.edu/yftp/emet1_2020/kmeans/akdata.RData"))

#ak_data_matrix <-  data.matrix(akdataf[ , 1:2])
akdataf <- akdataf[ , 1:2]
akdataf$wage_decile <- ntile(akdataf$logwage, 10)  

wagedecile_group_average <- akdataf %>% 
  group_by(wage_decile) %>% 
  summarise(educ = mean(educ))

g <- ggplot(wagedecile_group_average, aes(wage_decile, educ)) + 
  geom_point(color="#1B9E77", size=2) +
  ylab("Expected Years of Education") +
  xlab("Log of Wage Decile") + scale_x_discrete(limits=1:10)
g


educ_group_average <- akdataf %>% 
  group_by(educ) %>% 
  summarise(logwage_decile = mean(wage_decile))

g2 <- ggplot(educ_group_average, aes(educ, logwage_decile)) + 
  geom_point(color="#E7298A", size=2) +
  ylab("Expected Log of Wage Decile") +
  xlab("Years of Education") + scale_x_discrete(limits=seq(from=0, to=20, by=2))
g2
