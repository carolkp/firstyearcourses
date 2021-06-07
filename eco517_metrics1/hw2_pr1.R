library(pracma)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# Function to calculate expected utility given probability of 
# rain and utility in each state
calculate_curve <- function(x, utility_rain, utility_sun){
  y <- utility_rain*x + utility_sun*(1-x)
  return(y)
}

number_of_points <- 5000
umbrella_options <- c("A", "B", "C", "D", "E", "F", "None")

# Create table with utility values for each choice
utility <- data.frame("Umbrella"=umbrella_options, "Rain" = 
                        c(1, 0.6, 0.5, 0.3, 0.25, 0.1, -1), 
                      "Sun" = c(1, 2, 2.6, 3, 3, 4, 6))

plots_values <- matrix(ncol=nrow(utility), nrow=number_of_points)
colnames(plots_values) <- umbrella_options

# Grid on interval [0,1] corresponding to probabilities of rain
x_points <- linspace(0, 1, n= number_of_points)

#For each choice, calculate expected utility for all x values
for (row in 1:nrow(utility)){
  plots_values[,row] <- calculate_curve(x_points, utility[row, 2], 
                                        utility[row, 3])
}
plots_df <- as.data.frame(plots_values)

#Find max expected utility for each value of x
plots_df$max <- apply(plots_values[, 1:nrow(utility)], 1, max)
plots_df$x_values <- x_points

# Plot results
plots_df <- melt(plots_df ,  id.vars = 'x_values', variable.name = 'umbrella')

my_dark_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
                     "#E6AB02", "#A6761D", "#000000")

g <- ggplot(plots_df, aes(x_values, value)) + geom_line(aes(colour=umbrella), 
                                                        size=0.8) + 
    xlab("Probability of Rain") + 
  ylab("Utility for each choice") + 
  scale_colour_manual(values=my_dark_palette, name = "Umbrella") +
  ggtitle("Expected Utility of Each Choice vs. \n Probability of Rain")
g

# See which choices were optimal at some point
print("The umbrella choices that maximize expected utility for some 
      value of probability of rain are:")
print(unique(colnames(plots_values)[apply(plots_values,1,which.max)]))

