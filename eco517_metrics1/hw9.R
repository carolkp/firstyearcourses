## PSET 9 - Metrics I

## Problem 1
# Compare coefficient test using standard std. error and heteroskedasticity robust 
rm(list = ls())
library(car)
library(sandwich)

set.seed(1)
n = 200
simulations = 1000

# Generate data 
x = matrix(data = runif(n*simulations, min = 0, max = 1), ncol = simulations, nrow = n)
v = matrix(data = rnorm(n*simulations, mean = 0, sd = 1), ncol = simulations, nrow = n)
eps = x*v
y = x + eps

reject_usual <- numeric(simulations)
pvalue_usual <- numeric(simulations)
reject_robust <- numeric(simulations)
pvalue_robust <- numeric(simulations)

for (i in 1:simulations){
  df <- data.frame(cbind(y[, i], x[, i], eps[, i]))
  names(df) <- c("dependent", "independent", "error")
  reg <- lm(dependent ~ 0 + independent, data = df)
  test <- linearHypothesis(reg, "independent=1", vcov =vcovHC(reg, type = "const"))
  test_robust <- linearHypothesis(reg, "independent=1", vcov =vcovHC(reg, type = "HC0"))
  pvalue_usual[i] <- test[2,4]
  reject_usual[i] <- ifelse((test[2,4] < 0.05), TRUE, FALSE)
  pvalue_robust[i] <- test_robust[2,4]
  reject_robust[i] <- ifelse((test_robust[2,4] < 0.05), TRUE, FALSE)
}

mean(reject_usual)
mean(reject_robust)


# Problem 2

rm(list = ls())
set.seed(0)
samples <- 100
deg_freedom <- samples - 2
j_values <- c(1, 10, 100, 1000, 10000)
number_regs <- length(j_values) + 1

pvalue_usual <- numeric(number_regs)
pvalue_robust <-numeric(number_regs)
std_error <- numeric(number_regs)
t_statistic_homosk <- numeric(number_regs)
t_statistic_heterosk <- numeric(number_regs)

x <- rnorm(samples, mean = 0, sd = 1)
eps <- rnorm(samples, mean = 0, sd = 1)
y <- 1 + x + eps

df <- data.frame(cbind(y, x, eps))
names(df) <- c("dependent", "independent", "error")
regression <- lm(dependent ~ independent, data = df)

#cov_matrix_homosk = vcovHC(regression, type = "const")
cov_matrix_heterosk = vcovHC(regression, type = "HC0")
std_error[number_regs] <- coef(summary(regression))[2,2]

t_statistic_homosk[number_regs] <- abs(regression$coefficients[2] - 1)/std_error[number_regs]
t_statistic_heterosk[number_regs] <- abs(regression$coefficients[2] - 1)/sqrt(cov_matrix_heterosk[2,2])

pvalue_usual[number_regs] <- 2*pt(t_statistic_homosk[number_regs], regression$df.residual, lower.tail = FALSE)
pvalue_robust[number_regs] <- 2*pt(t_statistic_heterosk[number_regs], regression$df.residual, lower.tail= FALSE)


for (i in 1:length(j_values)){
  x[1] <- j_values[i]
  y[1] <- 1 + x[1] + eps[1]
  df <- data.frame(cbind(y, x, eps))
  names(df) <- c("dependent", "independent", "error")
  regression <- lm(dependent ~ independent, data = df)
  
  std_error[i] <- coef(summary(regression))[2,2]
  cov_matrix_heterosk = vcovHC(regression, type = "HC0")
  print(regression$coefficients[2])
  t_statistic_homosk[i] <- abs(regression$coefficients[2] - 1)/std_error[i]
  t_statistic_heterosk[i] <- abs(regression$coefficients[2] - 1)/sqrt(cov_matrix_heterosk[2,2])
  
  pvalue_usual[i] <- 2*pt(t_statistic_homosk[i], regression$df.residual, lower.tail = FALSE)
  pvalue_robust[i] <- 2*pt(t_statistic_heterosk[i], regression$df.residual, lower.tail= FALSE)
}

print("Pvalues:")
print(pvalue_usual)
print(pvalue_robust)


print("t-Statistics:")
print(t_statistic_homosk)
print(t_statistic_heterosk)