# PSET 11 - Metrics - Part 2
# Efficient GMM (multiple and single equation)
rm(list = ls())
set.seed(6)
library(gmm)

# Numberof eq
M <- 3
# Number of obs
n <- 1000

# Coefficients
delta0 <- rep(1,M)
delta1 <- rep(1,M)
alpha1 <- rep(1,M)
alpha2 <- rep(1,M)
beta <- matrix(0, 3, M)
alpha3 <- matrix(0, 3, M)


w_1 <- rnorm(n, mean = 0, sd = 1)
w_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)

v_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
v_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
v_3 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)

x_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
x_3 <-matrix(data = runif(n*M, min = 0, max = 3), ncol = 3, nrow = n)
x_2 <- w_1 + w_2

u_1 <- rnorm(n, mean = 0, sd = 1)
u_2 <- rnorm(n, mean = 0, sd = 1)

z <- x_1*alpha1 + x_2*alpha2 + x_3%*%alpha3 + v_1 + v_3 + u_2
y <- delta0 + z*delta1 + x_3%*%beta + v_1 + v_3 + u_1

# Data: y, z, x_1, x_2, x_3
data <- data.frame(cbind(y, z, x_1, x_2, x_3))

for (m in 1:M){
  names(data)[m] <- paste0("y",m)
  names(data)[m + M] <- paste0("z",m)
  names(data)[m + 2*M] <- paste0("x1",m)
  names(data)[m + 3*M] <- paste0("x2",m)
}
x3names <- paste0("x3", 1:3, sep="")
ncolumns <- dim(data)[2]
names(data)[(ncolumns-2):ncolumns] <- x3names

#4a) Regress y on the x's with single equation
y1 <- y[,1]
y2 <- y[,2]
y3 <- y[,3]
singleeq1 <- gmm(y1 ~ z[,1] + x_3, ~ x_1[,1] + x_2[, 1] + x_3, type="twoStep", wmatrix="optimal")
summary(singleeq1)
singleeq2 <- gmm(y2 ~ z[,2] + x_3, ~ x_1[,2] + x_2[, 2] + x_3, type="twoStep", wmatrix="optimal")
summary(singleeq2)
singleeq3 <- gmm(y3 ~ z[,3] + x_3, ~ x_1[,3] + x_2[, 3] + x_3, type="twoStep", wmatrix="optimal")
summary(singleeq3)


# 4b

eq1 <- y1 ~ z[,1] + x_3
eq2 <- y2 ~ z[,2] + x_3
eq3 <- y3 ~ z[,3] + x_3
instrum <- list( ~ x_1[,1] + x_2[, 1] + x_3,  ~ x_1[,2] + x_2[, 2] + x_3, ~ x_1[,3] + x_2[, 3] + x_3 )
formulas <- list(eq1, eq2, eq3)

multeq <- sysGmm(formulas, instrum, wmatrix="optimal")
summary(multeq)

# 5
simulations = 500

coefs_single <- matrix(0, simulations, 15)
# pvalues_single <- matrix(0, simulations, 15)

coefs_mult <- matrix(0, simulations, 15)
# pvalues_mult <- matrix(0, simulations, 15)

for (s in 1:simulations){
  if (s %% 10 == 0){
    print(s)
  }
  # generate data
  w_1 <- rnorm(n, mean = 0, sd = 1)
  w_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  
  v_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  v_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  v_3 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  
  x_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  x_3 <-matrix(data = runif(n*M, min = 0, max = 3), ncol = 3, nrow = n)
  x_2 <- w_1 + w_2
  
  u_1 <- rnorm(n, mean = 0, sd = 1)
  u_2 <- rnorm(n, mean = 0, sd = 1)
  
  z <- x_1*alpha1 + x_2*alpha2 + x_3%*%alpha3 + v_1 + v_3 + u_2
  y <- delta0 + z*delta1 + x_3%*%beta + v_1 + v_3 + u_1
  
  
  y1 <- y[,1]
  y2 <- y[,2]
  y3 <- y[,3]
  
  # Single equation gmm
  singleeq1 <- gmm(y1 ~ z[,1] + x_3, ~ x_1[,1] + x_2[, 1] + x_3, type="twoStep", wmatrix="optimal")
  singleeq2 <- gmm(y2 ~ z[,2] + x_3, ~ x_1[,2] + x_2[, 2] + x_3, type="twoStep", wmatrix="optimal")
  singleeq3 <- gmm(y3 ~ z[,3] + x_3, ~ x_1[,3] + x_2[, 3] + x_3, type="twoStep", wmatrix="optimal")
  # save coefficients and pvalues
  coefs_single[s, ] <- cbind(singleeq1$coefficients, singleeq2$coefficients, singleeq3$coefficients)
  # pvalues_single[s, ] <- cbind(summary(singleeq1)$coefficients[,"Pr(>|t|)"], summary(singleeq2)$coefficients[,"Pr(>|t|)"], summary(singleeq3)$coefficients[,"Pr(>|t|)"])

  # Multiple equation gmm
  eq1 <- y1 ~ z[,1] + x_3
  eq2 <- y2 ~ z[,2] + x_3
  eq3 <- y3 ~ z[,3] + x_3
  instrum <- list( ~ x_1[,1] + x_2[, 1] + x_3,  ~ x_1[,2] + x_2[, 2] + x_3, ~ x_1[,3] + x_2[, 3] + x_3 )
  formulas <- list(eq1, eq2, eq3)
  
  multeq <- sysGmm(formulas, instrum, wmatrix="optimal")
  
  coefs_mult[s, ] <- cbind(multeq$coefficients$System_1, multeq$coefficients$System_2, multeq$coefficients$System_3)
  # pvalues_mult[s, ] <- cbind(summary(multeq)$coefficients[[1]][,"Pr(>|t|)"], summary(multeq)$coefficients[[2]][,"Pr(>|t|)"], summary(multeq)$coefficients[[3]][,"Pr(>|t|)"])
  
  }

true_coefs <- c(1, 1, 0, 0, 0)
for (i in 1:5){
  print("Single Eq. RMSE:")
  print(sqrt(mean((coefs_single[,i] - true_coefs[i])^2)))
  print("")
  print("Multiple Eq. RMSE:")
  print(sqrt(mean((coefs_mult[,i] - true_coefs[i])^2)))
}


# 6
rm(list = ls())
# Number of eq
M <- 25
n <- 1000
simulations <- 500
ncolumns <- 4*M + 3

data <-data.frame(matrix(0, n, ncolumns))

for (m in 1:M){
  names(data)[m] <- paste0("y",m)
  names(data)[m + M] <- paste0("z",m)
  names(data)[m + 2*M] <- paste0("x1",m)
  names(data)[m + 3*M] <- paste0("x2",m)
}
x3names <- paste0("x3", 1:3, sep="")
ncolumns <- dim(data)[2]
names(data)[(ncolumns-2):ncolumns] <- x3names

# True Coefficients
delta0 <- rep(1,M)
delta1 <- rep(1,M)
alpha1 <- rep(1,M)
alpha2 <- rep(1,M)
beta <- matrix(0, 3, M)
alpha3 <- matrix(0, 3, M)

# Create matrices to save coefficients
coefs_single <- matrix(0, simulations, 5)
pvalues_single <- matrix(0, simulations, 5)

coefs_mult <- matrix(0, simulations, 5)
pvalues_mult <- matrix(0, simulations, 5)

for (s in 1:simulations){
  if (s %% 10 == 0){
    print(s)
  }
  # generate data
  w_1 <- rnorm(n, mean = 0, sd = 1)
  w_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  
  v_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  v_2 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  v_3 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  
  x_1 <- matrix(data = rnorm(n*M, mean = 0, sd = 1), ncol = M, nrow = n)
  x_3 <-matrix(data = runif(n*3, min = 0, max = 3), ncol = 3, nrow = n)
  x_2 <- w_1 + w_2
  
  u_1 <- rnorm(n, mean = 0, sd = 1)
  u_2 <- rnorm(n, mean = 0, sd = 1)
  
  z <- x_1*alpha1 + x_2*alpha2 + x_3%*%alpha3 + v_1 + v_3 + u_2
  y <- delta0 + z*delta1 + x_3%*%beta + v_1 + v_3 + u_1
  
  y1 <- y[,1]
  # y2 <- y[,2]
  # y3 <- y[,3]
  
  
  # Data: y, z, x_1, x_2, x_3
  data[1:M] <- y
  data[(M+1):(2*M)]<- z
  data[(2*M + 1):(3*M)] <- x_1
  data[(3*M + 1):(4*M)] <- x_2
  data[(ncolumns - 2):ncolumns] <- x_3
  
  for (m in 1:M){
    names(data)[m] <- paste0("y",m)
    names(data)[m + M] <- paste0("z",m)
    names(data)[m + 2*M] <- paste0("x1",m)
    names(data)[m + 3*M] <- paste0("x2",m)
  }
  x3names <- paste0("x3", 1:3, sep="")
  ncolumns <- dim(data)[2]
  names(data)[(ncolumns-2):ncolumns] <- x3names
  
  # Single equation gmm
  singleeq1 <- gmm(y1 ~ z1 + x31 + x32 + x33, ~ x11 + x21 + x31 + x32 + x33, type="twoStep", wmatrix="optimal", data = data)
  #singleeq2 <- gmm(y2 ~ z[,2] + x_3, ~ x_1[,2] + x_2[, 2] + x_3, type="twoStep", wmatrix="optimal")
  #singleeq3 <- gmm(y3 ~ z[,3] + x_3, ~ x_1[,3] + x_2[, 3] + x_3, type="twoStep", wmatrix="optimal")
  # save coefficients and pvalues
  coefs_single[s, ] <- singleeq1$coefficients
  pvalues_single[s, ] <- summary(singleeq1)$coefficients[,"Pr(>|t|)"]
  
  # Multiple equation gmm
  instruments <- list()
  formulas <- list()
  
  rhs2 <-paste(x3names, collapse="+")
  
  for (m in 1:M){
    lhs <- paste0("y", m, " ~ ")
    rhs1 <- paste0("z", m, "+")
    
    eq <- as.formula(paste0(lhs, rhs1, rhs2))
    this_inst <- as.formula(paste0("~ x1", 1, "+", "x2", 1, "+", rhs2))
    
    formulas <- c(formulas, list(eq))
    instruments <- c(instruments, list(this_inst))
    }
  
  multeq <- sysGmm(formulas, instruments, wmatrix="optimal", data = data)
  # save coefficients and pvalues
  coefs_mult[s, ] <- multeq$coefficients$System_1
  pvalues_mult[s, ] <- summary(multeq)$coefficients[[1]][,"Pr(>|t|)"]
  
}

true_coefs <- c(1, 1, 0, 0, 0)
for (i in 1:5){
  print("Single Eq. RMSE:")
  print(sqrt(mean((coefs_single[,i] - true_coefs[i])^2)))
  print("")
  print("Multiple Eq. RMSE:")
  print(sqrt(mean((coefs_mult[,i] - true_coefs[i])^2)))
}