library("boot")

# The function to create a n*p data matrix X ----

GetMatrix <- function(n, p, c, D_distr) {
  
  # input n, p, c, D_distr
  # lambda1 = 1 + c*sqrt(p/n)
  # D_distr should be chosen from {1, 2, 3}
  # D_distr == 1 -> Dii ~ N(0, 1)
  # D_distr == 2 -> Dii ~ Unif(1/2, (sqrt(3)*sqrt(4-1/4))/2 - 1/4)
  # D_distr == 3 -> Dii ~ Exp(sqrt(2))
  
  # Z0 (n*p, i.i.d, normal)
  Z0 <- matrix(rnorm(n*p), nrow = n)
  
  # D (n*n, diag, D_distr)
  if (D_distr == 1) {
    Dii <- rnorm(n)
  } else if (D_distr == 2) {
    Dii <- runif(n, min = 1/2, max = (sqrt(3)*sqrt(4-1/4))/2 - 1/4)
  } else if (D_distr == 3) {
    Dii <- rexp(n, sqrt(2))
  }
  
  D <- diag(Dii)
  
  # V (p*p, i.i.d, normal)
  V <- matrix(rnorm(p*p), nrow = p)
  
  # Lambda (p*p, diag, lambda_1,1,1...)
  lambda1 <- 1 + c*sqrt(p/n)
  Lambda <- diag(c(lambda1, rep(1, p-1)))
  
  # The n*p data matrix X
  X <- D %*% Z0 %*% V %*% Lambda %*% solve(V)
  
  return(X)
  
}


# Create a sample matrix X with n = 300, p = 3, lambda_1 = 1 (c = 0), D_ii ~ Normal ----
n <- 300
p <- 3
c <- 0
X <- GetMatrix(n = n, p = p, c = c, D_distr = 1)


# Estimate bias and standard error using "boot" package ----

## The function to get the top eigenvalue of X's sample covariance matrix
GetLambdaOneHat <- function(X, indices) {
  
  data <- X[indices, ]
  lambda1_hat <- max(eigen(cov(data))$values)
  return(lambda1_hat)
  
}

GetLambdaOneHat(X)

## bootstrap estimate of lambda1 B = 200
B <- 200
boot(data = X, statistic = GetLambdaOneHat, R = B)

# Estimate bias and standard error step by step ----

ind_sample <- matrix(NA, nrow = n, ncol = B)
sample <- vector(mode = "list")  ## B bootstrap samples with sample size n
Lambda1_hat_boot <- NA  ## B bootstrap estimates of Lambda_hat
for (i in 1:B) {
  ind_sample[, i] <- sample(1:n, size = n, replace = T)
  sample[[i]] <- X[ind_sample[, i], ]
  Lambda1_hat_boot[i] <- GetLambdaOneHat(sample[[i]])
}

Lambda1_hat <- GetLambdaOneHat(X)  ## Lambda_hat of sample covariance matrix of X

bias_boot <- mean(Lambda1_hat_boot) - Lambda1_hat
std_err_boot <- sqrt( sum((Lambda1_hat_boot - mean(Lambda1_hat_boot))^2) / (B - 1))


# The function to perform a complete simulation ----

## input n, p, c, D_distr, B
## output a list of 
## $Lambda1_hat: Lambda_hat of sample covariance matrix of X
## $Lambda1_hat_boot: (vector) B bootstrap results of original Lambda1_hat 
## $bias_boot
## $std_err_boot
## $r: p/n
## $true_Lambda1
## $true_bias
## $true_std_err

PerfSimulation <- function(n, p, c, D_distr, B, ...) {
  
  X <- GetMatrix(n, p, c, D_distr)
  boot_result <- boot(data = X, statistic = GetLambdaOneHat, R = B)
  Lambda1_hat <- boot_result$t0
  Lambda1_hat_boot <- boot_result$t
  bias_boot <- mean(Lambda1_hat_boot) - Lambda1_hat
  std_err_boot <- sqrt( sum((Lambda1_hat_boot - mean(Lambda1_hat_boot))^2) / (B - 1))
  r <- p/n
  true_Lambda1 <- 1 + c*sqrt(r)
  true_bias <- true_Lambda1 - Lambda1_hat
  true_std_err <- "???"
  outcome <- list(Lambda1_hat = Lambda1_hat, 
                  Lambda1_hat_boot = Lambda1_hat_boot, 
                  bias_boot = bias_boot, 
                  std_err_boot = std_err_boot, 
                  r = r, 
                  true_Lambda1 = true_Lambda1, 
                  true_bias = true_bias, 
                  true_std_err = true_std_err)
  return(outcome)
  
}


# 500 simulations ----

## n = 300
## p = 3, 30, 100, 150
## c = 0, 3, 11, 50, 100
## Dii ~ N(0, 1), Dii ~ Unif(1/2, (sqrt(3)*sqrt(4-1/4))/2 - 1/4), Dii ~ Exp(sqrt(2))
## B = 200

## n = 300, p = 3, c = 0, Dii ~ N(0, 1), B = 200, 500 simulations
simu <- vector(mode = "list")
simu_vec <- matrix(NA, nrow = 500, ncol = 4)
for(i in 1:500) {
  simu[[i]] <- PerfSimulation(n = 300, p = 3, c = 0, D_distr = 1, B = 200)
  simu_vec[i, ] <- c(simu[[i]]$r, simu[[i]]$Lambda1_hat, simu[[i]]$bias_boot, simu[[i]]$std_err_boot)
}
colnames(simu_vec) <- c("r", "Lambda1_hat", "bias_boot", "std_err_boot")
head(simu_vec)
