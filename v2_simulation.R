library("boot")
library("tidyverse")

# Define the required functions ----

GetMatrix <- function(n, p, c, d.distr = 1) {
  # Creates an n*p data matrix x. x can have different structures. And the 
  # top eigenvalue of its sample covariance is lambda1 which is greater than 1  
  # while the rest eigenvalues are all 1.
  # 
  # Args:
  #   n: Number of rows of the created matrix x.
  #   p: Number of columns of the created matrix x.
  #   c: Determines the top eigenvalue of the sample covariance of x (lambda1)
  #      together with n and p. The relationship is lambda1 = 1 + c*sqrt(p/n).
  #   d.distr: An argument can be chosen from "1", "2" and "3" which represent 
  #      different distributions. "d" is a diagonal matrix which is multiplied  
  #      by other matrices to create the data matrix x. The parameter "d.distr"  
  #      is to specify the distribution of d's diagonal entries denoted by 
  #      "d.ii".
  #      If 1, d.ii will be i.i.d N(0, 1);
  #      if 2, d.ii will be i.i.d uniformly distributed from 1/2 to 
  #      (sqrt(3)*sqrt(4-1/4))/2 - 1/4 (around 1.427);
  #      if 3, d.ii will be i.i.d Exp(sqrt(2)).
  #      Default is "1".
  #      
  # 
  # Returns:
  #   An n*p data matrix x. x is the product of the following matrices in order:
  #     d: An n*n diagnial matrix. The distribution of its diagonal entries is  
  #        determined by the parameter "d.distr";
  #     z0: An n*p matrix i.i.d N(0, 1);
  #     v: The right eigenvectors of the SVD of an n*p matrix with enties i.i.d 
  #        N(0, 1);
  #     lambda: A diagonal matrix of eigenvalues lambda1 = 1 + c*sqrt(p/n) and 
  #             lambda2 = lambda3  = ... = lambdap = 1;
  #     t(v): The transpose of v.
  
  if (d.distr == 1) {
    d.ii <- rnorm(n)
  } else if (d.distr == 2) {
    d.ii <- runif(n, min = 1/2, max = (sqrt(3)*sqrt(4-1/4))/2 - 1/4)
  } else if (d.distr == 3) {
    d.ii <- rexp(n, sqrt(2))
  } else {
    stop("Please choose from '1', '2' and '3'.")
  }
  d <- diag(d.ii)
  z0 <- matrix(rnorm(n*p), nrow = n)
  v <- svd(matrix(rnorm(n*p), nrow = n))$v
  lambda1 <- 1 + c*sqrt(p/n)
  lambda <- diag(c(lambda1, rep(1, p-1)))
  x <- d %*% z0 %*% v %*% lambda %*% t(v)
  
  return(x)
}

GetTopEigenvalue <- function(x, indices) {
  # Get the top eigenvalue of a matrix's sample covariance matrix.
  #
  # Args:
  #   x: The data matrix whose sample covariance matrix's top eigenvalue is to
  #      to be calculated.
  #   indices: The argument to make sure that function "boot::boot()" can run.
  # 
  # Returns:
  #   The top eigenvalue of x's sample covariance matrix.
  
  data <- x[indices, ]
  lambda1.hat <- max(eigen(cov(data))$values)
  
  return(lambda1.hat)
}

PerfSimulation <- function(N, B, n, p, c, d.distr = 1) {
  # Perform the simulation for N times.
  # For each simulation, create an n*p matrix with function "GetMatrix()" and 
  # get the bootstrap estimate of the top eigenvalue of its sample covariance 
  # matrix as well as the bootstrap estimate of its bia and standard error.The 
  # number of the bootstrap samples for each simulation is B.
  #
  # Args:
  #   N: Number of simulation repetitions.
  #   B: Number of the bootstrap samples for each simulation.
  #   n, p, c, d.distr: The arguments from the function "GetMatix()".
  #
  # Returns:
  #   A data frame of:
  #     $lambda1.hat: The top eigenvalues of sample covariance matrices
  #                   (i.e. the statistic of interest) for these N simulations.
  #     $bias.boot: The bootstrap estimate of bias of the top eigenvalues
  #                 for these N simulations.
  #     $std.err.boot: The bootstrap estimate of standard errors of the 
  #                    top eigenvalues for these N simulations.
  #     $c: One of the arguments from "GetMatrix()". lambda1 = 1 + c*sqrt(p/n).
  #     $r: p/n (i.e. number of rows of the original data matrix divided by 
  #         number of columns).
  #     $true.lambda1: The true top eigenvalue.
  #     $true.bias: The "true" bias.
  #     $true.std.err: The "true" standard error.
  
  lambda1.hat <- NA
  bias.boot <- NA
  std.err.boot <- NA
  
  for (i in 1:N) {
    x.matrix <- GetMatrix(n, p, c, d.distr)
    boot.result <- boot(data = x.matrix, statistic = GetTopEigenvalue, R = B)
    lambda1.hat[i] <- boot.result$t0
    lambda1.hat.boot <- boot.result$t
    bias.boot[i] <- mean(lambda1.hat.boot) - lambda1.hat[i]
    std.err.boot[i] <- sqrt(
      sum((lambda1.hat.boot - mean(lambda1.hat.boot[i]))^2) / (B - 1)
    )
  }
  
  r <- p/n
  true.lambda1 <- 1 + c*sqrt(r)
  true.bias <- mean(lambda1.hat) - true.lambda1
  true.std.err <- sd(lambda1.hat)
  
  outcome <- data.frame(lambda1.hat = lambda1.hat,
                        bias.boot = bias.boot,
                        std.err.boot = std.err.boot,
                        c = c,
                        r = r,
                        true.lambda1 = true.lambda1,
                        true.bias = true.bias,
                        true.std.err = true.std.err)
  
  return(outcome)
}

# Simulations for different values of c and r when N = 500, B = 200, n = 300 ----
simu.result <- NA
for (c in c(0, 3, 11, 50, 100)) {
  for (p in c(3, 30, 100, 150)) {
      simu.result <- rbind(simu.result,
                           PerfSimulation(N = 200, B = 100, n = 300,
                                          p = p, c = c, d.distr = 1))
  }
}
simu.result <- simu.result[-1, ]

# plot the bias
simu.data.plot <- simu.result %>%
  mutate(r = as.character(r)) %>%
  mutate(c = as.character(c)) %>%
  mutate(bias.scaled = bias.boot/true.lambda1)
ggplot(simu.data.plot, aes(x = r, y = bias.scaled)) +
  geom_boxplot(aes(fill = c)) +
  stat_summary(fun = "median", geom = "line", aes(group = c), size = .5)
