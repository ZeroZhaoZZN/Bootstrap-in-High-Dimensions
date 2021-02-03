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
  #   d.distr: An argument can be chosen from "1", "2", "3" and "4" which 
  #      represent different distributions. "d" is a diagonal matrix which is   
  #      multiplied by other matrices to create the data matrix x. The parameter   
  #      "d.distr" is to specify the distribution of d's diagonal entries  
  #      which is denoted by "d.ii".
  #      If 1, d will be an identity matrix;
  #      If 2, d.ii will be i.i.d N(0, 1);
  #      if 3, d.ii will be i.i.d uniformly distributed from 1/2 to 
  #      (sqrt(3)*sqrt(4-1/4))/2 - 1/4 (around 1.427);
  #      if 4, d.ii will be i.i.d Exp(sqrt(2)).
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
    d.ii <- rep(1, n)
  } else if (d.distr == 2) {
    d.ii <- rnorm(n)
  } else if (d.distr == 3) {
    d.ii <- runif(n, min = 1/2, max = (sqrt(3)*sqrt(4-1/4))/2 - 1/4)
  } else if (d.distr == 4) {
    d.ii <- rexp(n, sqrt(2))
  } else {
    stop("Please choose from '1', '2', '3' and '4'.")
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
  # matrix as well as the bootstrap estimate of its bias and standard error.The 
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
  #     $bias.boot: The bootstrap estimates of bias of the top eigenvalues
  #                 for these N simulations.
  #     $std.err.boot: The bootstrap estimates of standard errors of the 
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


# Simulations for different values of c, r and d.distr ----
# when N = 500, B = 100, n = 300

# Z ~ Norm
simu.norm <- NA  
for (c in c(0, 3, 11, 50, 100)) {
  for (p in c(3, 30, 100, 150)) {
    simu.norm <- rbind(simu.norm,
                       PerfSimulation(N = 500, B = 100, n = 300,
                                      p = p, c = c, d.distr = 1))
  }
}
simu.norm <- simu.norm[-1, ]

# Z ~ Ellip Norm
simu.ellip.norm <- NA  
for (c in c(0, 3, 11, 50, 100)) {
  for (p in c(3, 30, 100, 150)) {
    simu.ellip.norm <- rbind(simu.ellip.norm,
                          PerfSimulation(N = 500, B = 100, n = 300,
                                         p = p, c = c, d.distr = 2))
  }
}
simu.ellip.norm <- simu.ellip.norm[-1, ]

# Z ~ Ellip Uniform
simu.ellip.unif <- NA  
for (c in c(0, 3, 11, 50, 100)) {
  for (p in c(3, 30, 100, 150)) {
    simu.ellip.unif <- rbind(simu.ellip.unif,
                          PerfSimulation(N = 500, B = 100, n = 300,
                                         p = p, c = c, d.distr = 3))
  }
}
simu.ellip.unif <- simu.ellip.unif[-1, ]

# Z ~ Ellip Exp
simu.ellip.exp <- NA  
for (c in c(0, 3, 11, 50, 100)) {
  for (p in c(3, 30, 100, 150)) {
    simu.ellip.exp <- rbind(simu.ellip.exp,
                          PerfSimulation(N = 500, B = 100, n = 300,
                                         p = p, c = c, d.distr = 4))
  }
}
simu.ellip.exp <- simu.ellip.exp[-1, ]


# Plot the bias and variance ----

labels <- c(expression(lambda[1]*"="*1),
            expression(lambda[1]*"="*1 + 3*sqrt(r)),
            expression(lambda[1]*"="*1 + 11*sqrt(r)),
            expression(lambda[1]*"="*1 + 50*sqrt(r)),
            expression(lambda[1]*"="*1 + 100*sqrt(r)))
# Z ~ Norm
norm.data.plot <- simu.norm %>%
  mutate(r = factor(r)) %>%
  mutate(c = factor(c)) %>%
  mutate(bias.ratio = abs(bias.boot)/true.bias) %>%
  mutate(variance.ratio = (std.err.boot/true.std.err)^2)

bias.ratio.lim1 <- boxplot.stats(norm.data.plot$bias.ratio)$stats[c(1, 5)]
variance.ratio.lim1 <- boxplot.stats(norm.data.plot$variance.ratio)$stats[c(1, 5)]

p.norm.bias <- ggplot(norm.data.plot, aes(x = r, y = bias.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = bias.ratio.lim1*6) +
  labs(title = expression(
    frac("Bootstrap estimate of bias", "True bias") * " against r (Z ~ Norm)"))
p.norm.bias

p.norm.var <- ggplot(norm.data.plot, aes(x = r, y = variance.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = variance.ratio.lim1*6) +
  labs(title = expression(
    frac("var(" * lambda[1]^{"*"} * ")", "var(" * widehat(lambda[1]) * ")") *  
      " against r (Z ~ Norm)"))
p.norm.var

# Z ~ Ellip Norm
ellip.norm.data.plot <- simu.ellip.norm %>%
  mutate(r = factor(r)) %>%
  mutate(c = factor(c)) %>%
  mutate(bias.ratio = abs(bias.boot)/true.bias) %>%
  mutate(variance.ratio = (std.err.boot/true.std.err)^2)

bias.ratio.lim2 <- boxplot.stats(ellip.norm.data.plot$bias.ratio)$stats[c(1, 5)]
variance.ratio.lim2 <- boxplot.stats(ellip.norm.data.plot$variance.ratio)$stats[c(1, 5)]

p.ellip.norm.bias <- ggplot(ellip.norm.data.plot, aes(x = r, y = bias.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = bias.ratio.lim2*1.8) +
  labs(title = expression(
    frac("Bootstrap estimate of bias", "True bias") * " against r (Z ~ Norm)"))
p.ellip.norm.bias

p.ellip.norm.var <- ggplot(ellip.norm.data.plot, aes(x = r, y = variance.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = variance.ratio.lim2*4) +
  labs(title = expression(
    frac("var(" * lambda[1]^{"*"} * ")", "var(" * widehat(lambda[1]) * ")") *  
      " against r (Z ~ Ellip Norm)"))
p.ellip.norm.var

# Z ~ Ellip Uniform

ellip.unif.data.plot <- simu.ellip.unif %>%
  mutate(r = factor(r)) %>%
  mutate(c = factor(c)) %>%
  mutate(bias.ratio = abs(bias.boot)/true.bias) %>%
  mutate(variance.ratio = (std.err.boot/true.std.err)^2)

bias.ratio.lim3 <- boxplot.stats(ellip.unif.data.plot$bias.ratio)$stats[c(1, 5)]
variance.ratio.lim3 <- boxplot.stats(ellip.unif.data.plot$variance.ratio)$stats[c(1, 5)]

p.ellip.unif.bias <- ggplot(ellip.unif.data.plot, aes(x = r, y = bias.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = bias.ratio.lim3*4.5) +
  labs(title = expression(
    frac("Bootstrap estimate of bias", "True bias") * " against r (Z ~ Norm)"))
p.ellip.unif.bias

p.ellip.unif.var <- ggplot(ellip.unif.data.plot, aes(x = r, y = variance.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = variance.ratio.lim3*7) +
  labs(title = expression(
    frac("var(" * lambda[1]^{"*"} * ")", "var(" * widehat(lambda[1]) * ")") *  
      " against r (Z ~ Ellip Uniform)"))
p.ellip.unif.var

# Z ~ Ellip Exp

ellip.exp.data.plot <- simu.ellip.exp %>%
  mutate(r = factor(r)) %>%
  mutate(c = factor(c)) %>%
  mutate(bias.ratio = abs(bias.boot)/true.bias) %>%
  mutate(variance.ratio = (std.err.boot/true.std.err)^2)

bias.ratio.lim4 <- boxplot.stats(ellip.exp.data.plot$bias.ratio)$stats[c(1, 5)]
variance.ratio.lim4 <- boxplot.stats(ellip.exp.data.plot$variance.ratio)$stats[c(1, 5)]

p.ellip.exp.bias <- ggplot(ellip.exp.data.plot, aes(x = r, y = bias.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = bias.ratio.lim4*1) +
  labs(title = expression(
    frac("Bootstrap estimate of bias", "True bias") * " against r (Z ~ Norm)"))
p.ellip.exp.bias

p.ellip.exp.var <- ggplot(ellip.exp.data.plot, aes(x = r, y = variance.ratio)) +
  geom_boxplot(aes(fill = c), outlier.size = 1, outlier.shape = NA) +
  stat_summary(fun = "median", 
               geom = "line",
               size = .5,
               aes(group = c, color = c),  
               position = position_dodge(0.75)) +
  scale_x_discrete(labels = c("0.01", "0.1", "0.3", "0.5")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = labels) +
  scale_fill_discrete(labels = labels) +
  coord_cartesian(ylim = variance.ratio.lim4*2) +
  labs(title = expression(
    frac("var(" * lambda[1]^{"*"} * ")", "var(" * widehat(lambda[1]) * ")") *  
    " against r (Z ~ Ellip Exp)"))
p.ellip.exp.var
