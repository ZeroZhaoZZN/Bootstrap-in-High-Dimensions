library("nimble")
rm(list = ls())

n1 <- 100
n2 <- 500
n3 <- 1000

k1 <- 0.01
k2 <- 0.1
k3 <- 0.3
k4 <- 0.5

p12 <- n1*k2     # 100*10, k = 0.1
p13 <- n1*k3     # 100*30, k = 0.3
p23 <- n2*k3     # 500*150, k = 0.3


set.seed(2000)

# Normal----------

# 100*10
n <- n1
p <- p12
x_norm1 <- vector(mode = "list")
for(i in 1:1000) {
  x_norm1[[i]] <- matrix(rnorm(n*p), nrow = n, ncol = p)
}

# 100*30
n <- n1
p <- p13
x_norm2 <- vector(mode = "list")
for(i in 1:1000) {
  x_norm2[[i]] <- matrix(rnorm(n*p), nrow = n, ncol = p)
}

# 500*150
n <- n2
p <- p23
x_norm3 <- vector(mode = "list")
for(i in 1:1000) {
  x_norm3[[i]] <- matrix(rnorm(n*p), nrow = n, ncol = p)
}


# Double Exponential----------

# 100*10
n <- n1
p <- p12
x_dexp1 <- vector(mode = "list")
for(i in 1:1000) {
  x_dexp1[[i]] <- matrix(rdexp(n*p), nrow = n, ncol = p)
}

# 100*30
n <- n1
p <- p13
x_dexp2 <- vector(mode = "list")
for(i in 1:1000) {
  x_dexp2[[i]] <- matrix(rdexp(n*p), nrow = n, ncol = p)
}

# 500*150
n <- n2
p <- p23
x_dexp3 <- vector(mode = "list")
for(i in 1:1000) {
  x_dexp3[[i]] <- matrix(rdexp(n*p), nrow = n, ncol = p)
}


# Elliptical----------

relli_exp <- function(n, p, ...) {
  lambda <- rexp(n)
  z <- matrix(rnorm(n*p), nrow = n, ncol = p)
  return(lambda*z)
}

relli_norm <- function(n, p, ...) {
  lambda <- rnorm(n)
  z <- matrix(rnorm(n*p), nrow = n, ncol = p)
  return(lambda*z)
}

relli_unif <- function(n, p, ...) {
  lambda <- runif(n)
  z <- matrix(rnorm(n*p), nrow = n, ncol = p)
  return(lambda*z)
}

# 100*10
n <- n1
p <- p12
x_elli_exp1 <- vector(mode = "list")
x_elli_norm1 <- vector(mode = "list")
x_elli_unif1 <- vector(mode = "list")

for(i in 1:1000) {
  x_elli_exp1[[i]] <- relli_exp(n, p, rate = sqrt(2))
  x_elli_norm1[[i]] <- relli_norm(n, p)
  x_elli_unif1[[i]] <- relli_unif(n, p, min = 0.5, max = 1.5)
}

# 100*30
n <- n1
p <- p13
x_elli_exp2 <- vector(mode = "list")
x_elli_norm2 <- vector(mode = "list")
x_elli_unif2 <- vector(mode = "list")

for(i in 1:1000) {
  x_elli_exp2[[i]] <- relli_exp(n, p, rate = sqrt(2))
  x_elli_norm2[[i]] <- relli_norm(n, p)
  x_elli_unif2[[i]] <- relli_unif(n, p, min = 0.5, max = 1.5)
}

# 500*150
n <- n2
p <- p23
x_elli_exp3 <- vector(mode = "list")
x_elli_norm3 <- vector(mode = "list")
x_elli_unif3 <- vector(mode = "list")

for(i in 1:1000) {
  x_elli_exp3[[i]] <- relli_exp(n, p, rate = sqrt(2))
  x_elli_norm3[[i]] <- relli_norm(n, p)
  x_elli_unif3[[i]] <- relli_unif(n, p, min = 0.5, max = 1.5)
}