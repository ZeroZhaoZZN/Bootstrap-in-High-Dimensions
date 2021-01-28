library("boot")

GetMatrix <- function(n, p, c, d.distr) {
  # Creates an n*p data matrix x. x can have different structures. And the 
  # top eigenvalue of its sample covariance is lambda1 which is greater than 1  
  # while the rest eigenvalues are all 1.
  # 
  # Args:
  #   n: Number of rows of the created matrix x.
  #   p: Number of columns of the created matrix x.
  #   c: Determines the top eigenvalue of the sample covariance of x (lambda1)
  #      together with n and p. The relationship is lambda1 = 1 + c*sqrt(p/n).
  #   d.distr: A parameter can be choosen from "1", "2" and "3" which represent 
  #      different distributions. "d" is a diagnoal matrix which is multiplied  
  #      by other matrices to create the data matrix x. The parameter "d.distr"  
  #      is to specify the distribution of d's diagonal entries denoted by "d.ii".
  #      If 1, d.ii will be i.i.d N(0, 1);
  #      if 2, d.ii will be i.i.d uniformly distributed from 1/2 to 
  #      (sqrt(3)*sqrt(4-1/4))/2 - 1/4 (around 1.427);
  #      if 3, d.ii will be i.i.d Exp(sqrt(2)).
  #      Default is "1".
  #      
  # Returns:
  #   The n*p data matrix x. x is the product of the following matrices in order:
  #     d: An n*n diagnial matrix. The distribution of its diagonal entries is  
  #      determined by the parameter "d.distr";
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
GetMatrix(n = 300, p = 3, c = 0, d.distr = 1)
