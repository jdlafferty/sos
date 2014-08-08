################################################################################
# 
# Convex Function Generator
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
# Version 1.0: August 7, 2014
# - A working version of convex/concave/zero additive data generator
#   with (at most) 10 non-zero components.
#
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

cvx.generator = function (n = 100, p = 5, sigma = 1, pattern = NULL,
                          num.convex = 2, num.sparse = 1, design = "uniform") {

  ########################################
  #
  # The main function.
  #
  # [Inputs]
  # n: sample size
  # p: dimension i.e. number of components
  #    this argument is ignored if 
  #    'pattern' is specified.
  # sigma: size of the Gaussian noise
  # pattern: pattern as a p-vector
  #          1 convex, -1 concave, 0 sparse
  # num.convex: number of convex components.
  #             this argument is ignored if 
  #             'pattern' is specified.
  # num.sparse: number of sparse components.
  #             this argument is ignored if
  #             'pattern' is specified.
  # design: the design of covariates.
  #         options are "uniform" (default),
  #         "regular", and "gaussian" or "normal".
  #
  # [Output]
  # a list consisting of...
  #   $pattern: a Boolean (0/1) matrix indicating
  #             whether each component is
  #             convex, concave, or zero
  #   $X: an n by p matrix, each column
  #         is the vector of data points 
  #         for each component
  #   $M: an n by p matrix, each column
  #       is a vector of values from
  #       each component
  #   $y: an n-vector of the output values.
  #       note that y = rowSums(M) + noise.
  #
  ########################################

  #---------------------------------------
  # Setup
  #---------------------------------------

  # Require packages
  if (!require("Matrix")) {
    stop ("Matrix package not installed")
  }
  if (!require("MASS")) {
    stop ("MASS package not installed")
  }

  # Make sure that the input pattern is valid, if given
  if (!is.null(pattern)) {
    if (!is.vector(pattern)) {
      stop ("Input 'pattern' must be a p-vector")
    } 
    if (!all(sapply(pattern, 
                    function (i) { (i == 1) || (i == 0) || (i == -1) } ))) {
      stop (paste("Input 'pattern' must contain only",
                "1 (convex), 0 (zero), or -1 (concave)"))
    }
    # Set variables
    p = length(pattern)
    num.convex = length(pattern[pattern==1])
    num.sparse = length(pattern[pattern==0])
  }
  # Choose random pattern if it is not given   
  else {
    pattern = sample(c(rep(1, num.convex), rep(0, num.sparse), 
                     rep(-1, p - num.convex - num.sparse)), p, replace = FALSE)
  }
  
  # make sure that the input design is valid
  design = tolower(design)
  if (!any(c(design == "uniform", design == "regular",
             design == "gaussian", design == "normal"))) {
    stop ("Input 'design' must be either uniform, regular, gaussian, or normal")    
  }

  #---------------------------------------
  # Convex functions on [-1, 1]
  #---------------------------------------

  num.functions = 10
  
  if (p - num.sparse > num.functions) {
    stop (paste("Sorry, the current implementation does not support",
                sprintf("more than %d non-zero components", num.functions)))
  }
  
  f = list()
  
  f[[1]] = function (x) { x^4 - 2*x }
  f[[2]] = function (x) { 2*x^4 + x }
  f[[3]] = function (x) { 3*x^2 + x }
  f[[4]] = function (x) { x^2 + 0.01*x }
  f[[5]] = function (x) { (x-5) * log(3*(x+4)) }
  f[[6]] = function (x) { exp(x^2 - x) }
  f[[7]] = function (x) { exp(x^2 + 0.25*x) }
  f[[8]] = function (x) { x^6 }
  f[[9]] = function (x) { -exp(-x^2 + 0.2*x) }
  f[[10]] = function (x) { -log(x+5) }
  
  # # Show functions
  # x = runif(n, -1, 1)
  # ord = order(x)
  # M = scale(sapply(1:num.functions, function (i) { sapply(x, f[[i]]) }))
  # colors = rainbow(num.functions)
  # plot(x, x, pch = 20)
  # for (i in 1:num.functions) { 
    # lines(x[ord], M[,i][ord], col = colors[i], lwd = 2) 
  # }
  # legend("top", as.character(1:num.functions), col = colors, lwd = 2)
  
  #---------------------------------------
  # Generate covariates
  #---------------------------------------
  
  # Uniform (default): Generate n*p points from Unif(-1, 1)
  if (design == "uniform") {
    X = Matrix(runif(n*p, -1, 1), nrow = n)
  } 
  # Regular: uniformly spaced points between -1 and 1,
  #          permuted randomly in each covariate
  else if (design == "regular") {
    perm = seq(-1, 1, length = n)
    X = Matrix(sapply(1:p, function (i) { sample(perm, n, replace = FALSE) }))
  } 
  # Gaussian/normal: Generate n points from a p-variate correlated Gaussian
  #                  where the covariance matrix is randomly generated < 1
  # (note: to avoid function value errors, we truncate the values 
  #        greater than 1 or less than -1)
  else {
    cov.half = Matrix(runif(p*p, -1, 1), nrow = p)
    cov = t(cov.half) %*% cov.half
    # scale the covariance by its F-norm, to keep things [-1, 1]
    cov = cov / norm(cov, type = "F")
    X = mvrnorm(n, mu = rep(0, p), Sigma = cov)
    
    # remove outliers
    X[X < -1] = -1
    X[X > 1] = 1
  }
  
  #---------------------------------------
  # Sample functions and compute M
  #---------------------------------------
  
  components = sample(1:num.functions, p, replace = FALSE)
  M = Matrix(0, nrow = n, ncol = p)
  c = 1
  
  for (i in 1:p) {
    if (pattern[i] != 0) {
      M[, i] = pattern[i] * scale(sapply(X[, i], f[[c]]))
      c = c + 1
    }
  }
  
  #---------------------------------------
  # Compute y
  #---------------------------------------
  
  noise = rnorm(n, 0, sigma)
  y = rowSums(M) + noise
  
  
  return (list(pattern = pattern,
               X = X,
               M = M,
               y = y))
}