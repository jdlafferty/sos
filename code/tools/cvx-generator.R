################################################################################
# 
# Convex Sequence Generator
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
# Version 1.1: August 13, 2014
# - Added parse.pattern() for converting pattern matrices into vectors.
# Version 1.2: August 18, 2014
# - Now generates random convex sequences instead of those from known functions,
#   in particular can generate unlimited number of convex sequences.
#   'steepest.slope' is a new parameter serving as the absolute bound for
#   the maximum slope.
#
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

cvx.generator = function (n = 100, p = 5, sigma = 1, pattern = NULL,
                          num.convex = 2, num.sparse = 1, 
                          steepest.slope = 10, design = "uniform") {

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
  # steepest.slope: the absolute maximum of
  #                 the slope in any sequence
  # design: the design of covariates.
  #         options are "uniform" (default),
  #         "regular", and "gaussian" or "normal".
  #
  # [Output]
  # a list consisting of...
  #   $pattern: pattern as a p-vector
  #             1 convex, -1 concave, 0 sparse
  #   $X: an n by p matrix, each column
  #         is the vector of data points 
  #         for each component
  #   $M: an n by p matrix, each column
  #       is a vector of true values from
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
  # Convex sequences on [-1, 1]
  #---------------------------------------

  M = Matrix(0, nrow = n, ncol = p)
  
  # n * (n-1) structured matrix that gives the centered values
  # given the differences
  P = function (n) {
    col = function (i) { c(rep(-(n-i), i), rep(i, n-i)) }
    return ((1/n) * Matrix(sapply(1:(n-1), col)))
  }
  
  for (j in 1:p) {
    if (pattern[j] != 0) {
      
      # Find the ordering of the jth component and sort
      ord = order(X[, j])
      x = sort(X[, j])
      
      # Compute the diagonal matrix of differences
      D = diag(x[2:n] - x[1:(n-1)])
      
      # Generate an increasing slope
      s = sample(steepest.slope, 1)
      beta = sort(runif(n-1, min = -s, max = s))
      
      # Find the convex sequence f and return to original index
      f = scale(P(n) %*% D %*% beta)
      M[, j] = pattern[j] * f[invPerm(ord)]
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


# Parser: from pattern matrix ((1,0), (0,0), (0,1)) 
#         to pattern vector (-1, 0, 1)
parse.pattern = function (pattern.matrix) {
  if (any(colSums(pattern.matrix) > 1)) {
    stop (paste("Error: Resulting pattern contains a component with",
                "both components active"))
  }
  return (as.vector(c(1, -1) %*% pattern.matrix))
}