################################################################################
# 
# Convexity Pattern Regression as Lasso
# (using Sabyasachi's new formulation)
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
# Version 1.0: August 15, 2014
# - A working version of the p-dimensional convexity pattern problem
#   with the lasso-like penalty on the difference in subgradients,
#   summed over the p components.
# - Initial observations: the result is very sensitive to the choice of lambda,
#   and the fit is less accurate (less linear pieces) than the full MISOCP.
#
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

# Change this to the location of your local repository (must end with /)
SOURCE_DIRECTORY = "~/Code/sos-convexity/sos/"

convexity.pattern.lasso = function (X, y, lambda = 1) {

  ########################################
  #
  # The main function.
  #
  # [Inputs]
  # X: n by p design matrix (covariates)
  # y: n-vector (outcomes, centered)
  # lambda: positive scalar (sparsity&smoothness parameter)
  #
  # [Output]
  # a list consisting of...
  #   $status: program and solution statuses
  #            as given by Rmosek
  #   $pattern: a Boolean matrix indicating
  #             whether each component is
  #             convex, concave, or zero
  #   $fit: an n by p matrix, each column
  #         is the vector of fitted values 
  #         for each component
  #   $f: an n by p matrix, values from
  #       the convex component
  #   $g: an n by p matrix, values from
  #       the concave component
  #   $MSE: mean squared error
  #   $r: raw output from Rmosek solver
  #
  ########################################

  #---------------------------------------
  # Setup
  #---------------------------------------

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }
  
  # Program size
  n = nrow(X) # number of points
  p = ncol(X) # dimension

  #---------------------------------------
  # Mixed Integer (0/1) SOCP
  # using Rmosek CQP
  #---------------------------------------

  # Helper function: Accessing the last element of a sequence
  last = function (x) {
    # x is a sequence
    return (x[length(x)])
  }

  # [Program variables]
  # Number of variables: 4np + n - 2p + 2
  # "Effective" number of variables: 2np + 1
  # Number of integer variables: 0 (convex program)
  # For i = 1, ..., n and j = 1, ..., p,
  # fitted values: f_ij, g_ij 
  # (ordering: f_11, f_12, ..., f_1p, f_21, ...)
  # scaled fitted values (auxiliary): u_ij, v_ij
  # pointwise errors (auxiliary): r_i
  # slopes for convexity (auxiliary): beta_ij, gamma_ij (i up to n-1)
  # RQUAD constant 1/2: s
  # objective (MSE) replacement: t
  f.index = seq(1, n*p) 
  g.index = last(f.index) + seq(1, n*p) 
  r.index = last(g.index) + seq(1, n) 
  beta.index = last(r.index) + seq(1, (n-1)*p) 
  gamma.index = last(beta.index) + seq(1, (n-1)*p)
  s.index = last(gamma.index) + 1
  t.index = last(s.index) + 1 
  num.vars = t.index

  # Helpers: Selecting variables
  # (Think of f_ij as the (i,j)th entry of an n by p matrix.)
  f.select.row = lapply(1:n, function (i) { seq(i*p-p+1, i*p) })
  f.select.col = lapply(1:p, function (j) { seq(j, n*p, p) })
  g.select.row = lapply(1:n, function (i) { last(f.index) + seq(i*p-p+1, i*p) })
  g.select.col = lapply(1:p, function (j) { last(f.index) + seq(j, n*p, p) })
  beta.select.col = lapply(1:p, function (j) { 
                              last(r.index) + seq(j, (n-1)*p, p) })
  gamma.select.col = lapply(1:p, function (j) { 
                              last(beta.index) + seq(j, (n-1)*p, p) })

  # Set up the program
  convexity.pattern = list(sense = "min")

  # Objective: (1/n) * t + lambda * penalty
  # where t = sum_i((y_i - sum_j(f_ij+g_ij))^2) and
  # penalty = sum_j (beta_{n-1,j} - beta_{1,j} + gamma_{1,j} - gamma_{n-1,j})
  # and the indices n-1, 1 actually indicate the indices of 
  # the second-largest and the smallest values from {x_1j, ..., x_nj}.
  # Note: MSE = (1/n) * t
  objective = rep(0, num.vars)

  objective[t.index] = 1/n
  
  penalty.beta = c(-lambda, rep(0, n-3), lambda)
  penalty.gamma = c(lambda, rep(0, n-3), -lambda)
  for (j in 1:p) {
  	ord = order(X[, j])
  	inv.ord = invPerm(ord)[-ord[n]]
    objective[beta.select.col[[j]]] = penalty.beta[inv.ord]
    objective[gamma.select.col[[j]]] = penalty.gamma[inv.ord]
  }
  
  convexity.pattern$c = objective

  # Affine constraint 1: auxiliary variables for residuals [no cost]
  # r_i = y_i - sum_j(f_ij + g_ij)
  
  # number of rows (affine contraints)
  n1 = n
  
  A1 = Matrix(0, nrow = n1, ncol = num.vars)
  for (i in 1:n) {
    A1[i, f.select.row[[i]]] = 1
    A1[i, g.select.row[[i]]] = 1
    A1[i, r.index[i]] = 1
  }

  # Affine constraint 2: convexity
  # [2*(n-1) affine no-cost equalities, 2*(n-2) affine inequalities]
  # In each dimension j = 1, ..., p, for SORTED (in j-th component) points x_i, 
  # part 1: (i = 1, ..., n-1)
  # f_{i+1} - f_i = beta_i * (x_{i+1} - x_i)
  # g_{i+1} - g_i = gamma_i * (x_{i+1} - x_i)
  # part 2: (i = 1, ..., n-2)
  # beta_i <= beta_{i+1}
  # gamma_i >= gamma_{i+1}
  
  # number of rows (affine contraints)
  n2.part1 = 2*(n-1)*p  
  n2.part2 = 2*(n-2)*p
  n2 = n2.part1 + n2.part2
  
  A2 = Matrix(0, nrow = n2, ncol = num.vars)
  
  for (j in 1:p) {
    # take the ordering (for sorting and putting back in order)
    ord = order(X[, j])
    # the following sequence x is sorted
    x = X[ord, j]
    
    # part 1
    diag.f = bandSparse(n-1, n, c(0, 1), list(rep(-1, n-1), rep(1, n-1)))
    diag.g = diag.f
    diag.beta = cBind(.sparseDiagonal(n-1, x[2:n] - x[1:(n-1)]),
                       Matrix(0, n-1, 1))
    diag.gamma = diag.beta
    
    rows.f = 2*(n-1)*(j-1) + seq(1, n-1)
    rows.g = 2*(n-1)*(j-1) + (n-1) + seq(1, n-1)
    # permute the (upper-)diagonal matrices back in order
    inv.ord = invPerm(ord)[-ord[n]] # remove the last dummy entry
    A2[rows.f, f.select.col[[j]]] = diag.f[, invPerm(ord)]
    A2[rows.f, beta.select.col[[j]]] = -diag.beta[, inv.ord]
    A2[rows.g, g.select.col[[j]]] = diag.g[, invPerm(ord)]
    A2[rows.g, gamma.select.col[[j]]] = -diag.gamma[, inv.ord]
    
    # part 2
    diag2.beta = bandSparse(n-2, n-1, c(0, 1), list(rep(-1, n-2), rep(1, n-2)))
    diag2.gamma = -diag2.beta
    
    rows.beta = 2*(n-1)*p + 2*(n-2)*(j-1) + seq(1, n-2)
    rows.gamma = 2*(n-1)*p + 2*(n-2)*(j-1) + (n-2) + seq(1, n-2)
    A2[rows.beta, beta.select.col[[j]]] = diag2.beta[, inv.ord]
    A2[rows.gamma, gamma.select.col[[j]]] = diag2.gamma[, inv.ord]
  }

  # Affine constraint 3: identifiability constraints [4*p equalities]
  # We want the fitted values of each component j to be centered
  # and also orthogonal to the data vector x_j. (inner product = 0)
  # This eliminates the subspace of vectors that are both convex and concave.

  # number of rows (affine contraints)
  n3 = 4*p

  A3 = Matrix(0, nrow = n3, ncol = num.vars)
  for (j in 1:p) {
    A3[0*p+j, f.select.col[[j]]] = rep(1, n)
    A3[1*p+j, g.select.col[[j]]] = rep(1, n)
    A3[2*p+j, f.select.col[[j]]] = X[, j]
    A3[3*p+j, g.select.col[[j]]] = X[, j]
  }

  # Affine constraint 4: auxiliary variable for Rmosek rotated cone [no cost]
  # s = 1/2 
  
  # number of rows (affine contraints)
  n4 = 1
  
  A4 = Matrix(0, nrow = n4, ncol = num.vars)
  A4[, s.index] = 1

  # Affine constraints combined with bounds
  convexity.pattern$A = rBind(A1, A2, A3, A4)
  convexity.pattern$bc = rbind(
    blc = c(y, 
            rep(0, n2.part1), rep(0, n2.part2), 
            rep(0, n3), 
            rep(1/2, n4)),
    buc = c(y, 
            rep(0, n2.part1), rep(Inf, n2.part2), 
            rep(0, n3), 
            rep(1/2, n4))
  )

  # Constraints on the program variables
  convexity.pattern$bx = rbind(
    blx = c(rep(-Inf, length(c(f.index, g.index))),
            rep(-Inf, length(r.index)),
            rep(-Inf, length(c(beta.index, gamma.index))),
            0, # s
            0), # t
    bux = c(rep(Inf, length(c(f.index, g.index))),
            rep(Inf, length(r.index)),
            rep(Inf, length(c(beta.index, gamma.index))),
            Inf, # s
            Inf) # t
  )

  # Conic constraint 2: quadratic objective (mean squared error)
  # Note: In Rmosek terms, this is 2*t*s >= sum(r^2)
  #       where s = 1/2. 
  convexity.pattern$cones = cbind(list("RQUAD", c(t.index, s.index, r.index)))

  # Solve the program using Rmosek!
  r = mosek(convexity.pattern)

  #---------------------------------------
  # Results
  #---------------------------------------

  # Outputs
  status = list(
    solution = r$sol$itr$solsta, 
    program = r$sol$itr$prosta
  )
  output = r$sol$itr$xx

  # helper: check if the fit is "essentially" zero
  is.zero = function (fit, tol = 1e-4) {
    return (max(abs(fit)) < n*tol)
  }
  # helper: returns the pattern matrix
  get.pattern = function(f, g) {
    p = matrix(as.integer(c(!is.zero(f), !is.zero(g))), nrow = 2)
    rownames(p) = c("convex", "concave")
    colnames(p) = c("pattern")
    return (p)
  }

  pattern = sapply(1:p, function (j) { 
    get.pattern(output[f.select.col[[j]]], output[g.select.col[[j]]])
  })
  colnames(pattern) = 1:p
  rownames(pattern) = c("convex", "concave")

    
  f = Matrix(output[f.index], ncol = p, byrow = TRUE)
  g = Matrix(output[g.index], ncol = p, byrow = TRUE)
  MSE = (1/n) * output[t.index]

  print(status)
  print(list(pattern = pattern))

  fit.component = function (j) {
    
    ########################################
    # Given a component index j and
    # the outputs of the above program,
    # returns the fit of
    # the j-th component.
    ########################################

    if (pattern["convex", j]) {
      # return only f for component j
      return (f[, j])
    } else if (pattern["concave", j]) {
      # return only g for component j
      return (g[, j])
    } else {
      return (rep(0, n))
    }
  }

  return (list(status = status,
               pattern = pattern,
               fit = sapply(1:p, fit.component),
               f = f,
               g = g,
               MSE = MSE,
               r = r))
}

# Source the example generator (cvx.generator)
source(paste(SOURCE_DIRECTORY, "code/tools/cvx-generator.R", sep = ""))
# Source the plotting function (plot.components)
source(paste(SOURCE_DIRECTORY, "code/tools/plot-components.R", sep = ""))

example = function (n = 300, pattern = c(-1, 1, 1, 0, -1), lambda = 0.1) {
	
	data = cvx.generator(n, pattern = pattern)
	result = convexity.pattern.lasso(data$X, data$y, lambda)
	plot.components(data$X, data$y, result$fit, 
	                parse.pattern(result$pattern), true.components = data$M)	
	
	return (list(pattern = result$pattern, fit = result$fit))
}

