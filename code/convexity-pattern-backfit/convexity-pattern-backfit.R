################################################################################
# 
# Convexity Pattern Regression: Backfitting Version
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
# Version 1.0: July 29, 2014
# [Run backfit.example1() for a demo (a plot per iteration), and 
#  backfit.example2() and backfit.example3() for more sample results.]
# - A working version of the convexity pattern regression
#   with a backfitting algorithm. Each 1-D regression is 
#   basically a mixed integer conic quadratic program
#   involving two 0-1 integer variables (one for convex, one for concave). 
# - Sparsity is induced by including the equivalent of the L0-norm,
#   as done in the full version. That is, each 1-D regression has 
#   the objective MSE + lambda*(z+w), where z and w are the two integer
#   variables for convex and concave components, respectively.
# - Performance: the algorithm seems to work, but with less stability
#   and with slower convergence rate when compared to the full MISOCP
#   version (empirically). With larger p, starting points seem to 
#   matter a lot, as there are cases of runs in which the size of update 
#   simply blows up.
# 
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

convexity.pattern.regression.1d = function (x, y, B = 10, lambda = 1) {
	
  ########################################
  #
  # 1-D convexity pattern regression function.
  #
  # [Inputs]
  # x: n-vector (covariates)
  # y: n-vector (SCALED outcomes)
  # B: positive scalar (smoothness parameter)
  # lambda: positive scalar (sparsity parameter)
  #
  # [Output]
  # a list consisting of...
  #   $status: program and solution statuses
  #            as given by Rmosek
  #   $pattern: a Boolean pair indicating
  #             (convexity, concavity)
  #             either (0,0), (0,1), or (1,0)
  #   $fit: an n-vector from the component
  #   $f: an n-vector from the convex component
  #   $g: an n-vector from the concave component
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
  n = length(x) # number of points

  # Sort the points and keep the order
  ord = order(x)
  x = x[ord]
  y = y[ord]

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
  # Number of variables: 7n + 1
  # "Effective" number of variables: 2n + 3
  # Number of integer variables: 2
  # For i = 1, ..., n and j = 1, ..., p,
  # fitted values: f_i, g_i 
  # scaled fitted values (auxiliary): u_i, v_i 
  # pointwise errors (auxiliary): r_i
  # slopes for convexity (auxiliary): beta_i, gamma_i (i up to n-1)
  # 0/1 integer variables: z, w
  # objective replacement: t
  f.index = seq(1, n)
  g.index = last(f.index) + seq(1, n) 
  u.index = last(g.index) + seq(1, n)
  v.index = last(u.index) + seq(1, n)
  r.index = last(v.index) + seq(1, n) 
  beta.index = last(r.index) + seq(1, n-1)
  gamma.index = last(beta.index) + seq(1, n-1)
  z.index = last(gamma.index) + 1
  w.index = last(z.index) + 1
  t.index = last(w.index) + 1 
  num.vars = t.index

  # Set up the program
  convexity.pattern = list(sense = "min")

  # Objective: t + lambda * (z + w)
  # where t = sqrt(sum((y_i - (f_i + g_i))^2))
  # Note that MSE = (1/n) * (t^2).
  convexity.pattern$c = c(rep(0, last(gamma.index)), lambda, lambda, 1)

  # Affine constraint 1: auxiliary variables [no cost]
  # r_i = y_i - (f_i + g_i), that is
  # f_i + g_i + r_i = y_i
  
  # number of rows (affine contraints)
  n1 = n
  
  A1 = Matrix(0, nrow = n1, ncol = num.vars)
  A1[, f.index] = .sparseDiagonal(n1)
  A1[, g.index] = .sparseDiagonal(n1)
  A1[, r.index] = .sparseDiagonal(n1)

  # Affine constraint 2: convexity [2*(n-1) + 2*(n-2) affine inequalities]
  # For sorted points x_i, 
  # part 1: (i = 1, ..., n-1)
  # f_{i+1} - f_i = beta_i * (x_{i+1} - x_i)
  # g_{i+1} - g_i = gamma_i * (x_{i+1} - x_i)
  # part 2: (i = 1, ..., n-2)
  # beta_i <= beta_{i+1}
  # gamma_i >= gamma_{i+1}
  
  # number of rows (affine contraints)
  n2.part1 = 2*(n-1)
  n2.part2 = 2*(n-2)
  n2 = n2.part1 + n2.part2
  
  A2 = Matrix(0, nrow = n2, ncol = num.vars)
  	
  # part 1
  diag.f = bandSparse(n-1, n, c(0, 1), list(rep(-1, n-1), rep(1, n-1)))
  diag.g = diag.f
  diag.beta = .sparseDiagonal(n-1, x[2:n] - x[1:(n-1)])
  diag.gamma = diag.beta
  	
  rows.f = seq(1, n-1)
  rows.g = (n-1) + seq(1, n-1)

  A2[rows.f, f.index] = diag.f
  A2[rows.f, beta.index] = -diag.beta
  A2[rows.g, g.index] = diag.g
  A2[rows.g, gamma.index] = -diag.gamma
  	
  # part 2
  diag2.beta = bandSparse(n-2, n-1, c(0, 1), list(rep(-1, n-2), rep(1, n-2)))
  diag2.gamma = -diag2.beta
  
  rows.beta = 2*(n-1) + seq(1, n-2)
  rows.gamma = 2*(n-1) + (n-2) + seq(1, n-2)
  A2[rows.beta, beta.index] = diag2.beta
  A2[rows.gamma, gamma.index] = diag2.gamma

  # Affine constraint 3: identifiability via centering [2 equalities]

  # number of rows (affine contraints)
  n3 = 2

  A3 = Matrix(0, nrow = n3, ncol = num.vars)
  A3[1, f.index] = 1
  A3[2, g.index] = 1

  # Affine constraint 4: choice of pattern [1 integer inequality]
  # 0 <= z + w <= 1, where z, w are integers

  # number of rows (affine contraints)
  n4 = 1

  A4 = Matrix(0, nrow = n4, ncol = num.vars)
  A4[, c(z.index, w.index)] = 1
  
  # Affine constraint 5: scaling the fits f_i and g_i
  #                      into the standard quadratic cone [no cost]
  # f_i = sqrt(n)*B*u_i
  # g_i = sqrt(n)*B*v_i

  # number of rows (affine contraints)
  n5 = 2*n  

  A5 = Matrix(0, nrow = n5, ncol = num.vars)
  A5[, 1:n5] = .sparseDiagonal(n5)
  A5[, (n5 + 1:n5)] = -sqrt(n)*B*.sparseDiagonal(n5)

  # Affine constraints combined with bounds
  convexity.pattern$A = rBind(A1, A2, A3, A4, A5)
  convexity.pattern$bc = rbind(
    blc = c(y, 
            rep(0, n2.part1), rep(0, n2.part2), 
            rep(0, n3), 
            rep(0, n4), 
            rep(0, n5)),
    buc = c(y, 
            rep(0, n2.part1), rep(Inf, n2.part2), 
            rep(0, n3), 
            rep(1, n4),
            rep(0, n5))
  )

  # Constraints on the program variables
  # Integers are either 0 or 1, other constraints are vacuous
  convexity.pattern$bx = rbind(
    blx = c(rep(-Inf, length(c(f.index, g.index, u.index, v.index))),
            rep(-Inf, length(r.index)),
            rep(-Inf, length(c(beta.index, gamma.index))),
            rep(0, length(c(z.index, w.index))), 
            0), # t
    bux = c(rep(Inf, length(c(f.index, g.index, u.index, v.index))),
            rep(Inf, length(r.index)),
            rep(Inf, length(c(beta.index, gamma.index))),
            rep(1, length(c(z.index, w.index))), 
            Inf) # t
  )

  # Conic constraint 1: convexity pattern
  cone.f = list("QUAD", c(z.index, u.index)
  )
  cone.g = list("QUAD", c(w.index, v.index)
  )
  # Conic constraint 2: quadratic objective (mean squared error)
  cone.obj = cbind(list("QUAD", c(t.index, r.index)))

  # Combined
  convexity.pattern$cones = cbind(cone.f, cone.g, cone.obj)

  # Specify integer variables: z and w
  convexity.pattern$intsub = c(z.index, w.index)

  # Solve the program using Rmosek!
  r = mosek(convexity.pattern, opts = list(verbose = -1))

  #---------------------------------------
  # Results
  #---------------------------------------

  # Outputs
  status = list(
    solution = r$sol$int$solsta, 
    program = r$sol$int$prosta
  )
  pattern = matrix(r$sol$int$xx[c(z.index, w.index)])
  colnames(pattern) = c("pattern")
  rownames(pattern) = c("convex", "concave")
  
  f = r$sol$int$xx[f.index][invPerm(ord)]
  g = r$sol$int$xx[g.index][invPerm(ord)]
  MSE = (1/n) * (r$sol$int$xx[t.index]^2)

  print(status)
  print(list(pattern = pattern))

  fit = rep(0, n)
  if (pattern["convex",]) {
  	fit = f
  } else if (pattern["concave",]) {
  	fit = g
  }
  
  return (list(status = status,
               pattern = pattern,
               fit = fit,
               f = f,
               g = g,
               MSE = MSE,
               r = r))	
}

example.1d = function (n = 100, sigma = 1, B = 10, lambda = 1) {

  ########################################
  # Testing the univariate convexity
  # pattern regression function.
  ########################################

  # Generate Synthetic Data
  # True convex function: f(x) = x^4 + 2x + 3
  f = function (x) { x^4 + 2*x + 3 }
  x = runif(n, -5, 5)
  epsilon = rnorm(n, 0, sigma) # Gaussian noise
  f.val = scale(f(x))
  y = f.val + epsilon
  
  result = convexity.pattern.regression.1d(x, y, B, lambda)
  
  par(mfrow = c(1,1))
  plot(x, y, main = "1-D Convexity Pattern Regression")
  ord = order(x)
  lines(x[ord], f.val[ord], lwd = 2, col = "darkgray")
  lines(x[ord], result$f[ord], lwd = 2, col = "red")
  lines(x[ord], result$g[ord], lwd = 2, col = "blue")
  points(x, result$fit, pch = 20, col = "purple")
  mtext(paste(result$status$solution, ", ", result$status$program, sep=""), 
        side = 3, adj = 0)
  mtext(paste("N: ", n, ", MSE: ", result$MSE, sep=""), side = 3, adj = 1)
}

convexity.pattern.backfit = function (X, y, B = 10, lambda = 1, 
                                       max.step = 30, tol = 1e-5) {
	
  ########################################
  #
  # The backfitting function using
  # 1-D convexity pattern regression.
  #
  # [Inputs]
  # X: n by p design matrix (covariates)
  # y: n-vector (outcomes, centered)
  # B: positive scalar (smoothness parameter)
  # max.step: positive integer (# iterations)
  # tol: positive scalar (tolerance for convergence)
  #
  # [Output]
  # a list consisting of...
  #   $pattern: a Boolean matrix indicating
  #             whether each component is
  #             convex, concave, or zero
  #   $fit: an n by p matrix, each column
  #         is the vector of fitted values 
  #         for each component
  #   $MSE: mean squared error
  #   $iter: number of total iterations
  #
  ########################################
  
  #---------------------------------------
  # Setup
  #---------------------------------------

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }

  # Variables
  n = nrow(X)
  p = ncol(X)
  new.fit = Matrix(0, n, p)
  new.pattern = matrix(0, 2, p)
  rownames(new.pattern) = c("convex", "concave")
  colnames(new.pattern) = 1:p
  new.MSE = 0
  
  # Set a random visiting order
  visit.order = sample(p, p, replace = FALSE)
  
  history = matrix(0, 2, max.step)
  rownames(history) = c("L2.distance", "MSE")
  colnames(history) = 1:max.step
  
  #---------------------------------------
  # Iteration
  #---------------------------------------
  
  for (iter in 1:max.step) {
  	
  	current.fit = new.fit
  	current.pattern = new.pattern
  	current.MSE = new.MSE
  	
  	for (j in visit.order) {
  	  
  	  # Fit the residuals
  	  if (p == 2) { # rowSums does not work...
  	  	residual = y - current.fit[,-j]
  	  } else {
  	    residual = y - rowSums(current.fit[,-j])
  	  }
  	  
  	  cat(sprintf("Iteration %d, Component %d:\n", iter, j))
  	  # Note: The resulting fit is centered because of
  	  #       affine constraint 3 in the 1-D program
      result = convexity.pattern.regression.1d(X[,j], residual, B, lambda)
  	  
  	  # Update
  	  new.fit[,j] = result$fit
  	  new.pattern[,j] = result$pattern
  	  
  	}
  	
  	
  	L2.distance = sqrt(sum((new.fit-current.fit)^2))
  	new.MSE = mean((y - rowSums(new.fit))^2)
  	history[,iter] = c(L2.distance, new.MSE)
  	
  	# stoping criterion: L2-distance between current and new fitted values 
  	#                    or change in mean squared error is small
  	if (current.MSE < new.MSE) {
  		cat("Warning: the objective (mean squared error) increased.")
  	}
  	else if (L2.distance <= tol || current.MSE - new.MSE <= tol) {
  		break
  	}
  	cat(sprintf("Iteration %d: L2-distance=%f\n", iter, L2.distance))
  	cat(sprintf("Mean Squared Error: %f\n", new.MSE))
  	print(list(history=history[,1:iter]))
  	cat("========================================\n")
  	
  	if (iter == max.step) {
  		cat("WARNING: Maximum steps reached before convergence.\n")
  	}
  	
  }
  
  return (list(fit = new.fit,
               pattern = new.pattern,
               MSE = mean((y - rowSums(new.fit))^2),
               iter = iter))

}

backfit.example1 = function (n = 100, sigma = 1, B = 10, lambda = 1, 
                          max.step = 30, tol = 1e-5) {
	
  ########################################
  #
  # A demo showing each iteration as a plot.
  # Data: The same as example1 in
  #       convexity.pattern.R. 
  #       p = 2; 1 convex and 1 concave.
  #
  ########################################
  
  #---------------------------------------
  # Setup
  #---------------------------------------

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }
  
  # Fix dimension = 2
  p = 2

  # True convex component: f(x) = 0.1*x^4 + 2x
  # True concave component: g(x) = -4*x^2 - x
  # Function values for each component are then scaled.
  f = function (x) { 0.1*x^4 + 2*x }
  g = function (x) { -4*x^2 - x }

  # X: n by p design matrix, y: scaled function values (n-vector)
  x1 = runif(n, -5, 5)
  x2 = runif(n, -5, 5)
  X = Matrix(c(x1, x2), nrow=n) 

  y1 = scale(sapply(x1, f))
  y2 = scale(sapply(x2, g))

  epsilon = rnorm(n, 0, sigma) # Gaussian noise
  y = y1 + y2 + epsilon

  # Variables
  new.fit = Matrix(0, n, p)
  new.pattern = matrix(0, 2, p)
  rownames(new.pattern) = c("convex", "concave")
  colnames(new.pattern) = 1:p
  new.MSE = 0
  
  # Set a random visiting order
  visit.order = sample(p, p, replace = FALSE)
  
  #---------------------------------------
  # Iteration
  #---------------------------------------
  
  for (iter in 1:max.step) {
  	
  	current.fit = new.fit
  	current.pattern = new.pattern
  	current.MSE = new.MSE
  	
  	for (j in visit.order) {
  	  
  	  # Fit the residuals
  	  if (p == 2) { # rowSums does not work...
  	  	residual = y - current.fit[,-j]
  	  } else {
  	    residual = y - rowSums(current.fit[,-j])
  	  }
  	  
  	  cat(sprintf("Iteration %d, Component %d:\n", iter, j))
  	  # Note: The resulting fit is centered because of
  	  #       affine constraint 3 in the 1-D program
      result = convexity.pattern.regression.1d(X[,j], residual, B, lambda)
  	  
  	  # Update
  	  new.fit[,j] = result$fit
  	  new.pattern[,j] = result$pattern
  	  
  	}
  	
  	L2.distance = sqrt(sum((new.fit-current.fit)^2))
  	new.MSE = mean((y - rowSums(new.fit))^2)
  	
  	par(mfrow = c(1,2))
  	
  	plot(X[,1], y, main = sprintf("Iteration %d: L2-distance = %f", 
  	                              iter, L2.distance), xlab = expression(x[1]))
  	points(X[,1], rowSums(new.fit), pch = 20, col = "purple")
  	ord1 = order(X[,1])
  	lines(X[,1][ord1], y1[ord1], lwd = 2, col = "darkgray")
  	lines(X[,1][ord1], new.fit[,1][ord1], lwd = 2, col = "red")
  	
  	plot(X[,2], y, main = sprintf("MSE = %f", mean((y - rowSums(new.fit))^2)),
  	     xlab = expression(x[2]))
  	points(X[,2], rowSums(new.fit), pch = 20, col = "purple")
  	ord2 = order(X[,2])
  	lines(X[,2][ord2], y2[ord2], lwd = 2, col = "darkgray")
  	lines(X[,2][ord2], new.fit[,2][ord2], lwd = 2, col = "blue")
  	
  	Sys.sleep(0.3)
  	
  	# stoping criterion: L2-distance between current and new fitted values 
  	#                    or change in mean squared error is small
  	if (current.MSE < new.MSE) {
  		cat("Warning: the objective (mean squared error) increased.")
  	}
  	else if (L2.distance <= tol || current.MSE - new.MSE <= tol) {
  		break
  	}
  	cat(sprintf("Iteration %d: L2-distance=%f\n", iter, L2.distance))
  	cat(sprintf("Mean Squared Error: %f\n", new.MSE))
  	cat("========================================\n")
  	
  	if (iter == max.step) {
  		cat("WARNING: Maximum steps reached before convergence.\n")
  	}
  	
  }
  
  return (list(fit = new.fit,
               pattern = new.pattern,
               MSE = mean((y - rowSums(new.fit))^2)))

}

backfit.example2 = function (n = 100, sigma = 1, B = 10, lambda = 1, 
                             max.step = 20, tol = 1e-5) {
  
  ########################################
  # A sample run with the same example2
  # as in convexity-pattern.R.
  ########################################

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }

  # Fix dimension = 5
  p = 5

  # True components (domain [-5, 5]):
  # Note: Function values for each component are then scaled.
  f = list()
  f[[1]] = function (x) { -2*x^2 - 10*x } # concave
  f[[2]] = function (x) { 0.0005*x^8 - 0.0028*x^5 + 0.014*x^2 } # convex
  f[[3]] = function (x) { 5 * (x-20) * log(0.1*(x+6.1)) } # convex
  f[[4]] = function (x) { 0 } # zero
  f[[5]] = function (x) { -60 * exp(0.05*x^2 - 0.1*x) } # concave

  # X: n by p design matrix, y: scaled function values (n-vector)
  X = Matrix(runif(n*p, -5, 5), nrow=n) 
  y.components = sapply(1:p, function (j) { 
  	vals = sapply(X[,j], f[[j]])
  	# scale unless all entries are zero
  	if (sum(abs(vals)) > 0) {
  		return (scale(vals))
  	} else {
  		return (vals)
  	}
  })
  y = rowSums(y.components) + rnorm(n, 0, sigma) # Gaussian noise
  
  # Run the main function
  result = convexity.pattern.backfit(X, y, B, lambda, max.step, tol)
  fit = rowSums(result$fit)

  # Sample plot in each component
  m = matrix(c(1,2,3,4,5,6,6,6,6,6), nrow = 2, ncol = 5, byrow = TRUE)
  layout(mat = m, heights = c(0.4, 0.2))
  for (j in 1:p) {
    par(mar = c(2,4,4,2))
    # sort by x_j for drawing lines
    ord = order(X[,j])
    plot(X[,j], y, 
         main = sprintf("(%d)", j), 
         xlab = expression(x[j]), ylab = expression(y))
    f_j = sapply(X[,j], f[[j]])
    points(X[,j], fit, pch = 20, col = "purple")
    lines(X[,j][ord], y.components[,j][ord], lwd = 2, col = "darkgray")
    if (result$pattern["convex",j]) {
      lines(X[,j][ord], result$fit[,j][ord], lwd = 2, col = "red")    	
    } else if (result$pattern["concave",j]) {
      lines(X[,j][ord], result$fit[,j][ord], lwd = 2, col = "blue")    	
    }
  }

  # Legend & Text
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  mtext(sprintf("#iterations: %d, N: %d, Scaled MSE: %f", 
                result$iter, n, result$MSE), side = 3, adj = 1)
  legend("top", inset = 0, 
         c("Data", "True component", "Convex component", 
           "Concave component", "Additive fit"), 
         col = c("black", "darkgray", "red", "blue", "purple"), 
         pch = c(21, 20, 20, 20, 20))
}

backfit.example3 = function (n = 100, sigma = 1, B = 10, lambda = 1, 
                             max.step = 10, tol = 1e-5) {
  
  ########################################
  # A sample run with the same example3
  # as in convexity-pattern.R.
  ########################################

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }
  
  if (!require("MASS")) {
  	stop ("MASS not installed.")
  }
  
  # Fix dimension = 50
  p = 50
  
  # X: n by p design matrix, y: function values (n-vector)
  
  # The covariates are sampled from a p-dimensional correlated Gaussian 
  # with mean 0, some random covariance.
  # the last two dimensions are zero.
  Sigma.half = Matrix(runif(p*p, -2, 2), nrow = p)
  Sigma = t(Sigma.half) %*% Sigma.half
  X = cBind(mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma))
  # sample again if any value is too large or too small
  while (max(abs(X)) > 50) {
  	cat("Re-sampling covariates from a correlated Gaussian...\n")
  	X = cBind(mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma))
  }
  cat("Obtained a sample of covariates from a correlated Gaussian.\n")

  # True components:
  # Note: Function values for each component are then scaled.
  f = list()
  f[[1]] = function (x) { -0.1*x^4 - 2*x } # concave
  f[[2]] = function (x) { 3*x^2 + x } # convex
  for (j in 3:p) {
  	f[[j]] = function (x) { 0 } # zero
  }

  # y: scaled function values (n-vector)
  y.components = sapply(1:p, function (j) { 
  	vals = sapply(X[,j], f[[j]])
  	# scale unless all entries are zero
  	if (sum(abs(vals)) > 0) {
  		return (scale(vals))
  	} else {
  		return (vals)
  	}
  })
  y = rowSums(y.components) + rnorm(n, 0, sigma) # Gaussian noise
  
  # Run the main function
  result = convexity.pattern.backfit(X, y, B, lambda, max.step, tol)
  fit = rowSums(result$fit)

  # Sample plot in each component
  p.display = 5
  m = matrix(c(1:p.display, rep(p.display+1, p.display)), 
              nrow = 2, ncol = p.display, byrow = TRUE)
  layout(mat = m, heights = c(0.4, 0.25))
  for (j in 1:p.display) {
    par(mar = c(2,4,4,2))
    # sort by x_j for drawing lines
    ord = order(X[,j])
    plot(X[,j], y, 
         main = sprintf("(%d)", j), 
         xlab = expression(x[j]), ylab = expression(y))
    f_j = sapply(X[,j], f[[j]])
    points(X[,j], fit, pch = 20, col = "purple")
    lines(X[,j][ord], y.components[,j][ord], lwd = 2, col = "darkgray")
    if (result$pattern["convex",j]) {
      lines(X[,j][ord], result$fit[,j][ord], lwd = 2, col = "red")    	
    } else if (result$pattern["concave",j]) {
      lines(X[,j][ord], result$fit[,j][ord], lwd = 2, col = "blue")    	
    }
  }

  # Legend & Text
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  mtext(sprintf("#Components: %d, #iterations: %d, N: %d, Scaled MSE: %f", 
                p, result$iter, n, result$MSE), side = 3, adj = 1)
  legend("top", inset = 0, 
         c("Data", "True component", "Convex component", 
           "Concave component", "Additive fit"), 
         col = c("black", "darkgray", "red", "blue", "purple"), 
         pch = c(21, 20, 20, 20, 20))
}