########################################
# 
# Convexity Pattern Regression
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# Version 1.0: July 25, 2014
#
########################################

convexity.pattern.regression <- function (X, y, B = 100) {

  ########################################
  #
  # The main function.
  #
  # [Inputs]
  # X: n by p design matrix (covariates)
  # y: n-vector (outcomes)
  # B: positive scalar (smoothness parameter)
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
  n = dim(X)[1] # number of points
  p = dim(X)[2] # dimension

  #---------------------------------------
  # Mixed Integer (0/1) SOCP
  # using Rmosek CQP
  #---------------------------------------

  # Helper function: Accessing the last element of a sequence
  last <- function (x) {
  	# x is a sequence
  	return (x[length(x)])
  }

  # [Program variables]
  # Number of variables: 4np + 2n + 2p
  # "Effective" number of variables: 2np + n + 2p
  # Number of integer variables: 2p
  # For i = 1, ..., n and j = 1, ..., p,
  # fitted values: f_ij, g_ij 
  # (ordering: f_11, f_12, ..., f_1p, f_21, ...)
  # scaled fitted values (auxiliary): u_ij, v_ij
  # pointwise errors (auxiliary): r_i
  # slopes for convexity: beta_ij, gamma_ij (i up to n-1)
  # 0/1 integer variables: z_j, w_j
  # objective replacement: t
  f.index <- seq(1, n*p) 
  g.index <- last(f.index) + seq(1, n*p) 
  u.index <- last(g.index) + seq(1, n*p) 
  v.index <- last(u.index) + seq(1, n*p) 
  r.index <- last(v.index) + seq(1, n) 
  beta.index <- last(r.index) + seq(1, (n-1)*p) 
  gamma.index <- last(beta.index) + seq(1, (n-1)*p)
  z.index <- last(gamma.index) + seq(1, p) 
  w.index <- last(z.index) + seq(1, p) 
  t.index <- last(w.index) + 1 
  num.vars <- t.index

  # Helpers: Selecting variables
  # (Think of f_ij as the (i,j)th entry of an n by p matrix.)
  f.select.row <- lapply(1:n, function (i) { seq(i*p-p+1, i*p) })
  f.select.col <- lapply(1:p, function (j) { seq(j, n*p, p) })
  g.select.row <- lapply(1:n, function (i) { last(f.index) + seq(i*p-p+1, i*p) })
  g.select.col <- lapply(1:p, function (j) { last(f.index) + seq(j, n*p, p) })
  beta.select.col <- lapply(1:p, function (j) { 
  	                          last(r.index) + seq(j, (n-1)*p, p) })
  gamma.select.col <- lapply(1:p, function (j) { 
  	                          last(beta.index) + seq(j, (n-1)*p, p) })
  # For conic constraint 1:
  u.select.col <- lapply(1:p, function (j) { last(g.index) + seq(j, n*p, p) })
  v.select.col <- lapply(1:p, function (j) { last(u.index) + seq(j, n*p, p) })

  # Set up the program
  convexity.pattern <- list(sense = "min")

  # Objective: t = sqrt(sum(y_i - sum(f_ij+g_ij)))^2))
  # Note that MSE = (1/n) * (t^2).
  convexity.pattern$c <- c(rep(0, t.index-1), 1)

  # Affine constraint 1: auxiliary variables [no cost]
  # r_i = y_i - sum_j(f_ij + g_ij)
  A1 <- Matrix(0, nrow = n, ncol = num.vars)
  for (i in 1:n) {
    A1[i, f.select.row[[i]]] <- 1
    A1[i, g.select.row[[i]]] <- 1
    A1[i, r.index[i]] <- 1
  }

  # Affine constraint 2: convexity [2*(n-1)*p + 2*(n-2)*p affine inequalities]
  # In each dimension j = 1, ..., p, for sorted points x_i, 
  # part 1: (i = 1, ..., n-1)
  # f_{i+1} - f_i >= beta_i * (x_{i+1} - x_i)
  # g_{i+1} - g_i <= gamma_i * (x_{i+1} - x_i)
  # part 2: (i = 1, ..., n-2)
  # beta_i <= beta_{i+1}
  # gamma_i >= gamma_{i+1}
  
  A2 <- Matrix(0, nrow = 2*(n-1)*p + 2*(n-2)*p, ncol = num.vars)
  
  for (j in 1:p) {
  	# take the ordering (for sorting and putting back in order)
  	ord <- order(X[, j])
  	# the following sequence x is sorted
  	x <- X[ord, j]
  	
  	# part 1
  	diag.f <- bandSparse(n-1, n, c(0, 1), list(rep(-1, n-1), rep(1, n-1)))
  	diag.g <- diag.f
  	diag.beta <- cBind(.sparseDiagonal(n-1, x[2:n] - x[1:(n-1)]),
  	                   Matrix(0, n-1, 1))
  	diag.gamma <- diag.beta
  	
  	rows.f <- 2*(n-1)*(j-1) + seq(1, n-1)
    rows.g <- 2*(n-1)*(j-1) + (n-1) + seq(1, n-1)
  	# permute the (upper-)diagonal matrices back in order
  	# (note: we remove the largest beta constraint as we should here.
  	# diag.beta has an extra 0 just for an inverse permutation
  	# we do NOT know which i is the largest for each j, so we just
  	# skip the entry and keep the total number of betas to be n-1
  	# for each j.)
  	inv.ord <- invPerm(ord)[-ord[n]] # remove the last dummy entry
  	A2[rows.f, f.select.col[[j]]] <- diag.f[, invPerm(ord)]
  	A2[rows.f, beta.select.col[[j]]] <- -diag.beta[, inv.ord]
  	A2[rows.g, g.select.col[[j]]] <- diag.g[, invPerm(ord)]
  	A2[rows.g, gamma.select.col[[j]]] <- -diag.gamma[, inv.ord]
  	
  	# part 2
  	diag2.beta <- bandSparse(n-2, n-1, c(0, 1), list(rep(-1, n-2), rep(1, n-2)))
  	diag2.gamma <- -diag2.beta
  	
  	rows.beta <- 2*(n-1)*p + 2*(n-2)*(j-1) + seq(1, n-2)
    rows.gamma <- 2*(n-1)*p + 2*(n-2)*(j-1) + (n-2) + seq(1, n-2)
    A2[rows.beta, beta.select.col[[j]]] <- diag2.beta[, inv.ord]
    A2[rows.gamma, gamma.select.col[[j]]] <- diag2.gamma[, inv.ord]
  }

  # Affine constraint 3: identifiability via centering [2*p equalities]
  A3 <- Matrix(0, nrow = 2*p, ncol = num.vars)
  for (j in 1:p) {
    A3[j, f.select.col[[j]]] <- 1
    A3[p+j, g.select.col[[j]]] <- 1
  }

  # Affine constraint 4: choice of pattern [p integer inequalities]
  # 0 <= z_j + w_j <= 1, 
  # where z_j, w_j are 0/1 integers
  A4 <- cBind(.sparseDiagonal(p), .sparseDiagonal(p))
  # Add zero columns for unused variables
  A4 <- cBind(Matrix(0, nrow = p, ncol = last(gamma.index)), A4,
              Matrix(0, nrow = p, ncol = 1))

  # Affine constraint 5: scaling the fits f_ij and g_ij
  #                      into the standard quadratic cone [no cost]
  # f_ij = sqrt(n)*B*u_ij
  # g_ij = sqrt(n)*B*v_ij
  num.fits <- 2*n*p
  A5 <- Matrix(0, nrow = num.fits, ncol = num.vars)
  A5[, 1:num.fits] <- .sparseDiagonal(num.fits)
  A5[, (num.fits + 1:num.fits)] <- -sqrt(n)*B*.sparseDiagonal(num.fits)

  # Affine constraints combined with bounds
  convexity.pattern$A <- rBind(A1, A2, A3, A4, A5)
  convexity.pattern$bc <- rbind(
    blc = c(y, 
            rep(0, 2*(n-1)*p + 2*(n-2)*p), 
            rep(0, 2*p), 
            rep(0, p), 
            rep(0, 2*n*p)),
    buc = c(y, 
            rep(0, 2*(n-1)*p), rep(Inf, 2*(n-2)*p), 
            rep(0, 2*p), 
            rep(1, p), 
            rep(0, 2*n*p))
  )

  # Constraints on the program variables
  # Integers are either 0 or 1, other constraints are vacuous
  convexity.pattern$bx <- rbind(
    blx = c(rep(-Inf, 4*n*p), # f, g, u, v
            rep(-Inf, n), # r
            rep(-Inf, 2*(n-1)*p), # beta, gamma
            rep(0, 2*p), # z, w
            0), # t
    bux = c(rep(Inf, 4*n*p), # f, g, u, v
            rep(Inf, n), # r
            rep(Inf, 2*(n-1)*p), # beta, gamma
            rep(1, 2*p), # z, w
            Inf) # t
  )

  # Conic constraint 1: convexity pattern+smoothness [2*p cones]
  cone.f <- sapply(1:p, 
    function (j) { list("QUAD", c(z.index[j], u.select.col[[j]])) }
  )
  cone.g <- sapply(1:p, 
    function (j) { list("QUAD", c(w.index[j], v.select.col[[j]])) }
  )
  # Conic constraint 2: quadratic objective (mean squared error)
  cone.obj <- cbind(list("QUAD", c(t.index, r.index)))

  # Combined
  convexity.pattern$cones <- cbind(cone.f, cone.g, cone.obj)

  # Specify integer variables: z_j and w_j
  convexity.pattern$intsub <- c(z.index, w.index)

  # Solve the program using Rmosek!
  r <- mosek(convexity.pattern)

  #---------------------------------------
  # Results
  #---------------------------------------

  # Outputs
  status <- list(
    solution = r$sol$int$solsta, 
    program = r$sol$int$prosta
  )
  pattern <- r$sol$int$xx[c(z.index, w.index)]
  pattern <- sapply(1:p, function(j) { c(pattern[j], pattern[p+j]) })
  colnames(pattern) <- 1:p
  rownames(pattern) <- c("convex", "concave")
  fit <- r$sol$int$xx[c(f.index, g.index)]
  MSE <- (1/n) * (r$sol$int$xx[t.index]^2)

  print(status)
  print(list(pattern = pattern))

  fit.component <- function(j) {
  	
    ########################################
    # Given a component index j and
    # the outputs of the above program,
    # returns the fit of
    # the j-th component.
    ########################################

    if (pattern["convex", j]) {
      # return only f for component j
      return (fit[f.select.col[[j]]])
    } else if (pattern["concave", j]) {
      # return only g for component j
      return (fit[g.select.col[[j]]])
    } else {
      return (rep(0, n))
    }
  }

  return (list(status = status,
               pattern = pattern,
               fit = sapply(1:p, fit.component),
               f = Matrix(r$sol$int$xx[f.index], ncol = p, byrow = TRUE),
               g = Matrix(r$sol$int$xx[g.index], ncol = p, byrow = TRUE),
               MSE = MSE,
               r = r))
}

#---------------------------------------
# Examples
#---------------------------------------

example1 <- function (n = 100, sigma = 10, B = 100) {
	
  ########################################
  # A sample run with outputs and plots.
  # NOTE: p=2.
  ########################################
  
  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }

  # Fix dimension = 2
  p = 2

  # True convex component: f(x) = 0.1*x^4 + 2x
  # True concave component: g(x) = -4*x^2 - x
  # Function values for each component are then centered.
  f <- function (x) { 0.1*x^4 + 2*x }
  g <- function (x) { -4*x^2 - x }

  # X: n by p design matrix, y: function values (n-vector)
  x1 <- runif(n, -5, 5)
  x2 <- runif(n, -5, 5)
  X <- Matrix(c(x1, x2), nrow=n) 

  y1 <- sapply(x1, f)
  y2 <- sapply(x2, g)
  y1 <- y1 - mean(y1)
  y2 <- y2 - mean(y2)

  epsilon <- rnorm(n, 0, sigma) # Gaussian noise
  y <- y1 + y2 + epsilon

  plot(x1, y, 
       main = expression(y == (0.1*x[1]^4 + 2*x[1]) + (-4*x[2]^2 - x[2]) + epsilon), 
       xlab=expression(x[1]))
  plot(x2, y, 
       main = expression(y == (0.1*x[1]^4 + 2*x[1]) + (-4*x[2]^2 - x[2]) + epsilon), 
       xlab=expression(x[2]))
       
  # Run the main function and print out the patterns
  result <- convexity.pattern.regression(X, y, B)
  f1 <- result$f[,1]
  f2 <- result$f[,2]
  g1 <- result$g[,1]
  g2 <- result$g[,2]

  # Sample plot in each component
  m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
  layout(mat = m, heights = c(0.4, 0.2))

  # Component 1
  par(mar = c(2,4,4,2))
  plot(x1, y, 
       main = "Convexity Pattern Regression (1)", xlab = expression(x[1]))
  # sort by order of x1 for drawing lines
  ord1 <- order(x1)
  lines(x1[ord1], f1[ord1], lwd = 2, col = "red") # convex comp
  lines(x1[ord1], g1[ord1], lwd = 2, col = "blue") # concave comp
  mtext(paste(result$status$solution, ", ", result$status$program, sep=""), 
        side = 3, adj = 0)
  mtext(paste("N: ", n, ", MSE: ", result$MSE, sep=""), side = 3, adj = 1)

  # Component 2
  par(mar = c(2,4,4,2))
  plot(x2, y, 
       main = "Convexity Pattern Regression (2)", xlab = expression(x[2]))
  # sort by the order of x2 for drawing lines
  ord2 <- order(x2)
  lines(x2[ord2], f2[ord2], lwd = 2, col = "red") # convex comp
  lines(x2[ord2], g2[ord2], lwd = 2, col = "blue") # concave comp
  mtext(paste(result$status$solution, ", ", result$status$program, sep=""), 
        side = 3, adj = 0)
  mtext(paste("N: ", n, ", MSE: ", result$MSE, sep=""), side = 3, adj = 1)

  # Legend
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", inset = 0, 
         c(expression(y == (0.1*x[1]^4 + 2*x[1]) + (-4*x[2]^2 - x[2]) + epsilon), 
           "Convex component", "Concave component"), 
         col = c("black", "red", "blue"), pch = c(21, 20, 20))

  return (result)
}

example2 <- function (n = 100, sigma = 10, B = 800) {
  
  ########################################
  # Another sample run with more components.
  ########################################

  # Import Rmosek
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }

  # Fix dimension = 5
  p = 5

  # True components (domain [-5, 5]):
  # Note: Function values for each component are then centered.
  f <- list()
  f[[1]] <- function (x) { -2*x^2 - 10*x } # concave
  f[[2]] <- function (x) { 0.0005*x^8 - 0.0028*x^5 + 0.014*x^2 } # convex
  f[[3]] <- function (x) { 5 * (x-20) * log(0.1*(x+6.1)) } # convex
  f[[4]] <- function (x) { 0 } # zero
  f[[5]] <- function (x) { -60 * exp(0.05*x^2 - 0.1*x) } # concave

  # X: n by p design matrix, y: function values (n-vector)
  X <- Matrix(runif(n*p, -5, 5), nrow=n) 

  y <- rowSums(sapply(1:p, function (j) {
  	                    f_j <- sapply(X[,j], f[[j]])
  	                    return (f_j - mean(f_j))
                        })) 
  y <- y + rnorm(n, 0, sigma) # Gaussian noise
         
  # Run the main function and print out the patterns
  result <- convexity.pattern.regression(X, y, B)
  fit <- rowSums(result$fit)

  # Sample plot in each component
  m <- matrix(c(1,2,3,4,5,6,6,6,6,6), nrow = 2, ncol = 5, byrow = TRUE)
  layout(mat = m, heights = c(0.4, 0.2))
  for (j in 1:p) {
    par(mar = c(2,4,4,2))
    # sort by x_j for drawing lines
    ord <- order(X[,j])
    plot(X[,j], y, 
         main = sprintf("Convexity Pattern Regression (%d)", j), 
         xlab = expression(x[j]))
    f_j <- sapply(X[,j], f[[j]])
    lines(X[,j][ord], (f_j - mean(f_j))[ord], lwd = 2, col = "darkgray")
    lines(X[,j][ord], result$f[,j][ord], lwd = 2, col = "red")
    lines(X[,j][ord], result$g[,j][ord], lwd = 2, col = "blue")
    points(X[,j], fit, pch = 20, col = "purple")
  }

  # Legend & Text
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  mtext(paste(result$status$solution, ", ", result$status$program, sep=""), 
          side = 3, adj = 0)
  mtext(paste("N: ", n, ", MSE: ", result$MSE, sep=""), side = 3, adj = 1)
  legend("top", inset = 0, 
         c("Data", "True component", "Convex component", 
           "Concave component", "Additive fit"), 
         col = c("black", "darkgray", "red", "blue", "purple"), 
         pch = c(21, 20, 20, 20, 20))
  
  return (result)
}
