################################################################################
# 
# Convexity Pattern Regression as a Convex Program: Backfitting
# (using Sabyasachi's new lasso-like idea)
#
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
#
# 
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

# Source the 1-D version
SOURCE_DIRECTORY = "~/Code/sos-convexity/sos/"
source(paste(SOURCE_DIRECTORY, 
             "code/convexity-pattern-cvx/convexity-pattern-cvx-1d.R", sep = ""))

convexity.pattern.backfit = function (X, y, lambda = 1, 
                                      max.step = 30, tol = 1e-5) {
	
  ########################################
  #
  # The backfitting function using
  # 1-D convexity pattern regression
  # as a convex program.
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
