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

# Source the 1-D version (convexity.pattern.cvx.1d.auto)
SOURCE_DIRECTORY = "~/Code/sos-convexity/sos/"
source(paste(SOURCE_DIRECTORY, 
             "code/convexity-pattern-cvx/convexity-pattern-cvx-1d.R", sep = ""))

convexity.pattern.cvx.backfit = function (X, y, step = 0.1, 
                                          max.step = 20, tol = 1e-5) {
  
  ########################################
  #
  # The backfitting function using
  # 1-D convexity pattern regression
  # as a convex program.
  #
  # [Inputs]
  # X: n by p design matrix (covariates)
  # y: n-vector (outcomes, centered)
  # step: positive scalar (accuracy of lambda)
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
  new.MSE = Inf
  
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
      
      # Note: The resulting fit is centered because of
      #       affine constraint 3 in the 1-D program
      result = convexity.pattern.cvx.1d.auto(X[,j], residual, step)
      
      # Update
      new.fit[,j] = result$fit
      new.pattern[,j] = result$pattern
      
    }
    
    
    L2.distance = sqrt(sum((new.fit-current.fit)^2))
    new.MSE = mean((y - rowSums(new.fit))^2)
    
    # stoping criterion: L2-distance between current and new fitted values 
    #                    or change in mean squared error is small
    if (current.MSE < new.MSE) {
      cat("Warning: the objective (mean squared error) increased.\n")
    }
    else if (L2.distance <= tol || current.MSE - new.MSE <= tol) {
      break
    }
    cat(sprintf("Iteration %d: L2-distance=%f\n", iter, L2.distance))
    cat(sprintf("Mean Squared Error: %f\n", new.MSE))
    cat("========================================\n")
    
    if (iter == max.step) {
      cat("Maximum steps reached before convergence.\n")
    }
    
  }
  
  return (list(fit = new.fit,
               pattern = new.pattern,
               MSE = mean((y - rowSums(new.fit))^2),
               iter = iter))

}

# Source the example generator (cvx.generator)
source(paste(SOURCE_DIRECTORY,
             "code/tools/cvx-generator.R", sep = ""))
             
