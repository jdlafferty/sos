################################################################################
# 
# Generic Plotting Function for Additive Components
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
# Version 1.0: August 7, 2014
# - A working version of an auxiliary plotting function for additive models.
#
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

plot.components = function (X, y, fit, pattern, true.components = NULL, 
                            display = NULL, file = NULL) {
  
  ########################################
  #
  # The main function.
  #
  # [Inputs]
  # X: n by p matrix of data points
  # y: n-vector of output values (truth + noise)
  # pattern: pattern as a p-vector
  #          1 convex, -1 concave, 0 sparse
  # fit: n by p matrix containing 
  #      the fitted values of each component
  # true.components: n by p matrix of true components.
  #                  this argument is optional.
  # display: list of indices for display.
  #          if not specified, the first 5 (at most)
  #          components are shown.
  # file: if specified, the plot is
  #       saved in file as pdf, jpeg, or png.
  #       must be a valid path ending with .pdf/.jpg/.png.
  #
  # [Output]
  # None, only plots are shown/saved.
  #
  ########################################

  #---------------------------------------
  # Setup
  #---------------------------------------

  # Require packages
  if (!require("Matrix")) {
    stop ("Matrix package not installed")
  }

  # Save in 'file' as pdf if specified
  if (!is.null(file)) {
  	extension = substring(file, nchar(file)-2, nchar(file))
  	switch(extension, pdf = pdf(file), jpg = jpeg(file), png = png(file))
  }

  # Variables
  n = nrow(X)
  p = ncol(X)
  y.hat = rowSums(fit)

  # Set list of covariates for display
  if (is.null(display)) {
    display = seq(1, max(p, 5))
  }

  #---------------------------------------
  # Plot
  #---------------------------------------  
  
  k = length(display)
  m = matrix(c(1:k, rep(k+1, k)), nrow = 2, byrow = TRUE)
  layout(mat = m, heights = c(0.4, 0.25))
  
  for (j in display) {
    
    par(mar = c(2, 4, 4, 2))
    
    ord = order(X[, j])
    
    plot(X[, j], y, main = sprintf("Component %d", j),
         xlab = expression(x[j]), ylab = expression(y))
    
    if (!is.null(true.components)) {
      lines(X[, j][ord], true.components[, j][ord], lwd = 3, col = "darkgray")
    }
    
    if (pattern[j] == 1) {
      lines(X[, j][ord], fit[, j][ord], lwd = 3, col = "red")
    }
    if (pattern[j] == -1) {
      lines(X[, j][ord], fit[, j][ord], lwd = 3, col = "blue")
    }
    
    points(X[, j], y.hat, pch = 20, col = "purple")
  }
  
  # Legend
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("top", inset = 0, 
         c("Data", "True component", "Convex component", 
           "Concave component", "Additive fit"), 
         col = c("black", "darkgray", "red", "blue", "purple"), 
         pch = c(21, NA, NA, NA, 20), lwd = c(NA, 3, 3, 3, NA))
  
  # Save to file if specified
  if (!is.null(file)) {
    dev.off()
  }
}
