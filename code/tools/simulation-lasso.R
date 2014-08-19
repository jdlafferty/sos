################################################################################
# 
# Convexity Pattern Problem: Simulation with the Lasso Version
# 
# Professor John Lafferty's Group
# YJ Choe
# The University of Chicago
#
# ------------------------------------------------------------------------------
#
# Version 1.0: August 15, 2014
# - Performs a simulation and gives plots of the estimated success rate
#   (in finding correct convexity patterns) and the average running time
#   for the Lasso solution.
# - Warning: Running the function with default values may take a long time.
#
# ** Please feel free to improve the code and leave a note here. **
#
################################################################################

# Change this to the location of your local repository (must end with /)
SOURCE_DIRECTORY = "~/Code/sos-convexity/sos/"

# Source the 1-D version (convexity.pattern.lasso)
source(paste(SOURCE_DIRECTORY, 
             "code/convexity-pattern-lasso/convexity-pattern-lasso.R", sep = ""))

pattern.simulation = function (n.range = c(16, 32, 64, 128, 256, 512), 
                               p.range = c(5, 10, 15, 20), 
                               sigma = 0.5,
                               sparsity = function (p) { ceiling (0.4 * p^0.75) },
                               num.repeats = 20,
                               parallel = FALSE) {

  ########################################
  #
  # The main function.
  #
  # [Inputs]
  # n.range: vector of testing sample sizes
  # p.range: vector of testing dimensions
  # sigma: a fixed value for the size
  #        of Gaussian noise to data
  # sparsity: a function of p indicating
  #           the ratio of sparse components
  # num.repeats: number of repeats in
  #              estimating the probability
  #              of success
  # parallel: if TRUE, calls the R parallel
  #           package and performs the
  #           simulation in parallel
  #
  # [Output]
  # None, gives the resulting plot.
  #
  ########################################

  #---------------------------------------
  # Simulation
  #---------------------------------------
  
  success.rate = list()
  average.time = list()
  
  # A line per dimension p
  for (p in p.range) {
  	
  	num.sparse = sparsity(p)
  	success.rate[[p]] = rep(0, length(n.range))
  	average.time[[p]] = rep(0, length(n.range))
  	
    for (n in n.range) {
      
      n.index = which(n == n.range)
      success = rep(0, num.repeats)
      running.time = rep(0, num.repeats)
      cat(sprintf("Testing n=%d and p=%d...", n, p))
      
      # Specify lambda
      lambda = sigma * sqrt(log(p)/n)
      cat(sprintf("(lambda=%.2f)\n", lambda))
      
      for (iter in 1:num.repeats) {
        
        data = cvx.generator(n, p, sigma, 
                             num.convex = sample(seq(0, p-num.sparse), 1),
                             num.sparse = num.sparse)
        true.pattern = data$pattern
        
        t0 = proc.time()
        result = tryCatch(convexity.pattern.lasso(data$X, data$y, 
                                                  lambda = lambda, verbose = 1),
                          error = function (e) NULL)
        t1 = proc.time()
        running.time[iter] = t1["elapsed"] - t0["elapsed"]
        
        if (is.null(result)) {
          cat("Warning: Algorithm failed.\n")
          pattern = rep(-2, p)
        }
        
        # Error if a pattern has a component with both components active
        pattern = tryCatch(parse.pattern(result$pattern), 
                           error = function (e) { 
                           	         cat("Warning: Both components active.\n")
                           	         rep(-2, p)
                           	       }) 

        if (all(pattern == true.pattern)) {
          success[iter] = 1
        }
        cat(sprintf("n=%d, p=%d, iteration=%d/%d: %s\n", n, p,
                    iter, num.repeats, ifelse(success[iter], "success", "failure")))
        
      }
      
      success.rate[[p]][n.index] = mean(success)
      average.time[[p]][n.index] = mean(running.time) 
      
    }
  }

  #---------------------------------------
  # Plot
  #---------------------------------------
  
  colors = c("red", "blue", "green2", "purple", "brown", "black", "magenta",
             "darkgray", "orange", "cyan3", "darkblue", "yellowgreen")
  # the following should not happen, but...
  if (length(p.range) > length(colors)) {
    colors = rep(colors, ceiling(length(p.range)/length(colors)))
  }

  layout(mat = matrix(c(1, 2, 3, 3), nrow = 2, byrow = T), heights = c(0.4, 0.25))
  
  # Plot 1: Success Rate
  plot(NA, main = "Success Rate",
       xlim = c(min(n.range), max(n.range)), ylim = c(0, 1),
       xlab = paste("Sample size", expression(n)),
       ylab = sprintf("Success Rate", num.repeats))
  for (p in p.range) {
    lines(n.range, success.rate[[p]], 
          col = colors[which(p==p.range)], type = "b", lwd = 2, pch = 20)
  }

  # Plot 2: Average Running Time
  max.time = max(sapply(p.range, function (p) { max(average.time[[p]]) }))
  plot(NA, main = "Average Running Time",
       xlim = c(min(n.range), max(n.range)), ylim = c(0, max.time),
       xlab = paste("Sample size", expression(n)),
       ylab = sprintf("Average Running Time", num.repeats))
  for (p in p.range) {
    lines(n.range, average.time[[p]], 
          col = colors[which(p==p.range)], type = "b", lwd = 2, pch = 20)
  }
  
  # Legend & Text
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("top", inset = 0,
         sapply(p.range, function (p) { paste("p=", p, sep = "") }),
         col = colors, lwd = 2, pch = 20)
  mtext(sprintf("Number of Trials: %d", num.repeats), side = 3, adj = 0)
  mtext(sprintf("Noise Level: %.2f", sigma), side = 3, adj = 1)

}