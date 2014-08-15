################################################################################
# 
# Convexity Pattern Problem, Lasso Version: Simulation
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

# Change this to the location of your local repository (must end with /)
SOURCE_DIRECTORY = "~/Code/sos-convexity/sos/"

# Source the 1-D version (convexity.pattern.cvx.1d.auto)
source(paste(SOURCE_DIRECTORY, 
             "code/convexity-pattern-cvx/convexity-pattern-cvx-1d.R", sep = ""))
             
# Source the backfitting version (convexity.pattern.cvx.backfit)
source(paste(SOURCE_DIRECTORY, 
             "code/convexity-pattern-cvx/convexity-pattern-cvx-backfit.R", 
             sep = ""))
             

pattern.simulation = function (n.range = seq(100, 1000, 100), 
                               p.range = seq(4, 20, 4), 
                               sigma = 0.5,
                               sparsity = function (p) { ceiling (0.4 * p^0.75) },
                               num.repeats = 200,
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
  # None, gives a resulting plot.
  #
  ########################################

  #---------------------------------------
  # Setup
  #---------------------------------------
  
  if (parallel) {
    worker = function (params) {    	
      data = cvx.generator(params$n, params$p, params$sigma, 
                           num.convex = sample(seq(0, params$p-params$num.sparse),
                                               1),
                           num.sparse = params$num.sparse)
      true.pattern = data$pattern
      
      t0 = proc.time()
      result = tryCatch(convexity.pattern.cvx.backfit(data$X, data$y, silent=T),
                        error = function (e) NULL)
      t1 = proc.time()
      
      if (is.null(result)) return (NULL)
      
      pattern = parse.pattern(result$pattern)
      
      return (list(success = all(pattern == true.pattern),
                   running.time = t1["elapsed"] - t0["elapsed"]))
    }
  }
  
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
      
      #cat(sprintf("Testing n=%d and p=%d...\n", n, p))
      
      # Parallel version
      if (parallel) {
        params = list(n = n, p = p, sigma = sigma, num.sparse = num.sparse)
        out.parallel = unlist(mclapply(rep(params, num.repeats), worker))
        next
      }
      
      for (iter in 1:num.repeats) {
        
        data = cvx.generator(n, p, sigma, 
                             num.convex = sample(seq(0, p-num.sparse), 1),
                             num.sparse = num.sparse)
        true.pattern = data$pattern
        
        t0 = proc.time()
        result = tryCatch(convexity.pattern.cvx.backfit(data$X, data$y, step=0.05,
                                                        max.step=100, silent=F),
                          error = function (e) NULL)
        t1 = proc.time()
        running.time[iter] = t1["elapsed"] - t0["elapsed"]
        
        if (is.null(result)) {
          cat("Warning: Algorithm failed.\n")
          next
        }
        pattern = parse.pattern(result$pattern)
        print (list(pattern=result$pattern, true.pattern=true.pattern))
        if (all(pattern == true.pattern)) {
          success[iter] = 1
        }
        cat(sprintf("n=%d, p=%d, iteration=%d/%d: %d\n", n, p,
                    iter, num.repeats, success[iter]))
        
      }
      
      success.rate[[p]][n.index] = mean(success)
      average.time[[p]][n.index] = mean(running.time) 
      
    }
  }

  #---------------------------------------
  # Plot
  #---------------------------------------
  
  colors = rainbow(length(p.range))

  par(mfrow = c(1, 2))
  
  plot(NA, main = "Success Rate",
       xlim = c(min(n.range), max(n.range)), ylim = c(0, 1),
       xlab = paste("Sample size", expression(n)),
       ylab = sprintf("Success Rate in %d Trials", num.repeats))
  for (p in p.range) {
    lines(n.range, success.rate[[p]], 
          col = colors[which(p==p.range)], type = "b")
  }

  plot(NA, main = "Average Running Time",
       xlim = c(min(n.range), max(n.range)), ylim = c(0, 100),
       xlab = paste("Sample size", expression(n)),
       ylab = sprintf("Average Running Time in %d Trials", num.repeats))
  for (p in p.range) {
    lines(n.range, average.time[[p]], 
          col = colors[which(p==p.range)], type = "b")
  }

  return (list(success.rate = success.rate,
               average.time = average.time))
}