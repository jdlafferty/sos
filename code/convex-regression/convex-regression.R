########################################
# 
# 1-D Convex Regression with Rmosek
# 
# YJ Choe
# The University of Chicago
#
# Version 1.0: July 19, 2014
#
########################################

#---------------------------------------
# Setup
#---------------------------------------

# Import Rmosek
if (!require("Rmosek")) {
  stop ("Rmosek not installed.")
}

# Generate Synthetic Data
# True convex function: f(x) = x^4 + 2x + 3
n <- 100
x <- sort(runif(n, -5, 5)) # Need sorted points
epsilon <- rnorm(n, 0, 20) # Gaussian noise
y <- (x**4 + 2*x + 3) + epsilon
plot(x, y, main = expression(y == x^4 + 2*x + 3 + epsilon))

#---------------------------------------
# Convex Programming using Rmosek CQP
#---------------------------------------

# [Program variables]
# fitted values: z_1, ..., z_n
# auxiliary variables: w_1, ... , w_n
# objective replacement: t
z.index = seq(1, n)
w.index = seq(n+1, 2*n)
t.index = 2*n + 1

# Set up the program
convex.regression <- list(sense = "min")

# Objective: t = sqrt(sum((Y-Z)^2))
# Note that MSE = (1/n) * (t^2).
convex.regression$c <- c(rep(0, 2*n), 1)

# Affine constraint 1: auxiliary variables [no cost]
# z_i - w_i = y_i (that is, w_i = z_i - y_i)
# for i = 1, ..., n
A1 <- cBind(Matrix(diag(n)), -Matrix(diag(n)))
A1 <- cBind(A1, Matrix(0, n, 1)) # unused variable t

# Affine constraint 2: convexity [n-2 affine inequalities]
# (x_{i+1}-x_i)*z_{i-1} - (x_{i+1}-x_{i-1})*z_i 
# + (x_i-x_{i-1})*z_{i+1} >= 0
# (that is, slope_i <= slope_{i+1})
# for i = 2, ..., n-1
diag1 <- x[3:n] - x[2:(n-1)]
diag2 <- x[1:(n-2)] - x[3:n]
diag3 <- x[2:(n-1)] - x[1:(n-2)]
A2 <- Matrix(diag(diag1, n-2, n)) + 
      cBind(Matrix(0, n-2, 1), diag(diag2, n-2, n-1)) +
      cBind(Matrix(0, n-2, 2), diag(diag3, n-2, n-2))
A2 <- cBind(A2, Matrix(0, n-2, n+1))

convex.regression$A <- rBind(A1, A2)
convex.regression$bc <- rbind(
  blc = c(y, rep(0, n-2)),
  buc = c(y, rep(Inf, n-2))
)

# Empty constraints on the program variables
convex.regression$bx <- rbind(
  blx = c(rep(-Inf, 2*n), 0),
  bux = rep(Inf, 2*n+1)
)

# Quadratic objective: mean squared error
# (replaced by an equivalent conic constraint)
convex.regression$cones <- cbind(
  list("QUAD", c(t.index, w.index))
)

# Solve the program using Rmosek!
r <- mosek(convex.regression)

#---------------------------------------
# Results and Plots
#---------------------------------------

# Outputs
status <- list(
  solution = r$sol$itr$solsta, 
  program = r$sol$itr$prosta
)
fit <- r$sol$itr$xx[z.index]
MSE <- (1/n) * (r$sol$itr$xx[t.index]^2)

# Plot
plot(x, y, main = "1-D Convex Regression using Rmosek")
lines(x, fit, lwd = 2, col = "red")
mtext(paste(status$solution, ", ", status$program, sep=""), 
      side = 3, adj = 0)
mtext(paste("N: ", n, ", MSE: ", MSE, sep=""), side = 3, adj = 1)
