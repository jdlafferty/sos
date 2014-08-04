################################################################################
# 
# 1-D Convex Regression with Rmosek
# 
# YJ Choe
# The University of Chicago
#
# Version 1.0: July 19, 2014
#
################################################################################

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
lambda = 25
#x <- sort(runif(n, -1, 1))# Need sorted points
x = seq(1:n)
x = x/n 
epsilon <- rnorm(n, 0, 1/3) # Gaussian noise
#theta = 10*(x^2) + 5*exp(x)
theta = 5*log(x + 0.2)  
y <- theta + epsilon
y <- y - rep(mean(y),n)
mat = matrix(0,n,2)
mat[,1] = rep(1,n)
mat[,2] = x
proj = mat%*%solve(t(mat)%*%mat)%*%t(mat)
y <- y - proj%*%y
plot(x, y, main = expression(y == x^2 + epsilon))

#---------------------------------------
# Convex Programming using Rmosek CQP
#---------------------------------------

# [Program variables]
# fitted values: z_1, ..., z_n
# auxiliary variables: w_1, ... , w_n
# objective replacement: t
f.index <- seq(1, n)
g.index <- seq(n+1, 2*n)
r.index <- seq(2*n + 1, 3*n)
t.index <- 3*n + 1
# Set up the program
convpattern.regression <- list(sense = "min")

# Objective: t = sqrt(sum((y-z)^2))
# Note that MSE = (1/n) * (t^2).
convpattern.regression$c <- c(rep(0, 3*n), 1)

A1 <- cBind(.sparseDiagonal(n),.sparseDiagonal(n),.sparseDiagonal(n))
# Add zero columns for unused variable t
A1 <- cBind(A1, Matrix(0, n, 1))


# Affine constraint 2: convexity [n-2 affine inequalities]
# (x_{i+1}-x_i)*z_{i-1} - (x_{i+1}-x_{i-1})*z_i 
# + (x_i-x_{i-1})*z_{i+1} >= 0
# (that is, slope_i <= slope_{i+1})
# for i = 2, ..., n-1
diag1 <- x[3:n] - x[2:(n-1)]
diag2 <- x[1:(n-2)] - x[3:n]
diag3 <- x[2:(n-1)] - x[1:(n-2)]
A2 <- bandSparse(n-2, n, c(0,1,2), list(diag1, diag2, diag3))
# Add zero columns for unused variables
#A2 <- cBind(A2, Matrix(0, n-2, n+1))
A2f <- cBind(A2, Matrix(0, n-2, 2*n + 1))
A2g <- cBind(Matrix(0, n-2, n), -A2 ,Matrix(0,n-2,n+1))
#Identifiability constraints
A3 = Matrix(0,4,3*n + 1)
A3[1,] = c(rep(1,n),rep(0,2*n + 1))
A3[2,] = c(rep(0,n),rep(1,n),rep(0,n+1))
A3[3,] = c(x,rep(0,2*n + 1))
A3[4,] = c(rep(0,n),x,rep(0,n+1))
##########################
pen <- rep(0,n)
pen[n] = 1/(x[n] - x[n-1])
pen[n-1] = -1/(x[n] - x[n-1])
pen[2] = -1/(x[2] - x[1])
pen[1] = 1/(x[2] - x[1])
pen = c(pen,-pen,rep(0,n+1))
convpattern.regression$A <- rBind(A1,A2f,A2g,A3,pen)

convpattern.regression$bc <- rbind(
  blc = c(y, rep(0, n-2), rep(0, n-2), rep(0,4),0),
  buc = c(y, rep(Inf, n-2), rep(Inf, n-2), rep(0,4),lambda)
)

# Empty constraints on the program variables
convpattern.regression$bx <- rbind(
  blx = c(rep(-Inf, 3*n), 0),
  bux = rep(Inf, 3*n + 1)
)


# Quadratic objective: mean squared error
# (replaced by an equivalent conic constraint)
convpattern.regression$cones <- cbind(
  list("QUAD", c(t.index,r.index))
)

# Solve the program using Rmosek!
r <- mosek(convpattern.regression)
fhat = r$sol$itr$xx[f.index]
ghat = r$sol$itr$xx[g.index]
fit = fhat + ghat
MSE <- (1/n) * (r$sol$itr$xx[t.index]^2)

# Plot
plot(x, y, main = "1-D Convex Regression using Rmosek")
lines(x,fit, lwd = 2, col = "red")
#mtext(paste(status$solution, ", ", status$program, sep=""), side = 3, adj = 0)
#mtext(paste("N: ", n, ", MSE: ", MSE, sep=""), side = 3, adj = 1)
plot(x, y, main = "1-D Convex Regression using Rmosek")
lines(x,fhat, lwd = 2, col = "green")
lines(x,ghat, lwd = 2, col = "blue")
MSEinc <- mean((y - fhat)^2)
MSEdec <- mean((y - ghat)^2)





