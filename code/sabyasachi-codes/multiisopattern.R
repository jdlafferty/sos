#This program attempts an idea similar to Liso of Meinshausen
# Import Rmosek
if (!require("Rmosek")) {
  stop ("Rmosek not installed.")
}
set.seed(690)
lambda = 3
n <- 100
p = 4
sigma = 1/3
X <- Matrix(runif(n*p, -1, 1), nrow=n) 
#epsilon <- rnorm(n, 0, ) # Gaussian noise
# Define an non decreasing function named myfn.
f <- list()
f[[1]] = function(x){
  #return(2*x)
  return(x^3 + 1)
}
f[[2]] = function(x){
  ifelse ((x > 0),-x^2,0)
}
f[[3]] = function(w){
  y = abs(w)
  ans = y^{0.2}
  ans = (w/y)*ans
  return(ans)
}

f[[4]] = function(w){
  ans = w + abs(w)
  return(-ans/2)
}
#Vectorize(myfn4)
y <- rowSums(sapply(1:p, function (j) {
  f_j <- sapply(X[,j], f[[j]])
  return (f_j - mean(f_j))
})) 
y <- y + rnorm(n, 0, sigma) # Gaussian noise
#plot(x, y, main = expression(y == x^4 + 2*x + 3 + epsilon))
#plot(x, y, main = "Isotonic sequence")
#---------------------------------------

# Convex Programming using Rmosek CQP
#---------------------------------------
last <- function (x) {
  # x is a sequence
  return (x[length(x)])
}

# [Program variables]
# fitted increasing values: f_1,..,f_n.
# fitted decreasing values: g_1,..,g_n.
# auxiliary variables: r_i = y_i - (f_i + g_i)
# objective replacement: t
f.index <- seq(1, n*p) 
g.index <- last(f.index) + seq(1, n*p) 
r.index <- last(g.index) + seq(1, n) 
t.index <- last(r.index) + 1 
num.vars <- t.index

# Helpers: Selecting variables
# (Think of f_ij as the (i,j)th entry of an n by p matrix.)
f.select.row <- lapply(1:n, function (i) { seq(i*p-p+1, i*p) })
f.select.col <- lapply(1:p, function (j) { seq(j, n*p, p) })
g.select.row <- lapply(1:n, function (i) { last(f.index) + seq(i*p-p+1, i*p) })
g.select.col <- lapply(1:p, function (j) { last(f.index) + seq(j, n*p, p) })
# Set up the program
isopattern.regression <- list(sense = "min")

# Objective: t = sqrt(sum((y-z)^2))
# Note that MSE = (1/n) * (t^2).
isopattern.regression$c <- c(rep(0, 2*n*p + n), 1)

# Affine constraint 1: auxiliary variables [no cost]
# r_i = y_i - sum_j(f_ij + g_ij)
A1 <- Matrix(0, nrow = n, ncol = num.vars)
for (i in 1:n) {
  A1[i, f.select.row[[i]]] <- 1
  A1[i, g.select.row[[i]]] <- 1
  A1[i, r.index[i]] <- 1
}

A2 <- Matrix(0, nrow = 2*(n-1)*p, ncol = num.vars)
for (j in 1:p) {
  # take the ordering (for sorting and putting back in order)
  ord <- order(X[, j])
  # the following sequence x is sorted
  x <- X[ord, j]
  diag.f <- bandSparse(n-1, n, c(0, 1), list(rep(-1, n-1), rep(1, n-1)))
  diag.g <- -diag.f
  rows.f <- 2*(n-1)*(j-1) + seq(1, n-1)
  rows.g <- 2*(n-1)*(j-1) + (n-1) + seq(1, n-1)
  A2[rows.f, f.select.col[[j]]] <- diag.f[, invPerm(ord)]
  A2[rows.g, g.select.col[[j]]] <- diag.g[, invPerm(ord)]
}


# Affine constraint 3: identifiability via centering [2*p equalities]
A3 <- Matrix(0, nrow = 2*p, ncol = num.vars)
for (j in 1:p) {
  A3[j, f.select.col[[j]]] <- 1
  A3[p+j, g.select.col[[j]]] <- 1
}

#Penalty constraints on each variable
A4 <- Matrix(0, nrow = 1, ncol = num.vars)
for (j in 1:p) {
  ord <- order(X[, j])
  A4[1,f.select.col[[j]][ord[n]]] = 1
  A4[1,f.select.col[[j]][ord[1]]] = -1
  A4[1,g.select.col[[j]][ord[n]]] = -1
  A4[1,g.select.col[[j]][ord[1]]] = 1
}
isopattern.regression$A <- rBind(A1,A2,A3,A4)

isopattern.regression$bc <- rbind(
  blc = c(y, rep(0, 2*(n-1)*p), rep(0, 2*p), 0),
  buc = c(y, rep(Inf, 2*(n-1)*p), rep(0, 2*p), lambda)
)

# Empty constraints on the program variables
isopattern.regression$bx <- rbind(
  blx = c(rep(-Inf, 2*n*p + n), 0),
  bux = rep(Inf, 2*n*p + n + 1)
)

# Quadratic objective: mean squared error
# (replaced by an equivalent conic constraint)
isopattern.regression$cones <- cbind(
  list("QUAD", c(t.index,r.index))
)

# Solve the program using Rmosek!
r <- mosek(isopattern.regression)

fhat = r$sol$itr$xx[f.index]
ghat = r$sol$itr$xx[g.index]
fit = fhat + ghat
MSE <- (1/n) * (r$sol$itr$xx[t.index]^2)

f1hat = fhat[f.select.col[[1]]]
f2hat = fhat[f.select.col[[2]]]
f3hat = fhat[f.select.col[[3]]]
f4hat = fhat[f.select.col[[4]]]

f1hat = f1hat[order(X[,1])]
f2hat = f2hat[order(X[,2])]
f3hat = f3hat[order(X[,3])]
f4hat = f4hat[order(X[,4])]


g1hat = ghat[g.select.col[[1]] - last(f.index)]
g2hat = ghat[g.select.col[[2]] - last(f.index)]
g3hat = ghat[g.select.col[[3]] - last(f.index)]
g4hat = ghat[g.select.col[[4]] - last(f.index)]

g1hat = g1hat[order(X[,1])]
g2hat = g2hat[order(X[,2])]
g3hat = g3hat[order(X[,3])]
g4hat = g4hat[order(X[,4])]






# Sample plot in each component
m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.4, 0.2))

# Component 1
par(mar = c(2,4,4,2))
plot(x1, y, 
     main = "Convexity Pattern Regression (1)", xlab = expression(x[1]))
# sort by order of x1 for drawing lines
ord1 <- order(x1)
true.f <- sapply(x1, f)
lines(x1[ord1], (true.f - mean(true.f))[ord1], lwd = 5, col = "darkgray")
lines(x1[ord1], f1[ord1], lwd = 2, col = "red") # convex comp
lines(x1[ord1], g1[ord1], lwd = 2, col = "blue") # concave comp
points(x1, fit, pch = 20, col = "purple")

# Component 2
par(mar = c(2,4,4,2))
plot(x2, y, 
     main = "Convexity Pattern Regression (2)", xlab = expression(x[2]))
# sort by the order of x2 for drawing lines
ord2 <- order(x2)
true.g <- sapply(x2, g)
lines(x2[ord2], (true.g - mean(true.g))[ord2], lwd = 5, col = "darkgray")
lines(x2[ord2], f2[ord2], lwd = 2, col = "red") # convex comp
lines(x2[ord2], g2[ord2], lwd = 2, col = "blue") # concave comp
points(x2, fit, pch = 20, col = "purple")

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
}
