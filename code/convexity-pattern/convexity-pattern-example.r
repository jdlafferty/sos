
soft.threshold = function(x,lambda) {
  return(sign(x) * (abs(x)-lambda) * (abs(x) > lambda))
}

generator1 = function(n=150,d=50,sigma2=1) {
  epsilon = rnorm(n,0,sigma2)
  M = matrix(NA,nrow=n,ncol=d)
  x = matrix(runif(d*n,-2.5,2.5),nrow=n,ncol=d)
  for (j in c(1:d)){
    if (j==1) { M[,j] = -x[,j]^2 - mean(-x[,j]^2) }
    if (j==2) { M[,j] = x[,j]^2 - mean(x[,j]^2) }
    if (j==3) { M[,j] = x[,j] }
    if (j==4) { M[,j] = exp(-x[,j])- mean(exp(-x[,j])) }
    if (j>4 ) { M[,j] = 0 }
  }
  for (j in c(1:d)) {
    x[,j] = (x[,j]-min(x[,j]))/(max(x[,j])-min(x[,j]))
  }
  M[,1:4]  = scale(M[,1:4])          # rescale the design matrix
  m        = apply(M[,c(1:d)],1,sum) # mean function
  y        = m + epsilon
  return(list(y=y,x=x,M=M))
}

kernel.smoothing.matrix = function(h, x) {
  n = length(x)
  S = matrix(0,n,n)
  for (i in 1:n) {
    tmp1  = x - x[i]
    w     = exp(-(tmp1)^2/(2*h^2))
    w     = w/sum(w)
    S[i,] = w
  }
  return(S)
}

local.linear.smoothing.matrix = function(h, x) {
  n = length(x)
  S = matrix(0,n,n)
  
  for (i in 1:n) {
    xdiff = x - x[i]
    w     = exp(-(xdiff)^2/(2*h^2))
    W     = diag(w/sum(w))
    X     = cbind(1,matrix(xdiff,n,1,byrow=TRUE))
    T     = solve(t(X) %*% W %*% X) %*% t(X) %*% W
    S[i,] = T[1,]
  }
  return(S)
}


plotM = function(y,x,M1,M) {
  d   = ncol(M)
  res = y - apply(M1,1,sum)
  par(mfrow=c(1,5))
  for (i in c(1:min(d,10))) {
    par(cex=1.0,cex.sub=1.0,cex.lab=1.0)
    ord   = order(x[,i])
    M1[,i]= M1[,i]- mean(M1[,i])
    M[,i] = M[,i]- mean(M[,i])
    l = sum(abs(M1[,i]))
    if (l == 0) {
      s = "zero"
    }
    else {
      s = sprintf("l1=%.2f", sum(abs(M1[,i])))
    }
    plot(x[,i][ord], M1[,i][ord],ylim=c(min(y),max(y)),
         type='l',lwd=2,xlab=paste('x',i,sep=''),
         ylab=paste('m',i,sep=''),col='blue',main=s)
    points(x[,i][ord],y[ord],lwd=1)
    lines(M[,i][ord]~x[,i][ord],lty=2,lwd=2,col=2)
  }

  num.zeros = 0
  for (i in 1:d) {
    M1[,i] = M1[,i] - mean(M1[,i])
    l1 = sum(abs(M1[,i]))
    if (l1 == 0) {
      num.zeros = num.zeros + 1
    }
  }
  cat(sprintf("Number of zeros: %d of %d\n", num.zeros, d))
}

backfit = function(x,y,h,lambda,maxStep = 100){
  d  = ncol(x)
  n  = length(y)
  M1 = matrix(0,n,d)
  alpha = rep(1,d)
  S  = matrix(0,d,n*n)

  cat("constructing smoothing matrices...\n")
  for (i in 1:d) {
    #W = kernel.smoothing.matrix(h[i],x[,i])
    W = local.linear.smoothing.matrix(h[i],x[,i])
    S[i,] = as.vector(W)
  }
  
  cat("lambda=",lambda,"\n")
  order = sample(d,d)  # visit the variables in a random order
  cat("Variable order: ", order[1:min(length(order),8)], "...\n")

  for (iter in c(1:maxStep)){
    curfit    = apply(M1 %*% diag(alpha), 1, sum)
    for(i in order) {
      Z = y - M1[,-i] %*% alpha[-i]
      Si = matrix(S[i,],n,n)
      f = as.vector(Si %*% Z)
      norm = sqrt(t(f) %*% f / n)
      M1[,i] = f / norm
      M1[,i] = M1[,i] - mean(M1[,i])
      alpha[i] = soft.threshold(norm, lambda)
    }
    newfit    = apply(M1 %*% diag(alpha),1,sum)
    tmp       = sqrt(sum((curfit-newfit)^2))
    if(tmp<1e-6){ break }
    cat(sprintf("Iteration %d: l2-distance=%f\n", iter, tmp))
    if(iter==maxStep){cat('maximum steps reached!','\n')}
  }
  
  return(list(M=(M1 %*% diag(alpha)), beta=alpha))
}


example1 = function(n=150,d=50,lambda=.4,sigma2=1,h0=.1) {
  out = generator1(n,d,sigma2)
  y = out$y
  x = out$x
  M = out$M
  n = length(y)
  d = ncol(x)

  h = rep(1,d)
  if (h0 != 0) { h = rep(h0,d) }
  else {h =rep(.8*h/(n^0.2),d)}  # the plug-in bandwidths
  cat("bandwidth: ",h[1],"\n")

  tmp = backfit(x,y,h,lambda,maxStep=100)
  Mc = tmp$M
  plotM(y,x,Mc,M)
  #return(list(y=y,x=x,M1=M1,M=M))
}

example1(n=150,d=5,lambda=.1,sigma2=1,h0=.2)

