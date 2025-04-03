###################################################
## CLASSIFICATION FUNCTIONS ##
###################################################


# Genetates n draws from a Dirichlet(alpha_1, ..., alpha_k) distribution
# Returns in a (k,n) matrix p
rdiric <- function(n, alpha, sumgam=F){
  k <- length(alpha)
  gam <- rgamma(k*n, rep(alpha,n))
  p <- matrix(gam, k, n)
  psum <- c(as.vector(rep(1, k)) %*% p)
  p <-  matrix(c(p) / c(t(matrix(rep(psum, k), n, k))), k, n)
  if(sumgam) return(list(p=p, psum=psum))
  else return(p)
}

###################################################

findprob <- function(alpha, N, nsamp){
  aN <- alpha + N
  phi <- array(0, dim=c(nsamp, ncol(aN), nrow(aN)))   # nsamp x M x T
  probTy <- matrix(0, nrow(aN), ncol(aN))
  for (i in (1:nrow(aN))) phi[,,i] <- t(rdiric(nsamp, aN[i,]))
  den <- apply(phi, c(1,2), sum)
  for (i in (1:nrow(aN))) probTy[i,] <- apply(phi[,,i] / den, 2, mean)
  dimnames(probTy) <- list(c("bulb", 
                             "car_window", 
                             "headlamp", 
                             "container",
                             "building_window"),
                           c(paste("cluster_", 1:ncol(aN), sep="")))
  return(probTy)
}


###################################################


