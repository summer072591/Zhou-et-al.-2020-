
#' ******************************************************************************
#' File was written to approximate the solution to composite quantile regression
#' ******************************************************************************
#' 
#' ******************************************************************************
#'             the effective score function for comp quantile
#' ******************************************************************************
Phii_quan <- function(z,b, psi, nu){
  left <- as.numeric(b*(nu %*% (psi -1)))
  right <- as.numeric(b*(nu %*% psi))
  
  phi <- sapply(z, FUN = efsc_quan, b = b, left, right)
  return(phi)
}

efsc_quan <- function(z1, b, left, right){
  
  if(z1 > right) phi <- right
  else if(z1 < left) phi <- left
  else phi <- z1
  
  return(phi)
}
#' ******************************************************************************
#' 
#' 
#' ******************************************************************************
#'                      the thresholding function
#' ******************************************************************************
#' 
eta <- function(x,theta){
  
  x <- sapply(x, FUN = FUN_eta, theta = theta)
  return(x)
}

FUN_eta <- function(x1, theta){
  if(x1 > theta)
    x1 <- x1 - theta
  else if (x1 < (-theta))
    x1 <- x1 + theta
  else
    x1 <- 0
  
  return(x1)
}

#' ******************************************************************************
#' derivative of the thresholding function w.r.t. the first argument
#' ******************************************************************************

etaprime <- function(x, theta) return(sapply(x, function(x1, theta) ifelse(x1 > theta | x1 < -theta, 1, 0),
                                             theta = theta))
#' ******************************************************************************
#' 
#' 
#' ******************************************************************************
#'                       calculating values for b_t
#' ******************************************************************************

b_quan_emp <- function(C, s, n, psi, nu, bpsi, h, b_start, b_end, bpoint){
  
  size <- length(C)
  grid.b <- seq(b_start, b_end, length.out = bpoint)
  
  search.eq <- as.numeric(sapply(grid.b, FUN = b_grid_prep, psi = psi, nu = nu, bpsi = bpsi, C = C, h = h)) -s/n
 
  b.p <- grid.b[which(search.eq > 0)[1]]
  b.n <- grid.b[which(search.eq > 0)[1] -1]
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
}


#' ******************************************************************************
#' 
#' 
b_quan_simu <- function(W, Z, sigma, s, n, psi, nu, bpsi, h, b_start, b_end, bpoint){
  
  C <- as.vector(sapply(W, function(x) x + sigma * Z))
  size <- length(W)*length(Z)
 
  grid.b <- seq(b_start, b_end, length.out = bpoint)
  search.eq <- as.numeric(sapply(grid.b, FUN = b_grid_prep, psi = psi, nu = nu, bpsi = bpsi, C = C, h = h)) - s/n
  
  
  b.p <- grid.b[which(search.eq > 0)[1]]
  b.n <- grid.b[which(search.eq > 0)[1] -1]
  
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
}



b_grid_prep <- function(b, psi, nu, bpsi, C, h){
  
  left <- as.numeric(b*(nu %*% (psi -1)))
  right <- as.numeric(b*(nu %*% psi))
  K <- length(psi)
  mid <- Phi_mid(psi, nu, K)
  size <- length(C)
  F1 <- sum(ifelse(C < right|C == right, 1, 0))/size
  F2 <- sum(ifelse(C < left|C == left, 1, 0))/size
  
  
  f1 <- right * 1/(size*h)* sum(dnorm(C, mean = bpsi[K] + right, sd = h)) ## estimate right end point
  f2 <- left * 1/(size*h) * sum(dnorm(C, mean = bpsi[1] + left, sd = h)) ## estimate left end point
  
  rec_seq <- numeric(K-1)
  for(A in 1:(K-1)){
    f_mid_lt_prep <- dnorm(C, mean = bpsi[A] + b*mid[A], sd = h)
    f_mid_rt_prep <- dnorm(C, mean = bpsi[A+1] + b*mid[A], sd = h)
    rec_seq[A] <- 1/(size*h) * sum(b*mid[A]*(f_mid_rt_prep - f_mid_lt_prep))
  }
  
  est <- F1 - F2 - f1 + f2 + sum(rec_seq) 
  
  return(est)
  
}

#' ******************************************************************************
#'               theoretical state evolution functions
#' ******************************************************************************

sig_quan = function(X,Z,tau,alpha, delta){
  
  temp <- sum(sapply(X, function(x, z, tau, alpha, delta) 1 / delta * (eta(x + tau * z,
                                                                           (tau * alpha)) - x)^2, 
                     z = Z, tau = tau, alpha = alpha, delta = delta))
  
  sig= (temp/length(X)/length(Z))
  return(sig)
}



tau_quan = function(W,Z,sig,b,w,delta, psi, bpsi, nu){
  
  temp <- sum(sapply(W, function(W, Z, sig, b, w, delta, psi, bpsi, nu) (delta/w*Phii_quan(W + sig*Z,b, psi, bpsi, nu))^2,
                     Z = Z, sig = sig, b = b, w = w, delta = delta, psi = psi, bpsi = bpsi, nu = nu))
  return((temp/length(W)/length(Z)))
  
} 

Phii_quan <- function(z, b, psi, bpsi, nu){
  K <- length(psi)
  lnd <- as.numeric(-nu %*% (1 - psi))
  rnd <- as.numeric(nu %*% psi)
  f <- Phi_mid(psi, nu, K)
  
  
  left_diff <- bpsi[1:(K-1)] + b*f
  right_diff <- bpsi[2:K] + b*f 
  
  left_non_diff <- c(bpsi[1] + b*lnd, bpsi[2:K] + b*f[1:(K-1)])
  right_non_diff <- c(bpsi[1:(K-1)] + b*f[1:(K-1)], bpsi[K] + b*rnd)
  
  
  
  phi <- sapply(z, FUN = Phii_general,
                b = b, psi = psi, bpsi = bpsi, nu = nu, K = K, lnd = lnd, rnd = rnd, 
                f = f, left_diff = left_diff, right_diff = right_diff, left_non_diff = left_non_diff,
                right_non_diff = right_non_diff)
  
  
  return(as.numeric(phi))
  
}

# Phi_mid calculates the points in the middel of the two end points
# this function corresponds to the function f(A)

Phi_mid <- function(psi, nu, K){
  
  tmp1 <- cumsum((nu * psi)[1:K-1])
  tmp2 <- rev((cumsum(rev(nu * (1 - psi))[1:K-1])))
  
  
  return(tmp1 - tmp2) 
}



Phii_general <- function(z, b, psi, bpsi, nu, K, lnd, rnd, f, left_diff, right_diff, left_non_diff, right_non_diff){
  
  
  pointer1 <- which(z >= left_non_diff & z <= right_non_diff)
  pointer2 <- which(z >= left_diff & z <= right_diff)
  
  
  if(length(pointer1)!=0){
    phii = z - bpsi[pointer1]
  }else if(length(pointer2)!=0){
    phii = b*(f[pointer2])
  }else if(z < bpsi[1] + b*lnd){
    phii = b*lnd
  }else if(z > bpsi[K] + b * rnd){

    phii = b*rnd
  }
  
  return(phii)
}

#' ******************************************************************************
##                     lower bound of the MSE 
#' ******************************************************************************

alpha_lower <- function(delta){
  
  fun <- function(alpha) (1 + alpha^2) * pnorm(-alpha) - alpha * dnorm(alpha) - delta/2 
  
  alpha.lower <- uniroot(fun, interval = c(0,5))$root
  
  return(alpha.lower)
}


#' ******************************************************************************
##                       the AMSE function
#' ******************************************************************************

amse <- function(tau2, alpha, x){
  
  tmp.gammas <- gammas(tau = sqrt(tau2), alpha, x)
  tmp <- tau2*tmp.gammas[1] + tmp.gammas[2]
  
  return(tmp)
}

#' ******************************************************************************
#' 
gammas <- function(tau, alpha, x){
  
  tmp.gamma.1 <- as.matrix(sapply(x, FUN = gamma.sub, alpha = alpha, tau = tau))
  
  gamma1 <- apply(tmp.gamma.1, 1, mean)
  
  gamma1[1] <- 1 + alpha^2 + gamma1[1]
  return((gamma1))
}

gamma.sub <- function(tau, alpha, x1){
  
  left <- alpha - x1/tau
  right <- alpha + x1/tau 
  
  tmp1 <- -right*dnorm(left) - left*dnorm(right) + (1 + alpha^2)*(pnorm(-right) - pnorm(left))
  
  tmp2 <- x1^2*(pnorm(left) - pnorm(-right))
  return(c(tmp1, tmp2))
}

#' ******************************************************************************

