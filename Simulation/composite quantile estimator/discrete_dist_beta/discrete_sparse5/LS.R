#' **************************************************************************************************
#' File was written to approximate solutions to l_1 regularization with single loss functions
#' 
#' Please only use "quant" for loss.fun
#' **************************************************************************************************


RAMP.single <- function(A, y, a, s, T, tol = 10^(-4), W, Z, psi = .5, x, h = 4, 
                 bpoint = 500, b.start = 0, b.end = 5.5, loss.fun = c("abs", "quant")){
  
  # xhat = reconstructAmp(A, y, T, tol, mse, verbose)
  #
  #   Arguments:
  #       A - Design matrix
  #       y - response
  #       T - max number of iterations (optional)
  #       tol - stopping criteria (optional)
  #       x - original vector used to print progress of MSE (optional), true value
  #       h - Bandwidth for the consistent estimator for f_z()
  #       psi - quantile level of single quantile loss function
  #       bpoint, b.start, b.end - grid search for the parameter b in the effective score function, 
  #                                bpoint is the total number of points set between b.start and b.end
  #       loss.fun - absolute deviation and quantile; 
  #                  but the absolute deviation loss function was not completed
  #                  only use "quant"!!!
  #       W, Z - random variables in the limit equations
  
  # n, p, s define
  if(missing(s)) s <- length(which(x != 0))
  n <- dim(A)[1]
  p <- dim(A)[2]
 
  x0 <- rep(0, p)
 
  
  omega <- s/p
  delta <- n/p

  ##------------------------------------choose b0--------------------------------------##
  
  sigma02 <- sum(x^2)/(p * delta)
  sigma02.theory <- mean(x^2)/delta
  
  
  
  if(loss.fun == "abs"){
    z0 <- y - A %*% x0
    
    
    b0 <- b_abs(z0, s, n, h, bpoint, b.start, b.end)
    if(length(b0) == 0) b0 <- NA

    div <- 1
    while(is.na(b0)){
        b0 <- b_abs(z0, s, n, h = h*0.1/div, bpoint, b.start, b.end)
        if(length(b0) == 0) b0 <- NA
        div <- div + 1
      
        if(div >= 6) break
    }
    
    if(is.na(b0)){
        b0 <- b_abs_simu(W,Z,sqrt(sigma02.theory), s, n, h, bpoint, b.start, b.end)
        if(length(b0) == 0) b0 <- NA
    }
    
    if(is.na(b0)) b0 <- 0.2
    
  
  
  
  ##-----------------------------------------------------------------------------------------------------
  
    tau02 <- (n/s)^2*mean((sapply(z0, function(x) Phi_abs(x, b0)))^2)   #### p to s
    tau02.theory <- tau_abs(W,Z,sqrt(sigma02.theory),b0,omega,delta =delta)
  
    theta0 <- a*sqrt(tau02)
    arg_eta <- as.vector(x0 + t(A) %*% (n/s*sapply(z0, function(x) Phi_abs(x, b0))))  # p to s
    x.new <- eta(arg_eta, theta = theta0)
  
  ##------------------------------------------------------------------------------------------------------  
  
  
    omega.seq <- bvec <- tau2 <- tau2.theory <- sigma2 <- sigma2.theory <- theta <- 
      mse.est2 <- mse <- mse.est <- mse.est3 <- rep(NA, T)
    omega.seq <- c(omega, omega.seq)
    bvec <- c(b0, bvec)
    tau2 <- c(tau02, tau2)
    tau2.theory <- c(tau02, tau2.theory)
    sigma2 <- c(sigma02, sigma2)
    sigma2.theory <- c(sigma02.theory, sigma2.theory)
    theta <- c(theta0, theta)
    z.old <- z0
   
    for(t in 2:(T+1)){
      
      onsage.sca <- etaprime(arg_eta, theta = theta[t-1])
    
    
      ord.z <- y - A %*% x.new
      omega.seq[t] <- (sum(x.new != 0)/p)
      z.new <- y - A %*% x.new + (sapply(z.old, function(x) Phi_abs(x,bvec[t-1])))*mean(onsage.sca)/omega

      sigma2[t] <- sum((x.new - x)^2)/(delta*p)
    
      sigma2.theory[t] <- sig(x,Z,sqrt(tau2.theory[t-1]),a, delta=delta)
     
    
    ##------------------------------------update b--------------------------------------##
    
   
    inc <- 0
    while(is.na(bvec[t])){
      bvec[t] <- b_abs(as.vector(z.new), s = s, n, h = h + inc*0.01, bpoint, b.start, b.end)
      if(length(bvec[t]) == 0) bvec[t] <- NA
      inc <- inc + 1
      if(inc >= 6) break
    }
   
   if(is.na(bvec[t])){
       bvec[t] <- b_abs_simu(W, Z, sqrt(sigma2.theory[t]), s, n, h, bpoint, b.start, b.end)
       if(length(bvec[t]) == 0) bvec[t] <- NA
   }
   
    if(is.na(bvec[t])) bvec[t] <- bvec[t - 1]
   
    ##------------------------------------ estimation-----------------------------------##
   
    arg_eta <- x.new + t(A) %*% (delta/omega*sapply(z.new, function(x) Phi_abs(x, bvec[t])))
    
    ##------------------------------------------------------------------------------------------------------
    
    tau2.theory[t] <- tau_abs(W,Z,sqrt(sigma2.theory[t]),bvec[t],omega,delta)
    
  
    tau.tmp <- 0
    for(i in 1:5){
      z.temp <- sample(z.new, n, replace = TRUE)
      tau2.tmp <- sapply(z.temp, function(x) (delta/omega)*Phi_abs(x, bvec[t]))
      tau.tmp <- tau.tmp + mean(tau2.tmp^2)
    }
    tau2[t] <- tau.tmp/5
    
    theta[t] <- a*sqrt(tau2[t])

    x.old <- x.new
    x.new <- eta(arg_eta, theta = theta[t])   
    mse[t] <- (sum((x.new - x)^2)/p)
    
    
    mse.est[t] <- tau2[t]*2*mean(etaprime(arg_eta, theta = theta[t])) - tau2[t] + mean((arg_eta - x.new)^2)
    
    
    ## theoretical expression
   
    mse.est3[t] <- amse(tau2[t], alpha = a, x)
    z.old <- z.new
    
    if(((sum((x.old - x.new)^2)/sum(x.old^2))) < tol){
      T_max <- t
      break
    }else{
      T_max <- t
    }
    }
    
    lambda <- a*sqrt(tau2[T_max])*omega/(bvec[T_max]*delta)
  }else if(loss.fun == "quant"){
    
    z0 <- y - A %*% x0
    b0 <- b_quant(z0, s, n, psi, h, bpoint, b.start, b.end)
    div <- 1
    while(is.na(b0)){
       b0 <- b_quant(z0, s, n, psi, h = h*0.5/div, bpoint, b.start, b.end)
        div <- div + 1
      
        if(div >= 2) break
             }
    
    if(is.na(b0)) b0 <- b_quant_simu(W,Z,sqrt(sigma02.theory), psi, s, n, h, bpoint, b.start, b.end)
    if(is.na(b0)) b0 <- 0.2
    
    ##-------------------------Estimate the coefficient vector---------------------------##
    
    
    ##------------------------------------------------------------------------------------------------------
   
    ## next is the possible solution to the state evolution tau^2

    tau02 <- (n/s)^2*mean((sapply(z0, function(x) Phi_quant(x, b0, psi)))^2)   #### p to s
    tau02.theory <- tau_quant(W,Z,sqrt(sigma02.theory),b0, psi, omega,delta =delta)
    
    theta0 <- a*sqrt(tau02)
    arg_eta <- as.vector(x0 + t(A) %*% (n/s*sapply(z0, function(x) Phi_quant(x, b0, psi))))  # p to s
    x.new <- eta(arg_eta, theta = theta0)
    
    ##------------------------------------------------------------------------------------------------------  
    
    tol.seq <- NA
    omega.seq <- bvec <- tau2 <- tau2.theory <- sigma2 <- sigma2.theory <- theta <- 
      mse.est2 <- mse <- mse.est <- mse.est3 <- rep(NA, T)
    omega.seq <- c(omega, omega.seq)
    bvec <- c(b0, bvec)
    tau2 <- c(tau02, tau2)
    tau2.theory <- c(tau02, tau2.theory)
    sigma2 <- c(sigma02, sigma2)
    sigma2.theory <- c(sigma02.theory, sigma2.theory)
    theta <- c(theta0, theta)
    z.old <- z0
    
    for(t in 2:(T+1)){
      onsage.sca <- etaprime(arg_eta, theta = theta[t-1])
      

    
      ord.z <- y - A %*% x.new
      omega.seq[t] <- (sum(x.new != 0)/p)
     
      
      z.new <- y - A %*% x.new + (sapply(z.old, function(x) Phi_quant(x,bvec[t-1], psi)))*mean(onsage.sca)/omega
      
     
      sigma2[t] <- sum((x.new - x)^2)/(delta*p)
     
      sigma2.theory[t] <- sig(x,Z,sqrt(tau2.theory[t-1]),a, delta=delta)
      
      
      ##------------------------------------update b--------------------------------------##
      bvec[t] <- b_quant(as.vector(z.new), s = s, n, psi, h, bpoint, b.start, b.end)
      inc <- 0
      while(is.na(bvec[t])){
        bvec[t] <- b_quant(as.vector(z.new), s = s, n, psi, h = h + 0.5*inc, bpoint, b.start, b.end)
        
        inc <- inc + 1
        if(inc >= 3) break
      }
      
      if(is.na(bvec[t])){
           bvec[t] <- b_quant_simu(W, Z, sqrt(sigma2.theory[t]), psi, s, n, h, bpoint, b.start, b.end)
       }
      
      if(is.na(bvec[t])) bvec[t] <- bvec[t - 1]
    

      ##------------------------------------ estimation-----------------------------------##
    
      arg_eta <- x.new + t(A) %*% (delta/omega*sapply(z.new, function(x) Phi_quant(x, bvec[t], psi)))
      
      
      ##------------------------------------------------------------------------------------------------------
      
      tau2.theory[t] <- tau_quant(W,Z,sqrt(sigma2.theory[t]),bvec[t],psi, omega,delta)
      
      
      tau2.tmp <- sapply(z.new, function(x) (delta/omega)*Phi_quant(x, bvec[t], psi))
      tau2[t] <- mean(tau2.tmp^2)
      
      
      theta[t] <- a*sqrt(tau2[t])
      
      x.old <- x.new
      x.new <- eta(arg_eta, theta = theta[t])   
      mse[t] <- (sum((x.new - x)^2)/p)
      
      
      mse.est[t] <- max(tau2[t]*2*mean(etaprime(arg_eta, theta = theta[t])) - tau2[t] + mean((arg_eta - x.new)^2), 0)
      
      
      ## theoretical expression
      
      mse.est3[t] <- amse(tau2[t], alpha = a, x)
      z.old <- z.new
      
      
      tol.tmp <- mean((x.old - x.new)^2)
      tol.seq <- c(tol.seq, tol.tmp)
      if(tol.tmp < tol){
        T_max <- t
        break
      }else{
        T_max <- t
      }
    }
    lambda <- a*sqrt(tau2[T_max])*omega/(bvec[T_max]*delta)
  
  }
  
  return(list(x.new = x.new, arg_eta = as.numeric(arg_eta), bvec = bvec[1:T_max], omega.seq = omega.seq[1:T_max], sigma2 = sigma2[1:T_max], 
              tau2 = tau2[1:T_max], tau2.theory = tau2.theory[1:T_max], 
              sigma2.theory = sigma2.theory[1:T_max], theta= theta[1:T_max], mse = mse[1:T_max], 
              mse.est = mse.est[1:T_max], #mse.est2 = mse.est2[1:T_max], 
              mse.est3 = mse.est3[1:T_max], tol.seq = tol.seq, lambda = lambda, T_max = T_max))
}




## effective score of absolute deviation
Phi_abs <- function(z, b){
  if(z > -b & z < b) tmp <- z
  else tmp <- b*sign(z)
   
  return(tmp)
}

## effective score of single quantile

Phi_quant <- function(z, b, psi){
  if(z > b*psi) tmp <- b*psi
  else if(z < b*(psi -1)) tmp <- b*(psi -1)
  else tmp <- z
  
  return(tmp)
}


## calulate the value for b
b_abs <- function(z, s, n, h, bpoint, b.start, b.end){
  

  size <- length(z)
  C <- z
  grid.b <- seq(b.start, b.end, length.out = bpoint)
  
  search.eq <- lapply(grid.b, FUN = b_grid_prep_abs, C = C, h = h)
  
  search.eq <- lapply(search.eq, function(x, size, h, n){
    mat <- matrix(x, nrow = 4, byrow = TRUE)
    tmp1 <- sum(mat[1, ])/size
    tmp2 <- sum(mat[2, ])/size
    tmp3 <- 1/(size*h) * sum(mat[3,])
    tmp4 <- 1/(size*h) * sum(mat[4,])
    return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
  
  search.eq <- as.numeric(matrix(search.eq))
  
  if(grid.b[1] > 0){
    b.n <- grid.b[which(search.eq < 0)[1]]
    b.p <- grid.b[which(search.eq < 0)[1] -1]
  }else{
    b.p <- grid.b[which(search.eq > 0)[1]]
    #b.n <- grid.b[tail(which(search.eq < 0), 1)]
    b.n <- grid.b[which(search.eq > 0)[1] -1]
  }
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
}


b_quant <- function(z, s, n, psi, h, bpoint, b.end, b.start){
  
 
   size <- length(z)
   C <- z
  
  grid.b <- seq(b.start, b.end, length.out = bpoint)
  
  search.eq <- lapply(grid.b, FUN = b_grid_prep_quant, psi = psi, C = C, h = h)
  
  search.eq <- lapply(search.eq, function(x, size, h, n){
    mat <- matrix(x, nrow = 4, byrow = TRUE)
    tmp1 <- sum(mat[1, ])/size
    tmp2 <- sum(mat[2, ])/size
    tmp3 <- 1/(size*h) * sum(mat[3,])
    tmp4 <- 1/(size*h) * sum(mat[4,])
    return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
  
  search.eq <- as.numeric(matrix(search.eq))

  if(search.eq[1] > 0){
    b.n <- grid.b[which(search.eq < 0)[1]]
    b.p <- grid.b[which(search.eq < 0)[1] -1]
  }else{
    b.p <- grid.b[which(search.eq > 0)[1]]
    b.n <- grid.b[which(search.eq > 0)[1] -1]
  }
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
  #}
}

b_quant_double <- function(z, s, n, psi, h, bpoint, b.end, b.start){
  

   size <- length(z)
   C <- z
  
  grid.b <- seq(b.start, b.end, length.out = bpoint)
  
  search.eq <- lapply(grid.b, FUN = b_grid_prep_quant, psi = psi, C = C, h = h)
  
  search.eq <- lapply(search.eq, function(x, size, h, n){
    mat <- matrix(x, nrow = 4, byrow = TRUE)
    tmp1 <- sum(mat[1, ])/size
    tmp2 <- sum(mat[2, ])/size
    tmp3 <- 1/(size*h) * sum(mat[3,])
    tmp4 <- 1/(size*h) * sum(mat[4,])
    return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
  
  search.eq <- as.numeric(matrix(search.eq))

  if(search.eq[1] > 0){
    b.n <- grid.b[which(search.eq < 0)[1]]
    b.p <- grid.b[which(search.eq < 0)[1] -1]
  }else{
    b.p <- grid.b[which(search.eq > 0)[1]]
    b.n <- grid.b[which(search.eq > 0)[1] -1]
  }
  
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else{
    
    b.tmp <- (b.p + b.n)/2
    
    grid.b <- seq(max(0, b.tmp - 0.5), b.tmp + .5, length.out = bpoint)
    
    search.eq <- lapply(grid.b, FUN = b_grid_prep_quant, psi = psi, C = C, h = h)
    
    search.eq <- lapply(search.eq, function(x, size, h, n){
      mat <- matrix(x, nrow = 4, byrow = TRUE)
      tmp1 <- sum(mat[1, ])/size
      tmp2 <- sum(mat[2, ])/size
      tmp3 <- 1/(size*h) * sum(mat[3,])
      tmp4 <- 1/(size*h) * sum(mat[4,])
      return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
    
    search.eq <- as.numeric(matrix(search.eq))
   
    if(search.eq[1] > 0){
      b.n <- grid.b[which(search.eq < 0)[1]]
      b.p <- grid.b[which(search.eq < 0)[1] -1]
    }else{
      b.p <- grid.b[which(search.eq > 0)[1]]
      b.n <- grid.b[which(search.eq > 0)[1] -1]
    }
    
    if(length(b.p) == 0 | length(b.n) == 0) return(NA)
    else return((b.p + b.n)/2)
  }
}



b_abs_simu = function(W, Z, sigma, s, n, h, bpoint, b.start, b.end){
  
  C <- as.vector(sapply(W, function(x) x + sigma * Z))
  size <- length(W)*length(Z)
  
  grid.b <- seq(b.start, b.end, length.out = bpoint)
  
  search.eq <- lapply(grid.b, FUN = b_grid_prep_abs, C = C, h = h)
  
  search.eq <- lapply(search.eq, function(x, size, h, n){
    mat <- matrix(x, nrow = 4, byrow = TRUE)
    tmp1 <- sum(mat[1, ])/size
    tmp2 <- sum(mat[2, ])/size
    tmp3 <- 1/(size*h) * sum(mat[3,])
    tmp4 <- 1/(size*h) * sum(mat[4,])
    return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
  
  search.eq <- as.numeric(matrix(search.eq))
  
  if(grid.b[1] > 0){
    b.n <- grid.b[which(search.eq < 0)[1]]
    b.p <- grid.b[which(search.eq < 0)[1] -1]
  }else{
    b.p <- grid.b[which(search.eq > 0)[1]]
    b.n <- grid.b[which(search.eq > 0)[1] -1]
  }
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
}

b_quant_simu = function(W, Z, sigma, psi, s,n, h, bpoint, b.start, b.end){
  
  C <- as.vector(sapply(W, function(x) x + sigma * Z))
  size <- length(W)*length(Z)
  
  grid.b <- seq(b.start, b.end, length.out = bpoint)
  
  search.eq <- lapply(grid.b, FUN = b_grid_prep_quant, psi = psi, C = C, h = h)
  
  search.eq <- lapply(search.eq, function(x, size, h, n){
    mat <- matrix(x, nrow = 4, byrow = TRUE)
    tmp1 <- sum(mat[1, ])/size
    tmp2 <- sum(mat[2, ])/size
    tmp3 <- 1/(size*h) * sum(mat[3,])
    tmp4 <- 1/(size*h) * sum(mat[4,])
    return(tmp1 - tmp2 - tmp3 + tmp4 - s/n)}, size = size, h = h, n = n)
  
  search.eq <- as.numeric(matrix(search.eq))
  
  if(grid.b[1] > 0){
    b.n <- grid.b[which(search.eq < 0)[1]]
    b.p <- grid.b[which(search.eq < 0)[1] -1]
  }else{
    b.p <- grid.b[which(search.eq > 0)[1]]
    b.n <- grid.b[which(search.eq > 0)[1] -1]
  }
  if(length(b.p) == 0 | length(b.n) == 0) return(NA)
  else return((b.p + b.n)/2)
}



b_grid_prep_abs <- function(b, C, h){
  right <- b
  left <- -b
  
  F1 <- ifelse(C < right|C == right, 1, 0)
  F2 <- ifelse(C < left|C == left, 1, 0)
  
  
  f1 <- right * dnorm(C, mean = right, sd = h) ## estimate f(b)
  f2 <- left * dnorm(C, mean = left, sd = h) ## estimate f(-b)
  
  return(c(F1, F2, f1, f2))
  
  
}

b_grid_prep_quant <- function(b, psi, C, h){
  right <- b*psi
  left <- b*(psi - 1)
  
  F1 <- ifelse(C < right|C == right, 1, 0)
  F2 <- ifelse(C < left|C == left, 1, 0)
  
  
  f1 <- right * dnorm(C, mean = right, sd = h) ## estimate f(b)
  f2 <- left * dnorm(C, mean = left, sd = h) ## estimate f(-b)
  
  return(c(F1, F2, f1, f2))
  
  
}



sig = function(X,Z,tau,a,delta,...){
  
  temp <- sum(sapply(X, function(x, z, tau, a, delta) 1 / delta * (eta(x + tau * z,
                                                                           (tau * a)) - x)^2, 
                     z = Z, tau = tau, a = a, delta = delta))
  
  sig= (temp/length(X)/length(Z))
  return(sig)
}


tau_abs = function(W,Z,sig,b,omega,delta){
  temp = 0
  for(i in 1:length(W)){
    for(j in 1:length(Z))
      temp = temp + (delta/omega*Phi_abs(W[i]+sig*Z[j],b))^2
  }
    return((temp/length(W)/length(Z)))
} 


tau_quant = function(W,Z,sig,b,psi,omega,delta){
  temp = 0
  for(i in 1:length(W)){
    for(j in 1:length(Z))
      temp = temp + (delta/omega*Phi_quant(W[i]+sig*Z[j],b, psi))^2
  }
  return((temp/length(W)/length(Z)))
} 
#########################################################################
### here we consider the AMSE function
#########################################################################

amse <- function(tau2, alpha, x){
  
  tmp.gammas <- gammas(tau = sqrt(tau2), alpha, x)
  tmp <- tau2*tmp.gammas[1] + tmp.gammas[2]
  
  return(tmp)
}

##### gamma1 implementation
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


#####


##***********************************************************************
## RAMP search the optimal alpha -- golden search
## return the optimal alpha
##***********************************************************************

gold.RAMP.single <- function(A, y, x, tol, T, s, Z, W, psi, h,
                      bpoint = 500, b.start = 0, b.end = 5.5, loss.fun = loss.fun,
                      alpha.range, epsilon = 0.01){
  
  ## epsilon is the required final length of the interval which the minimizer stays in
  gr <- (sqrt(5) + 1) / 2
  
  left <- alpha.range[1]
  right <- alpha.range[2]
  tmp.left = right - (right - left) / gr
  tmp.right = left + (right - left) / gr 
  
  while(abs(right - left) > epsilon){
    
    x.est <- RAMP.single(A, y, a = tmp.left, s = s, tol = tol, T = T, psi = psi, 
                  W = W, Z = Z, x = x, bpoint = bpoint, b.start = b.start, 
                  b.end = b.end, loss.fun = loss.fun, h = h)
    #print(x.est)
    mse.left <- mean(tail(x.est$mse.est, n = 1L))
    x.est <- RAMP.single(A, y, a = tmp.right, s = s, tol = tol, T = T, psi = psi, 
                  W = W, Z = Z, x = x, bpoint = bpoint, b.start =b.start, 
                  b.end = b.end, loss.fun = loss.fun, h = h)
    mse.right <- mean(tail(x.est$mse.est, n = 1L))
    
    if(mse.left < mse.right) right <- tmp.right
    else left <- tmp.left
    
    tmp.left <- right - (right - left) / gr
    tmp.right <- left + (right - left) / gr
    
  }
  
  return((left + right)/2)
  
}

##
grid.RAMP.single <- function(A, y, x, tol, T, s, Z, W, psi, h,
                             bpoint = 500, b.start = 0, b.end = 5.5, loss.fun = loss.fun,
                             alpha.range, epsilon = 0.01){
  alpha.grid <- seq(alpha.range[1], alpha.range[2], epsilon)
  len.grid <- length(alpha.grid)
  
  mse.grid <- numeric(len.grid)
  for(alpha in alpha.grid){
    x.est <-RAMP.single(A, y, a = alpha, s = s, T = T, psi = psi, 
                        W = W, Z = Z, x = x, bpoint = bpoint, b.start = b.start, 
                        b.end = b.end, loss.fun = loss.fun, h = h)
    mse.grid[which(alpha == alpha.grid)] <- mean(tail(x.est$mse.est, n = 1L))
  }
  
  
  alpha.opt <- alpha.grid[which.min(mse.grid)]
  mse.min <- min(mse.grid)
  
  return(list(alpha.opt = alpha.opt, mse.min = mse.min, alpha.grid = alpha.grid, mse.grid = mse.grid))
}


