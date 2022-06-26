
#' ******************************************************************************
#' File was written to approximate the solution to composite quantile regression
#' ******************************************************************************



RAMP <- function(A, y, x, alpha, T, tol = 10^(-5), s, Z, W, psi, bpsi, nu, h = 5, 
                 b_start = 0 , b_end = 6.5, bpoint = 500, provide_true = FALSE){
  # xhat = reconstructAmp(A, y, T, tol, mse, verbose)
  #s
  #   Arguments:
  #       A - Design matrix
  #       y - response
  #       T - max number of iterations (optional)
  #       tol - stopping criteria (optional)
  #       x - original vector used to print progress of MSE (optional), true value
  #       sigma2 - the state evolution par using the estimated coefficient x from the previous step
  #       sigma2.theory - the state evolution par using independent r.v. Z and x
  #       tau2 - the state evolution par using the adjusted residuals
  #       tau2.theory - the state evolution par using independent r.v. Z and true noise W
  #       provide_true - if we provide the true sparsity, bpsi based on the empirical quantile of W;
  #                      if FALSE, we use the numerical optimization to obtain estimations of bpsi and s
  #                      and the true value of x is only used for calculation MSE's
  #       h - Bandwidth for the consistent estimator of f_z()
  #       bpoint - number of points in the grid for b_t
  
  
  # n, p, s define
  n <- dim(A)[1]
  p <- dim(A)[2]
  x0 <- rep(0, p)
  
  
  
  if(provide_true){
    if(missing(bpsi)) bpsi <- quantile(W, probs = psi)
    
    if(missing(s)) s <- length(which(x != 0))
  }
  else{
    if(missing(s) | missing(bpsi)){
        source("pwcqr.R")
        s.A <- scale(A, TRUE, FALSE)
        s.y <- drop(y- mean(y))
        pcqr.tmp <- penCQR(s.y, s.A, tau = tau, type = "cv")

        if(missing(s)) s <- length(which(pcqr.tmp$betas!=0))
       
        if(missing(bpsi)) bpsi <- pcqr.tmp$btau
    }
  }
  
  
  omega <- s/p
  delta <- n/p
  
  # Initial estimate, x0 is the coefficient and z0 is the residual
  #x0 <- rep(0, n)    
  z0 <- y - A %*% x0
  
  ##------------------------------------choose b0--------------------------------------##
  
  sigma02 <- sum(x0^2)/(p * delta)
  sigma02.theory <- mean(x0^2)/delta
  
  
  b0 <- b_quan_emp(z0, s, n, psi, bpsi, nu, h = h, b_start, b_end, bpoint)
  div <- 1
  while(is.na(b0)){
    b0 <- b_quan_emp(z0, s, n, psi, bpsi, nu, h = h*0.5/div, b_start, b_end, bpoint)
    div <- div + 1
    
    if(div >= 2) break
  }
  if(is.na(b0)){
      b0 <- b_quan_simu(W,Z,sqrt(sigma02.theory), s, n, psi, nu, bpsi, h, b_start, b_end, bpoint)
  }
  
  if(is.na(b0)) b0 <- 0.2
  

  
  
  ##------------------------------------------------------------------------------------------------------
 
  tau02 <- (n/s)^2*mean(Phii_quan(x0, b0, psi, bpsi, nu)^2)   #### p to s
  tau02.theory <- tau_quan(W,Z,sqrt(sigma02.theory),b0, w = omega, delta, psi, bpsi, nu)

  theta0 <- alpha*sqrt(tau02)
  arg_eta <- as.vector(x0 + t(A) %*% (n/s*Phii_quan(z0, b0, psi, bpsi, nu))) # p to s
  x.old <- x0
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
    z.new <- y - A %*% x.new + Phii_quan(z.old,bvec[t-1], psi, bpsi, nu)*mean(onsage.sca)/omega
  
    ##---------------------------------------------------------------------------------------------
   
    sigma2[t] <- sum((x.new - x)^2)/(delta*p)
   
    sigma2.theory[t] <- sig_quan(x, Z, sqrt(tau2.theory[t-1]),alpha = alpha, delta=delta)
    
    
    ##------------------------------------update b--------------------------------------##
    
    inc <- 0
    while(is.na(bvec[t])){
       bvec[t] <- b_quan_emp(as.vector(z.new), s, n, psi, nu, bpsi, h + inc*0.5, b_start, b_end, bpoint)
       inc <- inc + 1
       if(inc >= 2) break
      }
    
    if(is.na(bvec[t])){
        bvec[t] <- b_quan_simu(W, Z, sqrt(sigma2.theory[t]), s, n, psi, nu, bpsi, h, b_start, b_end, bpoint)
    }
    
    if(is.na(bvec[t])) bvec[t] <- bvec[t - 1]
    
    ##------------------------------------ estimation-----------------------------------##

    ## error W has 2nd order moment sigma^2
    
    arg_eta <- x.new + t(A) %*% (delta/omega*Phii_quan(z.new, bvec[t], psi, bpsi, nu))
    
    ##------------------------------------------------------------------------------------------------------
    
    tau2.theory[t] <- tau_quan(W,Z,sqrt(sigma2.theory[t]),bvec[t],w = omega, delta, psi, bpsi, nu)
    
   
    tau2.tmp <- (delta/omega)*Phii_quan(z.new, bvec[t], psi, bpsi, nu)
   
    tau2[t] <- mean(tau2.tmp^2)
   
    ##------------------------------------------------------------------------------------------
    
    theta[t] <- alpha*sqrt(tau2[t])
    
    x.old <- x.new
    x.new <- eta(arg_eta, theta = theta[t]) 
   
    
    mse[t] <- (sum((x.new - x)^2)/p)
    
    
    mse.est[t] <- tau2[t]*2*mean(etaprime(arg_eta, theta = theta[t])) - tau2[t] +
                        mean((x.new - arg_eta)^2)
    
    
    ## theoretical expression
    mse.est2[t] <- sum(sapply(x, function(x1, tau, theta) ((eta(x1 + sqrt(tau)*Z, theta) - x1)^2), 
                              tau = tau2[t], theta = theta[t]))/(p*n)
    mse.est3[t] <- amse(tau2[t], alpha, x)
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
  lambda <- alpha*sqrt(tau2[T_max])*omega/(bvec[T_max]*delta)
  
  return(list(x.new = x.new, bvec = bvec[1:T_max], omega.seq = omega.seq[1:T_max], sigma2 = sigma2[1:T_max], 
              tau2 = tau2[1:T_max], tau2.theory = tau2.theory[1:T_max], 
              sigma2.theory = sigma2.theory[1:T_max], theta= theta[1:T_max], mse = mse[1:T_max], 
              mse.est = mse.est[1:T_max], mse.est2 = mse.est2[1:T_max], 
              mse.est3 = mse.est3[1:T_max], lambda = lambda, T_max = T_max, tol.seq = tol.seq))
}


##*************************************************************************************************************
