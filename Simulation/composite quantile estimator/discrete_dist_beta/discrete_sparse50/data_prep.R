#' **************************************************************************************************
#' File was written to generate simulation data 
#' Design.mat is for generate design matrix where "std.amp" complies with assumption A1
#' beta is implemented for generating coefficient vectors but was not used in the final simulation version
#' **************************************************************************************************



Design.mat <- function(n, p, corr = c("std.amp", "correlated")){
  
  if(corr == "std.amp") A.vec <- rnorm(n*p, mean = 0, sd = sqrt(1/n))
  else{
    corr <- matrix(rep(0, p^2), c(p, p))
    for (i in seq(1, p, 1))
      for (j in seq(1, p, 1)) corr[i, j] <- .3^((abs(i-j)))
      
      A.vec <- as.vector(rmvnorm(n, mean = rep(0, p), corr))
  }
  
  A.vec <- A.vec - mean(A.vec)
  A.vec <- (A.vec/sd(A.vec))*(1/sqrt(n))
  A <- matrix(A.vec, c(n, p))
  
  return(A)
}


beta <- function(p = p, omega = omega, dist = c("norm", "bin")){
  ## s is prob of taking non-zero coefficients
  if(dist == "bin"){
    betatrue <- rep(0,p) 
    betatrue <- sample(c(-1, 1, 0), p, replace = TRUE, 
                       prob =c(omega/2, omega/2, 1 - omega))
  }else{
    betatrue <- c(rnorm(p * omega), rep(0, p*(1 - omega)))
  }
  #betatrue <- c(rep(1, p*0.064), rep(-1, p*.064), rep(0, p*.872))
  betatrue <- betatrue[sample(1:p, p, replace = FALSE)]
  
  return(betatrue)
}
