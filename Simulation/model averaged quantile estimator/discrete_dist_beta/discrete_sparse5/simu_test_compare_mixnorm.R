source("LS.R")
source("pwcqr.R")
source("RAMP_main.R")
source("RAMP_search_tune.R")
source("supplementary_CQR_RAMP.R")
source("data_prep.R")
source("mixnormal.R")
source("main_body_of_the_simu.R")

##**************************************************************************
## simulation setup
##**************************************************************************
WrapFull <- function(s){
  set.seed(s)
  
  beta.original <- c(sample(c(-1, 1), size = 5, replace = TRUE, prob = c(0.5, 0.5)), rep(0, 495))
  # beta.original <- c(rnorm(50), rep(0, 450))
  
  beta.original <- sample(beta.original, size = 500, replace = FALSE)
  
  #set.seed(s)
  
  #p <- dim(beta.original)[1]
  p <- length(beta.original)
  n <- floor(0.5*p)
  delta <- n/p
  A <- Design.mat(n = n, p = p, corr = "std.amp")
  
  # set.seed(s)
  #noise <- rt(n, df = 3)
  noise <- mixnormal(n = n, mean1 = 0, mean2 = 5, sd1 = 1, sd2 = 3, mix_prob = 0.5)
  #noise <- rnorm(n)      
  
  noise <- noise - mean(noise)
  noise <- noise / sd(noise) * 0.2
  y <- A %*% beta.original + noise
  #matrix(noise, c(n, p))
  
  seed <- s
  
  # quantiles <- c(5, 7, 9)
  #  for(k in quantiles){
  k <- 3
  
  tau <- (1 : k) / (k + 1)
  #tau <- 0.5
  test_coef_dist(seed, beta.original = beta.original, K = k, y = y, n = n, p = p, A = A, delta = delta, tau = tau, err.dist = "mixnorm")
  

}
