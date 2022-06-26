source("LS.R")
source("pwcqr.R")
source("RAMP_main.R")
source("RAMP_search_tune.R")
source("supplementary_CQR_RAMP.R")
source("data_prep.R")
source("mixnormal.R")

#'****************************************************************************************
#'load the first for model averaged 
#'load the second for composite
source("main_body_of_the_simu.R")
#source("main_body_of_the_simu_cqr.R")


#'****************************************************************************************
#' Replace the file in the folders for simulation
#' 
#'****************************************************************************************
WrapFull <- function(s){
    
     set.seed(1)
   
     beta.original <- as.numeric(read.table("wave_sound_78.txt", sep = " ", header = FALSE))
    

     p <- length(beta.original)
     n <- floor(0.5*p)
     delta <- n/p
     A <- Design.mat(n = n, p = p, corr = "std.amp")
     
     set.seed(s)
     
     #' change noise distribution
     # noise <- rt(n, df = 3)
     # noise <- mixnormal(n = n, mean1 = 0, mean2 = 5, sd1 = 1, sd2 = 3, mix_prob = 0.5)
     noise <- rnorm(n)      

     noise <- noise - mean(noise)
     noise <- noise / sd(noise) * 0.03
     y <- A %*% beta.original + noise

     seed <- s

     # quantiles <- c(5, 7, 9)
     # for(k in quantiles){
     k <- 3
        
     tau <- (1 : k) / (k + 1)
        #tau <- 0.5
     test_coef_dist(seed, beta.original = beta.original, K = k, y = y, n = n, p = p, A = A, delta = delta, tau = tau, err.dist = "norm")
  
      
     #   }

}

