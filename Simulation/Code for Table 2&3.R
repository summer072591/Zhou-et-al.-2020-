
#'******************************************************************************************************
#'   code here is used to first filter out the results where 
#'   any of the k'th component did not reach the determined convergence level
#'   
#'   then calculate the MSE of the zero parts, nonzero parts, and the whole vectors
#'******************************************************************************************************


beta.original.all <- read.table("beta_truenorm.txt", sep = " ", fill = TRUE, 
                                header = FALSE, col.names = 1:501)


#' convergence check
iter <- read.table("tol_emp_ramp3norm_cqr_.txt", sep = " ", fill = TRUE, 
                   header = FALSE, col.names = 1:52)


#' can change the tolerance
check.converge.tmp <- apply(iter, 1, function(x) any(x <= 10^(-4.4), na.rm = TRUE))

includ <- iter[check.converge.tmp, 1]

includ <- includ[!duplicated(includ)]
##************************************************************************************************************




comb.est.mean.nonzero <- comb.est.mean.zero <- 
  comb.est.whole <- truepositive <- truenegative <- accuracy <- 
  matrix(rep(0, 4*500), c(500, 4)) 

seed.ramp <- as.numeric(read.table("coef_est_ramp3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.opt <- as.numeric(read.table("coef_est_opt3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.lasso <- as.numeric(read.table("lasso_beta3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.eq <- as.numeric(read.table("coef_est_eq3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.beta <- as.numeric(beta.original.all[, 1])




for(s in includ[1:500]){
  
  beta.original <- as.numeric(beta.original.all[which(s == seed.beta), -1][1, ])
  
  non.zero.index <- which(beta.original != 0)
  zero.index <- which(beta.original == 0)
  beta.original.non.zero <- beta.original[non.zero.index]
  

  
  
  sound.ramp.wave.coef.full <- as.numeric(
    read.table("coef_est_ramp3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.ramp), -1][1, ])
  sound.opt.wave.coef.full <- as.numeric(
    read.table("coef_est_opt3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.opt), -1][1, ])
  
  sound.lasso.wave.coef.full <- as.numeric(
    read.table("lasso_beta3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.lasso), -1][1, ])
  
  sound.eq.wave.coef.full <- as.numeric(
    read.table("coef_est_eq3norm_cqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.eq), -1][1, ])
  
  
  ## non-zero 
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full[non.zero.index]
  sound.opt.wave.coef <- sound.opt.wave.coef.full[non.zero.index]
  
  sound.lasso.wave.coef <- sound.lasso.wave.coef.full[non.zero.index]
  
  sound.eq.wave.coef <- sound.eq.wave.coef.full[non.zero.index]
  
  
  
  comb.est <- cbind((sound.ramp.wave.coef - beta.original.non.zero)^2, 
                    (sound.opt.wave.coef - beta.original.non.zero)^2,
                    (sound.eq.wave.coef - beta.original.non.zero)^2, 
                    (sound.lasso.wave.coef - beta.original.non.zero)^2)
  truepositive.tmp <- c(length(which(sound.ramp.wave.coef !=0)), 
                        length(which(sound.opt.wave.coef != 0)),
                        length(which(sound.eq.wave.coef != 0)), 
                        length(which(sound.lasso.wave.coef != 0))) 
  
  comb.est.mean.nonzero[which(s == includ), ] <- apply(comb.est, 2, mean)
  truepositive[which(s == includ), ] <- truepositive.tmp / length(beta.original.non.zero)
  
  
  
  
  
  
  
  
  ## all-zero fraction of the vector
  beta.original.zero <- beta.original[zero.index]
  
  
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full[zero.index]
  sound.opt.wave.coef <- sound.opt.wave.coef.full[zero.index]
  
  sound.lasso.wave.coef <- sound.lasso.wave.coef.full[zero.index]
  
  sound.eq.wave.coef <- sound.eq.wave.coef.full[zero.index]
  
  comb.est <- cbind((sound.ramp.wave.coef - beta.original.zero)^2, (sound.opt.wave.coef - beta.original.zero)^2,
                    (sound.eq.wave.coef - beta.original.zero)^2, (sound.lasso.wave.coef - beta.original.zero)^2)
  
  
  truenegative.tmp <- c(length(which(sound.ramp.wave.coef ==0)), 
                        length(which(sound.opt.wave.coef == 0)),
                        length(which(sound.eq.wave.coef == 0)), 
                        length(which(sound.lasso.wave.coef == 0))) 
  
  
  comb.est.mean.zero[which(s == includ), ] <- apply(comb.est, 2, mean)
  truenegative[which(s == includ), ] <- truenegative.tmp / length(beta.original.zero)
  accuracy[which(s == includ), ] <- (truepositive.tmp + truenegative.tmp) / length(beta.original)
  
  
  ## the whole vector
  
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full
  sound.opt.wave.coef <- sound.opt.wave.coef.full
  
  sound.lasso.wave.coef <- sound.lasso.wave.coef.full
  
  sound.eq.wave.coef <- sound.eq.wave.coef.full
  
  comb.est <- cbind((sound.ramp.wave.coef - beta.original)^2, (sound.opt.wave.coef - beta.original)^2,
                    (sound.eq.wave.coef - beta.original)^2, (sound.lasso.wave.coef - beta.original)^2)
  
  
  comb.est.whole[which(s == includ), ] <- apply(comb.est, 2, mean)
  
  
}



#'Table 1
apply(comb.est.mean.zero, 2, mean, na.rm = TRUE)

apply(comb.est.mean.nonzero, 2, mean, na.rm = TRUE)

apply(comb.est.whole, 2, mean, na.rm = TRUE)

#'Table 2
apply(truepositive, 2, mean, na.rm = TRUE)
apply(truenegative, 2, mean, na.rm = TRUE)
apply(accuracy, 2, mean, na.rm = TRUE)


##************************************************************************************************************

#' the following code is written for first revision of Zhou et al. (2019) submitted to EJS
#' in this file, we want to compare performance of single quantile estimator
#' to composite quantile estimator
#' outputs correpond to the values in Table 2 & 3 for single quantile estimator \beta_{0.5}

##************************************************************************************************************

#' for composite estimator, we consider the first 500 converged
#' 
#' convergence check
iter <- read.table("tol_emp_ramp3norm_cqr_.txt", sep = " ", fill = TRUE, header = FALSE, col.names = 1:52)


#' can change the tolerance
check.converge.tmp <- apply(iter, 1, function(x) any(x <= 10^(-4.4), na.rm = TRUE))

includ <- iter[check.converge.tmp, 1]

includ <- includ[!duplicated(includ)]
##************************************************************************************************************



#' choose single quantile estimator using the same seed numbers 

beta.original.all <- read.table("beta_truenorm.txt", sep = " ", fill = TRUE, 
                                header = FALSE, col.names = 1:501)

iter.5 <- read.table("tol_seq_ramp0.5norm_maqr_.txt", sep = " ", fill = TRUE, 
                     header = FALSE, col.names = 1:52)


check.converge.tmp <- apply(iter.5, 1, function(x) any(x <= 10^(-6), na.rm = TRUE))
includ.5 <- iter.5[check.converge.tmp, 1]

seed.5 <- as.numeric(read.table("betas_ramp0.5norm.txt", sep = " ", 
                                header = FALSE, col.names = 1:501, fill = TRUE)[, 1])




comb.est.mean.nonzero <- comb.est.mean.zero <- 
  comb.est.whole <- truepositive <- truenegative <- accuracy <- numeric(500)

seed.beta <- as.numeric(beta.original.all[, 1])


for(s in includ[1:500]){
  
  beta.original <- as.numeric(beta.original.all[which(s == seed.beta), -1][1, ])
  
  non.zero.index <- which(beta.original != 0)
  zero.index <- which(beta.original == 0)
  beta.original.non.zero <- beta.original[non.zero.index]
  
  
  
  
  sound.ramp.wave.coef.full <- as.numeric(
    read.table("betas_ramp0.5norm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.5), -1][1, ])
  
  
  ## non-zero 
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full[non.zero.index]
  
  
  
  comb.est <- (sound.ramp.wave.coef - beta.original.non.zero)^2
  
  truepositive.tmp <- length(which(sound.ramp.wave.coef !=0))
  
  comb.est.mean.nonzero[which(s == includ)] <- mean(comb.est)
  
  truepositive[which(s == includ)] <- truepositive.tmp / length(beta.original.non.zero)
  
  
  
  
  
  
  
  
  ## all-zero fraction of the vector
  beta.original.zero <- beta.original[zero.index]
  
  
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full[zero.index]
  
  
  comb.est <- (sound.ramp.wave.coef - beta.original.zero)^2
  
  truenegative.tmp <- length(which(sound.ramp.wave.coef ==0))
  
  comb.est.mean.zero[which(s == includ)] <- mean(comb.est)
  truenegative[which(s == includ)] <- truenegative.tmp / length(beta.original.zero)
  accuracy[which(s == includ)] <- (truepositive.tmp + truenegative.tmp) / length(beta.original)
  
  
  ## the whole vector
  
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full
  
  comb.est <- (sound.ramp.wave.coef - beta.original)^2
  
  comb.est.whole[which(s == includ)] <- mean(comb.est)
  
  
}



#'Table 1
mean(comb.est.mean.zero, na.rm = TRUE)

mean(comb.est.mean.nonzero, na.rm = TRUE)

mean(comb.est.whole, 2, na.rm = TRUE)

#'Table 2
mean(truepositive, na.rm = TRUE)
mean(truenegative, na.rm = TRUE)
mean(accuracy, na.rm = TRUE)


