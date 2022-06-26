
#'******************************************************************************************************
#'   code here is used to first filter out the results where 
#'   any of the k'th component did not reach the determined convergence level
#'   
#'   then calculate the MSE of the zero parts, nonzero parts, and the whole vectors
#'******************************************************************************************************

#' read in total iterations 
#' change accordingly to the files names for different 

## random generated intercept

beta.original.all <- read.table("beta_truenorm.txt", sep = " ", fill = TRUE, 
                                header = FALSE, col.names = 1:501)

#'******************************************************************************************************
#' read in tolerance by iteration

iter.25 <- read.table("tol_seq_ramp0.25norm_maqr_.txt", sep = " ", fill = TRUE, 
                      header = FALSE, col.names = 1:52)
iter.5 <- read.table("tol_seq_ramp0.5norm_maqr_.txt", sep = " ", fill = TRUE, 
                     header = FALSE, col.names = 1:52)
iter.75 <- read.table("tol_seq_ramp0.75norm_maqr_.txt", sep = " ", fill = TRUE, 
                      header = FALSE, col.names = 1:52)



#' check convergence by tolerance
#' can change tolerance, here the tolerance level is 10^(-6)

check.converge.tmp <- apply(iter.25, 1, function(x) any(x <= 10^(-6), na.rm = TRUE))
includ.25 <- iter.25[check.converge.tmp, 1]
check.converge.tmp <- apply(iter.5, 1, function(x) any(x <= 10^(-6), na.rm = TRUE))
includ.5 <- iter.5[check.converge.tmp, 1]
check.converge.tmp <- apply(iter.75, 1, function(x) any(x <= 10^(-6), na.rm = TRUE))
includ.75 <- iter.75[check.converge.tmp, 1]



#'select converged 
includ <- includ.25[which(includ.25 %in% includ.5 & includ.25 %in% includ.75)]

includ <- includ[!duplicated(includ)]

##************************************************************************************************************




comb.est.mean.nonzero <- comb.est.mean.zero <- 
  comb.est.whole <- truepositive <- truenegative <- accuracy <- 
  matrix(rep(0, 4*500), c(500, 4)) 

seed.ramp <- as.numeric(read.table("modavrg_beta_rampnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.opt <- as.numeric(read.table("modavrg_beta_optnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.lasso <- as.numeric(read.table("lasso_betanorm_maqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.eq <- as.numeric(read.table("modavrg_beta_eqnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[, 1])
seed.beta <- as.numeric(beta.original.all[, 1])



for(s in includ[1:500]){
  
  beta.original <- as.numeric(beta.original.all[which(s == seed.beta), -1])
  
  non.zero.index <- which(beta.original != 0)
  zero.index <- which(beta.original == 0)
  beta.original.non.zero <- beta.original[non.zero.index]
  
  #beta.original <- as.numeric(read.table("wave_sound_6575.txt", sep = " ", header = FALSE))
  
  ## the whole vector
  
  
  sound.ramp.wave.coef.full <- as.numeric(
    read.table("modavrg_beta_rampnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.ramp), -1])
  sound.opt.wave.coef.full <- as.numeric(
    read.table("modavrg_beta_optnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.opt), -1])
  
  sound.lasso.wave.coef.full <- as.numeric(
    read.table("lasso_betanorm_maqr.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.lasso), -1])
  
  sound.eq.wave.coef.full <- as.numeric(
    read.table("modavrg_beta_eqnorm.txt", sep = " ", header = FALSE, col.names = 1:501, fill = TRUE)[which(s == seed.eq), -1])
  
  
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
  
  
  
  
  
  
  ## non-zero fraction of the vector
  beta.original.zero <- beta.original[zero.index]
  
  #beta.original <- as.numeric(read.table("wave_sound_6575.txt", sep = " ", header = FALSE))
  
  
  
  sound.ramp.wave.coef <- sound.ramp.wave.coef.full[zero.index]
  sound.opt.wave.coef <- sound.opt.wave.coef.full[zero.index]
  
  sound.lasso.wave.coef <- sound.lasso.wave.coef.full[zero.index]
  
  sound.eq.wave.coef <- sound.eq.wave.coef.full[zero.index]
  
  
  comb.est <- cbind((sound.ramp.wave.coef - beta.original.zero)^2, (sound.opt.wave.coef - beta.original.zero)^2,
                    (sound.eq.wave.coef - beta.original.zero)^2, (sound.lasso.wave.coef - beta.original.zero)^2)
  
  
  truenegative.tmp <- c(length(which(sound.ramp.wave.coef == 0)), 
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



#' final results reported in the tables 
#' Table 1
apply(as.matrix(comb.est.mean.zero), 2, mean, na.rm = TRUE)
apply(comb.est.mean.nonzero, 2, mean, na.rm = TRUE)
apply(comb.est.whole, 2, mean, na.rm = TRUE)

#'Table 2
apply(truepositive, 2, mean, na.rm = TRUE)
apply(truenegative, 2, mean, na.rm = TRUE)
apply(accuracy, 2, mean, na.rm = TRUE)


