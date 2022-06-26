
#' ******************************************************************************
#' File contains the main simulation framework
#' 
#' ******************************************************************************



test_coef_dist <- function(seed, K, beta.original, y, n, p, A, delta, tau, err.dist){

   s.A <- scale(A, TRUE, FALSE)
   s.y <- drop(y- mean(y))
  # alpha.lower <- alpha_lower(delta)

   alpha.range <- c(1.3, 2.3)
      

   h.cqr <- 20
 
  
  
  #' ******************************************************************************
  #'                           composite quantile 
  #' ******************************************************************************
  
  #t1 <- Sys.time()

   pcqr.tmp <- penCQR(s.y, s.A, tau = tau, type = "cv")
  
  # t2 <- Sys.time()
   pcqr.tmp.eq <- penCQR(s.y, s.A, tau = tau, weight = rep(1/K, K), type = "cv")
  # t3 <- Sys.time()

  s.est.opt <- length(which(pcqr.tmp$betas!=0))
  cqr.betas.1.tmp <- pcqr.tmp$betas
  cqr.betas.eq <- pcqr.tmp.eq$betas
  #cqr.betas.1 <- rbind(cqr.betas.1, pcqr.tmp$betas)
  opt.weight.optmz <- pcqr.tmp$weight
  pcqr.btaus <- pcqr.tmp$btau
  pcqr.btaus.eq <- pcqr.tmp.eq$btau
  Z <- rnorm(n)
  
  
  W <- y - A %*% cqr.betas.1.tmp



   #test_tabu <- Tabu_weight_search(A, y, x, T = 10, s = 64, Z = Z, W = W,
  #                                   psi = tau, bpsi = pcqr.btaus, start.nu = opt.weight.optmz, h = 2.5, 
  #                                   b_start = 0 , b_end = 5, bpoint = 200, alpha_start = alpha.lower, 
  #                                   alpha_end = 2.1, method = "gs",
  #                                  rad = .1, rad_points = 2, centroid_keep = FALSE,
  #                                   tabu.step = 8)
  #time.tabu.end <- Sys.time()
  ## ********************* starting from the optimal weight **********************************************
   test_tabu_est <- Tabu_weight_search(A, y, x =  beta.original, T = 50, tol = 10^(-6),
                                    s = s.est.opt, Z = Z, W = W,
                                    psi = tau, bpsi = pcqr.btaus,
                                    start.nu =opt.weight.optmz,
                                    h = h.cqr,
                                    b_start = 0 , b_end = 7, bpoint = 4000,
                                    alpha_start = alpha.range[1],
                                    alpha_end = alpha.range[2],
                                    rad = .1, rad_points = 3, centroid_keep = FALSE,
                                    method = "gs",
                                    tabu.step = 5)

   x.est.cqr <- RAMP(A, y, x =  beta.original, alpha = test_tabu_est$alpha.opt, T = 50, tol = 10^(-6),
                                    s = s.est.opt,
                                    Z = Z, W = W, psi = tau, bpsi = pcqr.btaus,
                                    nu = test_tabu_est$weight.opt,
                                    h = h.cqr, b_start = 0 , b_end = 3, bpoint = 2000,
                                    provide_true = FALSE)
   cqr.betas.2.tmp <- x.est.cqr$x.new


   x.est.cqr <- RAMP(A, y, x =  beta.original, alpha = test_tabu_est$alpha.opt, T = 50, tol = 10^(-6),
                                    s = s.est.opt,
                                    Z = Z, W = W, psi = tau, bpsi = pcqr.btaus,
                                    nu = opt.weight.optmz,
                                    h = h.cqr, b_start = 0 , b_end = 3, bpoint = 2000,
                                    provide_true = FALSE)

   cqr.betas.1.tmp <- x.est.cqr$x.new

   ## ********************** starting from the equal weight *********************************************

   x.est.cqr <- RAMP(A, y, x =  beta.original, alpha = test_tabu_est$alpha.opt, T = 50, tol = 10^(-6),
                                    s = s.est.opt,
                                    Z = Z, W = W, psi = tau, bpsi = pcqr.btaus.eq,
                                    nu = rep(1/K, K),
                                    h = h.cqr, b_start = 0 , b_end = 3, bpoint = 2000,
                                    provide_true = FALSE)

 
   cqr.betas.eq <- x.est.cqr$x.new
 
   lasso.fit.beta <- cv.glmnet(A, y, alpha = 1, intercept = FALSE)
      

      
   write.table(t(c(seed, as.numeric(coef(lasso.fit.beta, lambda = lasso.fit.beta$lambda.min)[-1]))),
              file = paste("lasso_beta", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)



  ##*************************************************************************************************************

  

  write.table(t(c(seed, test_tabu_est$alpha.opt)), 
              file = paste("opt_alpha_ramp", K, err.dist, "_cqr", ".txt", sep = ""), 
              append = TRUE, col.names = FALSE, row.names = FALSE)


  write.table(t(c(seed, opt.weight.optmz)),
              file = paste("opt_weight_optm", K, err.dist, "_cqr", ".txt", sep = ""), 
              append = TRUE, col.names = FALSE, row.names = FALSE)

  write.table(t(c(seed, test_tabu_est$weight.opt)),
              file = paste("opt_weight_ramp", K, err.dist, "_cqr", ".txt", sep = ""), 
              append = TRUE, col.names = FALSE, row.names = FALSE)
    
  write.table(t(c(seed, x.est.cqr$bvec)),
              file = paste("bvec_ramp", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(seed, x.est.cqr$T_max)),
              file = paste("max_iteration", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
 
  write.table(t(c(seed, x.est.cqr$mse.est)),
              file = paste("mse_est_ramp", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(seed, x.est.cqr$mse)),
              file = paste("mse_emp_ramp", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)

  write.table(t(c(seed, x.est.cqr$tol.seq)),
              file = paste("tol_emp_ramp", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
              
            


  write.table(t(c(seed, mean((cqr.betas.1.tmp - cqr.betas.2.tmp)^2),
                  mean((cqr.betas.1.tmp - cqr.betas.eq)^2), mean((cqr.betas.2.tmp - cqr.betas.eq)^2))),
              file = paste("mse_dif_", K, err.dist, "_cqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)    
  
  write.table(t(c(seed, cqr.betas.1.tmp)), 
             file = paste("coef_est_opt", K, err.dist, "_cqr", ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
  write.table(t(c(seed, cqr.betas.2.tmp)), 
             file = paste("coef_est_ramp", K, err.dist, "_cqr", ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
  write.table(t(c(seed, cqr.betas.eq)),
             file = paste("coef_est_eq", K, err.dist, "_cqr", ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)

  write.table(t(c(seed, beta.original)), file = paste("beta_true", err.dist, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)



  
  
}
