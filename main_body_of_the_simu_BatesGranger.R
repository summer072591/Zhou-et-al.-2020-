
#' ******************************************************************************
#' File contains the main simulation framework
#' 
#' ******************************************************************************

test_coef_dist <- function(seed, beta.original, K, y, n, p, A, delta, tau, err.dist){
      set.seed(seed)
  
      h.maqr <- 20
    
    
     if(missing(delta)) delta <- n / p
  
      
    #alpha.lower <- alpha_lower(delta)
     alpha.range <- c(1.3, 2.3)
      
     s.A <- scale(A, TRUE, FALSE)
      
      
      s.y <- drop(y- mean(y))
      beta0 <- pwcqr(s.y, s.A, 0.5)$betas
      r <- s.y - s.A %*% beta0
      
      if(K == 1) opt.weight.1.tmp <- 1
      else opt.weight.1.tmp <- quantile.weight.weighted.QP(r, tau)
      S <- maqr.betas <- maqr.btaus <- NULL
      
      
      for(l in 1:K){
        maqr.tmp <- penCQR(s.y, s.A, tau = tau[l], weight = 1, beta0 = beta0, type = "cv")
        
        maqr.betas <- cbind(maqr.betas, maqr.tmp$betas)
        write.table(t(c(seed, maqr.tmp$betas)),
              file = paste("betas_opt", tau[l], err.dist, ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)

        maqr.btaus <- cbind(maqr.btaus, maqr.tmp$btau)
        S <- c(S, length(which(maqr.tmp$betas != 0)))
      }
      
      if(K == 1) maqr.betas.1.tmp <- maqr.betas.eq.tmp <- maqr.betas
      else {
      maqr.betas.1.tmp <- maqr.betas %*% opt.weight.1.tmp

      }
       

      W <- y - A %*% maqr.betas.1.tmp
      maqr.betas.ramp <- maqr.argera <- alpha.opt <- NULL
      for(l in 1:K){
        c.y <- y - maqr.btaus[l]
        Z <- rnorm(n)
        alpha_opt <- try(gold.RAMP.single(A, c.y, x = beta.original, tol = 10^(-8), T = 50, s = S[l],
                                          Z = Z, W = W, psi = tau[l],
                                          h = h.maqr, bpoint = 1000, b.start = 0, b.end = 5, loss.fun = "quant",
                                          alpha.range = alpha.range, epsilon = 0.01))
        
        if(!inherits(alpha_opt, "try-error") == FALSE) alpha_opt <- tail(alpha.opt, n = 1L)
        if(length(alpha_opt) == 0) alpha_opt <- 1.8
        
        x.est <- try(RAMP.single(A, c.y, a = alpha_opt, s = S[l],
                                          T = 50, tol = 10^(-8), W = W, Z = Z,
                                          psi = tau[l], x = beta.original, h = h.maqr, bpoint = 1000,
                                          b.start = 0, b.end = 5, loss.fun = "quant"))
        
        
        maqr.betas.ramp <- cbind(maqr.betas.ramp, x.est$x.new)
        maqr.argera <- cbind(maqr.argera, x.est$arg_eta)
        alpha.opt <- c(alpha.opt, alpha_opt)

     write.table(t(c(seed, x.est$x.new)),
              file = paste("betas_ramp", tau[l], err.dist, ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)

     write.table(t(c(seed, x.est$tol.seq)),
              file = paste("tol_seq_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
     write.table(t(c(seed, x.est$bvec)),
              file = paste("bvec_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
     write.table(t(c(seed, x.est$T_max)),
              file = paste("max_iteration", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
     write.table(t(c(seed, x.est$mse.est)),
              file = paste("mse_est_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
     write.table(t(c(seed, x.est$mse)),
              file = paste("mse_emp_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)
     write.table(t(c(seed, x.est$mse.est3)),
              file = paste("mse_est3_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)

    write.table(t(c(seed, x.est$tau2)),
              file = paste("tau2_ramp", tau[l], err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)



      }
      
      
     ## the second type of optimal weight
     if(K == 1) {
         opt.weight.2.tmp <- 1
         maqr.betas.2.tmp <- maqr.betas.ramp
         maqr.betas.1.tmp <- maqr.betas.eq.tmp <- maqr.betas.ramp

     }
     else {
        if(all(!is.na(maqr.argera)) == FALSE)  maqr.betas.2.tmp <- maqr.betas.1.tmp
        else{
        
        Sigma <- diag(apply(maqr.argera, 2, var))
        
        ## the two versions give out very close results
        
        m_inv <- try(solve(Sigma))
        I_k <- rep(1, K)
        
        ## the optimal weight without restriction and use generalized inverse instead of the inverse
        
        
        w_opt <- solve.QP(Dmat = Sigma, dvec=rep(0,K),
                          Amat=cbind(matrix(1,ncol=1,nrow=K), diag(K)), 
                          bvec=c(1, rep(0,K)), meq=1)

        opt.weight.2.tmp <- w_opt$solution   
        opt.weight.2.tmp <- round(opt.weight.2.tmp/sum(opt.weight.2.tmp), 8)


 
        # opt.weight.2 <- rbind(opt.weight.2, opt.weight.2.tmp)
        maqr.betas.2.tmp <- maqr.betas.ramp %*% opt.weight.2.tmp

       maqr.betas.1.tmp <- maqr.betas.ramp %*% opt.weight.1.tmp
       maqr.betas.eq.tmp <- maqr.betas.ramp %*% rep(1/K, K)
              }
      }
    
     lasso.fit.beta <- cv.glmnet(A, y, alpha = 1, intercept = FALSE)
      
     write.table(t(c(seed, as.numeric(coef(lasso.fit.beta, lambda = lasso.fit.beta$lambda.min)[-1]))),
              file = paste("lasso_beta", err.dist, "_maqr", ".txt", sep = ""),
              append = TRUE, col.names = FALSE, row.names = FALSE)


   
     write.table(t(c(seed, alpha.opt)), 
              file = paste("modavrg_alpha_opt_ramp", err.dist, ".txt", sep = ""), append = TRUE, 
              col.names = FALSE, row.names = FALSE)
     
      write.table(t(c(seed, maqr.betas.1.tmp)), file = paste("modavrg_beta_opt", err.dist, ".txt", sep = ""), append = TRUE, 
              col.names = FALSE, row.names = FALSE)

      write.table(t(c(seed, maqr.betas.2.tmp)), file = paste("modavrg_beta_ramp", err.dist, ".txt", sep = ""), append = TRUE, 
              col.names = FALSE, row.names = FALSE)

      write.table(t(c(seed, maqr.betas.eq.tmp)), file = paste("modavrg_beta_eq", err.dist, ".txt", sep = ""), append = TRUE, 
              col.names = FALSE, row.names = FALSE)
              
      write.table(t(c(seed, beta.original)), file = paste("beta_true", err.dist, ".txt", sep = ""), append = TRUE,
              col.names = FALSE, row.names = FALSE)

      write.table(t(c(seed, opt.weight.1.tmp)), file = paste("modavrg_opt_weight_1", err.dist, ".txt", sep = ""), append = TRUE, 
                  col.names = FALSE, row.names = FALSE)

      write.table(t(c(seed, opt.weight.2.tmp)), file = paste("modavrg_opt_weight_2", err.dist, ".txt", sep = ""), append = TRUE, 
                  col.names = FALSE, row.names = FALSE)

  
    }
    
