##*************************************************************************************************************
## RAMP cross validation for searching the optimal alpha -- grid search
##*************************************************************************************************************

cv.RAMP <- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi,
                    nu, h, b_start = 0 , b_end = 5, bpoint = 200, provide_true = FALSE, 
                    alpha.grid, folds = 5){
  len.grid <- length(alpha.grid)
  lambda.grid <- numeric(len.grid)
  n <- dim(A)[1]
  p <- dim(A)[2]
  index <- sample(1:n, n, replace = FALSE)
  valid_size <- round(n/folds)  ## validation set size
  cv.mse <- matrix(rep(0, folds*len.grid), c(folds, len.grid))
  
  
  
  for(fold in 1:folds){
    valid_seq <- index[(1 + (fold - 1)*valid_size) : (fold*valid_size)]   ## validation indexes
    train.A <- A[-valid_seq, ]
    train.y <- y[-valid_seq]
    valid.A <- A[valid_seq, ]
    valid.y <- y[valid_seq]
    
    for(alpha in alpha.grid){
      x.est <- RAMP(train.A, train.y, x, alpha = alpha, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi,
                    nu, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true)
      
      cv.mse[fold, which(alpha == alpha.grid)] <- mean(tail(x.est$mse, n = 3L))
    }
  }
  cv.mse.mean <- apply(cv.mse, 2, mean)
  
  alpha.opt <- alpha.grid[which.min(cv.mse.mean)]
  
  return(list(alpha.opt = alpha.opt, cv.mse = cv.mse, cv.mse.mean = cv.mse.mean))
}


##***********************************************************************
## RAMP search the optimal alpha -- golden search
## return the optimal alpha
##***********************************************************************
#alpha.range <- c(1.3, 2.1)
gold.RAMP <- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi,
                      nu, h, b_start = 0 , b_end = 5, bpoint = 200, 
                      provide_true = FALSE, alpha.range, epsilon = 0.01){
  
  ## epsilon is the required final length of the interval which the minimizer stays in
  gr <- (sqrt(5) + 1) / 2
  
  left <- alpha.range[1]
  right <- alpha.range[2]
  tmp.left = right - (right - left) / gr
  tmp.right = left + (right - left) / gr 
  
  while(abs(right - left) > epsilon){
    
    x.est <- RAMP(A, y, x, alpha = tmp.left, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi,
                  nu, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint)
    mse.left <- mean(tail(x.est$mse, n = 3L))
    x.est <- RAMP(A, y, x, alpha = tmp.right, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi,
                  nu, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true)
    mse.right <- mean(tail(x.est$mse, n = 3L))
    
    if(mse.left < mse.right) right <- tmp.right
    else left <- tmp.left
    
    tmp.left <- right - (right - left) / gr
    tmp.right <- left + (right - left) / gr
    
   # print(c(left, right))
  }
  
  return((left + right)/2)
  
  
}


#gold.RAMP(A, y, x, tol = 0.05, T = 10, s = 64, Z = Z, W = W, psi = c(.1, .5, .7), bpsi,
#          nu = rep(1/3,3), h = 1.1, 
#          b_start = 0 , b_end = 5, bpoint = 200, alpha.range = alpha.range, epsilon = 0.1)

##***********************************************************************

##***********************************************************************
## RAMP search the optimal alpha -- grid search
## return the optimal alpha and the corresponding mse
##***********************************************************************
alpha.grid<- c(1.3, 2.1)
grid.RAMP <- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi,
                      nu, h, b_start = 0 , b_end = 5, bpoint = 200, 
                      provide_true = FALSE, alpha.grid){
  
  len.grid <- length(alpha.grid)
  
  mse.grid <- numeric(len.grid)
  for(alpha in alpha.grid){
    x.est <- RAMP(A, y, x, alpha = alpha, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi, 
                  nu, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true)
    
    mse.grid[which(alpha == alpha.grid)] <- mean(tail(x.est$mse, n = 3L))
  }
  
  
  alpha.opt <- alpha.grid[which.min(mse.grid)]
  mse.min <- min(mse.grid)
  
  return(list(alpha.opt = alpha.opt, mse.min = mse.min, alpha.grid = alpha.grid, mse.grid = mse.grid))
}






## lambda from the optimization method 
#optimize.lambda <- (pcqr.tmp$lambda)*320

#mse.optimize <- mean((pcqr.tmp$betas - x)^2)



##***********************************************************************
##the function used for searching the nearby points for the weight vector 
##***********************************************************************

nearby_search <- function(nu, rad, rad_points, centroid_keep = FALSE){
  ## nu is the weight vector
  ## rad is the distance restriction denoting how far it deviates from the current weight 
  ##                                              rad is entry-wise
  ## if rad is not specified, choose rad = smallest weight/3
  
  ## Create a Data Frame from All Combinations of the entry-wise nearby points
  if(missing(rad)) rad <- min(nu)/3
  
  K <- length(nu)
  var_nu <- nu[1:(K - 1)]
  fix_nu <- nu[K]
  rad_p <- seq(0 + rad/rad_points, rad, length.out = rad_points)
  
  var_left <- sapply(rad_p, function(x) var_nu - x)
  var_right <- sapply(rad_p, function(x) var_nu + x)
  
  mat.prep <- t(cbind(var_left, var_nu, var_right))
  row.names(mat.prep) <- NULL
  mat.tmp <- as.matrix(expand.grid(data.frame(mat.prep)))  ## the expanded combinations 
  weight_search <- cbind(mat.tmp, 1 - apply(mat.tmp, 1, function(x) sum(x))) 
  constraint <- apply(weight_search, 1, function(x) all(x >= 0 & x <= 1))
  weight_search <- weight_search[constraint, ]
  ## add the last column to the weight matrix
  
  center <- which(apply(weight_search, 1, function(x) sum(abs(x - nu)) < 10^(-6)))
  
  if(centroid_keep) weight_search <- rbind(center, weight_search[-center, ])  ## the center goes to the top 
  else weight_search <- weight_search[-center, ]
  
  return(weight_search)
}


##*********************************************************************
## Tabu search with golden section search
##*********************************************************************

Tabu_weight_search.gs <- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi,
                                  start.nu, h = 1.1, b_start = 0 , b_end = 5, bpoint = 200, provide_true = FALSE,
                                  alpha.range = alpha.range, 
                                  epsilon = 0.01, tabu.step = 5, rad = .1, rad_points = 3, centroid_keep = FALSE){
  
  centroid.history <- centroid <- start.nu
  alpha.opt <- gold.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi = psi, bpsi = bpsi,
                         nu = centroid, h = h, 
                         b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true,
                         alpha.range = alpha.range, epsilon = epsilon)
  tmp <- RAMP(A, y, x, alpha = alpha.opt, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi, 
              nu = centroid, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true)
  mse.centroid <- mean(tail(tmp$mse, n = 3L))
  tabu.points <- centroid 
  step <- 0
  
  while(step < tabu.step){
    
    surrounding <- nearby_search(nu = centroid, rad = rad, rad_points = rad_points, centroid_keep = FALSE)
    if(step > 0) repeat.check <- apply(surrounding, 1, function(x) all(colSums(abs(x - t(tabu.points))) > 10^(-6)))
    else repeat.check <- rep(TRUE, dim(surrounding)[1])
    
    
    if(length(repeat.check) == 0) break
    else{surrounding <- as.matrix(surrounding[repeat.check, ])}
    
    nearby.point.size <- dim(surrounding)[1]
    
    if(is.null(nearby.point.size) | nearby.point.size == 0) break
    else{
        tabu.points <- rbind(tabu.points, surrounding)
    
        alpha.opt <- mse.ramp <- numeric(nearby.point.size)
        for(nu_p in 1:nearby.point.size){
      
            alpha.opt[nu_p] <- gold.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi = psi, bpsi = bpsi,
                                     nu = surrounding[nu_p, ], h = h, 
                                     b_start = b_start , b_end = b_end, bpoint = bpoint, 
                                     provide_true = provide_true, alpha.range = alpha.range, epsilon = epsilon)
            tmp <- RAMP(A, y, x, alpha = alpha.opt[nu_p], T = T, tol = tol, s, Z = Z, W = W, psi = psi, bpsi = bpsi,
                     nu = surrounding[nu_p, ], h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, 
                     provide_true = provide_true)
      #s.true <- length(which(tmp$x.new!=0))
            mse.ramp[nu_p] <- mean(tail(tmp$mse, n = 3L))
          }
    
    ## optimal weight record
    centroid <- surrounding[which.min(mse.ramp), ]
    centroid.history <- rbind(centroid.history, centroid)
    mse.centroid <- c(mse.centroid, min(mse.ramp))
    step <- step + 1
    }
  }
  weight.opt <- centroid.history[which.min(mse.centroid),]
  mse.opt <- mse.centroid[which.min(mse.centroid)]
  
  #print(c(weight.opt, mse.opt))
  return(list(weight.opt = weight.opt, mse.opt = mse.opt, centroid.history = centroid.history, 
              mse.centroid = mse.centroid))
}

#time.tabu.start <- Sys.time()
#test_tabu <- Tabu_weight_search.gs(s.A, s.y, tol = 10^(-6), T = 5, s = 64, Z = Z, 
#                                   psi = c(.1, .5, .7), start.nu = opt.weight.optmz, h = 1.1, 
#                                   b_start = 0 , b_end = 5, bpoint = 200, alpha.range = alpha.range, 
#                                   epsilon = 0.3, tabu.step = 5, rad = .1, rad_points = 2, centroid_keep = FALSE)
#time.tabu.end <- Sys.time()




##*********************************************************************
## Tabu search with grid search
##*********************************************************************

Tabu_weight_search.grid<- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi,
                                   start.nu, h, b_start = 0 , b_end = 5, bpoint = 200, 
                                   provide_true = FALSE, alpha.grid, 
                                   tabu.step = 5, rad = .1, rad_points = 3, centroid_keep = FALSE){
  
  centroid.history <- centroid <- start.nu
  grid.search <- grid.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi = psi, bpsi = bpsi, 
                           nu = centroid, h = h, 
                           b_start = b_start , b_end = b_end, bpoint = bpoint, provide_true = provide_true,
                           alpha.grid = alpha.grid)
  
  mse.centroid <- grid.search$mse.min
  tabu.points <- centroid 
  step <- 0
  
  while(step < tabu.step){
    
    surrounding <- nearby_search(nu = centroid, rad = rad, rad_points = rad_points, centroid_keep = FALSE)
    if(step > 0) repeat.check <- apply(surrounding, 1, function(x) all(colSums(abs(x - t(tabu.points))) > 10^(-6)))
    else repeat.check <- rep(TRUE, dim(surrounding)[1])
    
    if(length(repeat.check) == 0) break
    else{surrounding <- surrounding[repeat.check, ]}
    
    #print(surrounding)
    nearby.point.size <- dim(surrounding)[1]
    
    if(is.null(nearby.point.size) | nearby.point.size == 0) break
    else{
        tabu.points <- rbind(tabu.points, surrounding)
        alpha.opt <- mse.ramp <- numeric(nearby.point.size)
        for(nu_p in 1:nearby.point.size){
      
           tmp <- grid.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi = psi, bpsi = bpsi, 
                          nu = surrounding[nu_p, ], h = h, 
                          b_start = b_start, b_end = b_end, bpoint = bpoint, provide_true = provide_true,
                          alpha.grid = alpha.grid)
      
      #s.true <- length(which(tmp$x.new!=0))
      
          mse.ramp[nu_p] <- tmp$mse.min
            }
    
    ## optimal weight record
       centroid <- surrounding[which.min(mse.ramp), ]
       centroid.history <- rbind(centroid.history, centroid)
       mse.centroid <- c(mse.centroid, min(mse.ramp))
       step <- step + 1
        }
  }
  weight.opt <- centroid.history[which.min(mse.centroid),]
  mse.opt <- mse.centroid[which.min(mse.centroid)]
  return(list(weight.opt = weight.opt, mse.opt = mse.opt, centroid.history = centroid.history, 
              mse.centroid = mse.centroid))
}

##******************************************************************************************
## fix the value of alpha and search only for the weight
##******************************************************************************************
Tabu_weight_search <- function(A, y, x, tol = 0.05, T, s, Z, W, psi, bpsi, 
                               start.nu, h = 1.1, b_start = 0 , b_end = 5, bpoint = 200, 
                               provide_true = FALSE, alpha_start, alpha_end, method = c("grid", "gs"),
                               tabu.step = 5, rad = .1, rad_points = 3, centroid_keep = FALSE, ...){
  
  centroid.history <- centroid <- start.nu
  
  if(method == "grid"){
     alpha.grid <- seq(alpha_start, alpha_end, len = 0.01)
     grid.search <- grid.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi = psi, bpsi = bpsi,
                              nu = centroid, h = h, 
                              b_start = b_start, b_end = b_end, bpoint = bpoint, 
                              provide_true = provide_true, alpha.grid = alpha.grid)
     alpha.opt <- grid.search$alpha.opt
     mse.centroid <- grid.search$mse.min
  }else{
     alpha.range <- c(alpha_start, alpha_end)
     alpha.opt <- gold.RAMP(A, y, x, tol = tol, T = T, s = s, Z = Z, W = W, psi, bpsi,
                              nu = centroid, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, 
                              provide_true = FALSE, alpha.range = alpha.range, epsilon = 0.01)
      tmp <- RAMP(A, y, x, alpha = alpha.opt, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi, 
                           nu = centroid, h = h, b_start = b_start , b_end = b_end, bpoint = bpoint, 
                           provide_true = provide_true)
      mse.centroid <- mean(tail(tmp$mse, n = 3L))
  }
  #tmp <- RAMP(A, y, alpha = 2, T = T, tol = tol, s, Z = Z, psi, 
  #            nu = centroid, h = h, b_start = b_start , 
  #            b_end = b_end, bpoint = bpoint)
  
  #mse.centroid <- mean(tail(tmp$mse, n = 3L))
  
  tabu.points <- centroid 
  step <- 0
  
  #print(tabu.points)
  while(step < tabu.step){
    
    surrounding <- nearby_search(nu = centroid, rad = rad, rad_points = rad_points, centroid_keep = FALSE)
    if(step > 0) repeat.check <- apply(surrounding, 1, function(x) all(colSums(abs(x - t(tabu.points))) > 10^(-6)))
    else repeat.check <- rep(TRUE, dim(surrounding)[1])
    
    
    if(length(repeat.check) == 0) break
    else surrounding <- as.matrix(surrounding[repeat.check, ])
    
    #print(surrounding)
    tabu.points <- rbind(tabu.points, surrounding)
    nearby.point.size <- dim(surrounding)[1]
    mse.ramp <- numeric(nearby.point.size)
    for(nu_p in 1:nearby.point.size){
      
      tmp <- RAMP(A, y, x, alpha = alpha.opt, T = T, tol = tol, s, Z = Z, W = W, psi, bpsi = bpsi, 
                  nu = surrounding[nu_p, ], h = h, b_start = b_start , 
                  b_end = b_end, bpoint = bpoint, provide_true = provide_true)
      
      #s.true <- length(which(tmp$x.new!=0))
      
      mse.ramp[nu_p] <- mean(tail(tmp$mse, n = 3L))
    }
    
    ## optimal weight record
    centroid <- surrounding[which.min(mse.ramp), ]
    centroid.history <- rbind(centroid.history, centroid)
    mse.centroid <- c(mse.centroid, min(mse.ramp))
    step <- step + 1
    
  }
  weight.opt <- centroid.history[which.min(mse.centroid),]
  mse.opt <- mse.centroid[which.min(mse.centroid)]
  return(list(alpha.opt = alpha.opt, weight.opt = weight.opt, mse.opt = mse.opt, centroid.history = centroid.history, 
              mse.centroid = mse.centroid))
}

##************************************************************************************************************
