# File: pwcqr.R
# penalized weighted CQR
# Author: Jelena Bradic and Weiwei Wang

# type: how to select lambda
# if using CV, beta0 should be the default value ( i.e. LAD)
penCQR <- function(y,x,tau,weight,lambda,beta0,type=c("cv","sic"),nfold=5,s0=0,s1=0.0025)
{      ## s0, s1 are the starting and ending points of the tuning paramter sequence
  call <- match.call();
  
  n <- nrow(x);
  p <- ncol(x);
  K <- length(tau);

  if(missing(beta0)){beta0 <- pwcqr(y,x,0.5)$betas}
  r <- y-x %*% beta0;

  if(missing(weight))
  {
    weight <- quantile.weight(r,tau);
  }

  #print("OK")
  
  # choose lambda
  if(missing(lambda))
  {
    vtmp <- c(0,1);
    while((s1-s0>1e-3) & (max(vtmp)-min(vtmp)>1e-2))
      {
        tmp <- seq(s0,s1,(s1-s0)/4);
        
        if(type=="cv"){
          vtmp <- cv.pwcqr(y,x,tau,weight,beta0,tmp,nfold,FALSE)$cv;
        }else if(type=="sic"){
          vtmp <- sic.pwcqr(y,x,tau,weight,beta0,tmp)$sic;
        }
      #  print(tmp)
       # print(vtmp)
        pos <- which.min(vtmp);
        if(pos==1){s1 <- tmp[2];}
        else if(pos==5){s0 <- tmp[4];}
        else {s0=tmp[pos-1];s1 <- tmp[pos+1];}
      }
    lambda <- s0;
  }

  # estimation
  tmp <- pwcqr(y,x,tau,weight,n*scad(lambda, beta0));
  c(tmp,list(lambda=lambda));
}


check.ftn <- function(y,x,beta,btau,tau,weight)
{
  K <- length(tau);
  n <- length(y);
  r <- y-x%*%beta;
  A <- matrix(r, ncol=K,nrow=n) - matrix(btau,ncol=K,nrow=n,byrow=TRUE);
  mean(A%*%(weight*tau) - (A*(A<0))%*%weight);
}


cv.pwcqr <- function(y,x,tau,weight,beta0,s,K=5,plot=TRUE)
{
  n <- length(y);
  p <- ncol(x);
  
  all.folds <-  split(sample(1:n), rep(1:K, length = n));
  residmat <- matrix(0, length(s), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    for (j in 1:length(s)){
      tmp <- pwcqr(y[-omit],x[-omit,],tau=tau, weight=weight,
                   penalty=scad(s[j],beta0)*length(y[-omit]))
      beta.fit <- tmp$betas;
      btau.fit <- tmp$btau;
      residmat[j, i] <- check.ftn(y[omit],x[omit,],beta.fit,btau.fit,tau,weight);
    }
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  if(plot){plot(s,cv)}
  object <- list(s = s, cv = cv, cv.error = cv.error)
}


sic.pwcqr <- function(y,x,tau,weight,beta0,s)
{
  K <- length(tau);
  n <- length(y);
  sic <- rep(0,length(s));
  df <- sic;
  for(i in 1:length(s)){
    fit <- pwcqr(y,x,tau,weight,n*scad(s[i],beta0));
    
    r <- y-x%*%fit$beta;
    elbowset <- matrix(r,nrow=n,ncol=K) - matrix(fit$btau,nrow=n,ncol=K,byrow=T);
    df[i] <- sum(abs(elbowset)<1e-6);    
    sic[i] <- log(check.ftn(y,x,fit$beta,fit$btau,tau,weight)/K)+log(n*K)/2/n/K*df[i];
  }
  list(s=s,sic=sic,df=df);
}



#
# arg min sum_{i=1}^n sum_{k=1}^K w_k rho_k(y_i-X_i*beta-b_k)
#         + sum_{j=1}^p penalty_j |beta_j|
#
pwcqr <- function(y,x,tau,weight=rep(1,length(tau)),penalty=rep(0,ncol(x)))
{
  call <- match.call();
  
  w <- penalty;
    
    require("lpSolve")   ## Linear and Integer Programming
    call <- match.call()
    K<-length(tau)
    p<-dim(x)[2] 
    n<-length(y)  
    cvec<-c(rep(c(0,0,rep(1,n)),K),w,w)/(K*n)
    cvec <- cvec * c(rep(weight,rep(n+2, K)), rep(1,2*p));
    bvec<-c(tau[1]*y,(tau[1]-1)*y) 
    A1<-cbind(rep(1,n),rep(-1,n))
    A2<-rbind(diag(n),diag(n))
    Aleft<-cbind(rbind(tau[1]*A1,(tau[1]-1)*A1),A2)   
    Aright<-rbind(tau[1]*cbind(x,-x),(tau[1]-1)*cbind(x,-x))     
    if(K>=2){
      for(j in 2:K){
        Aleft<-expandmatrix(Aleft,cbind(rbind(tau[j]*A1,(tau[j]-1)*A1),A2) )
        Aright<-rbind(Aright,rbind(tau[j]*cbind(x,-x),(tau[j]-1)*cbind(x,-x)))
        bvec<-c(bvec,c(tau[j]*y,(tau[j]-1)*y) )
      }
    }
    Amat<-cbind(Aleft,Aright)
    lp.out<-lp(direction="min",objective.in=cvec,const.mat=Amat,const.dir=rep(">=",length(bvec)),const.rhs=bvec)  
    Beta<-lp.out$solution
    btau<-tau
    for (j in 1:K){
      a6<-1+(j-1)*(2+n)
      a7<-j*(2+n)
      bbb<-Beta[(a6:a7)]
      btau[j]<-bbb[1]-bbb[2]
    }
    betas1<-rev(rev(Beta)[1:(2*p)]) 
    betas<-betas1[(1:p)]-betas1[((1+p):(2*p))]
    
    object <- list(call = call, tau=tau, penalty=penalty,  betas=betas, btau=btau, status=lp.out$status, weight=weight)
    object
  }

expandmatrix<-function(A,B){
  n1<-dim(A)[1]
  p1<-dim(A)[2]
  n2<-dim(B)[1]
  p2<-dim(B)[2]
  A<-cbind(A,matrix(0,n1,p2))
  B<-cbind(matrix(0,n2,p1),B)
  return(rbind(A,B))
}

## section 4.2 Eq.25  
quantile.weight <- function(r,tau)   # r is the residual and tau is the quantile 
{
  K <- length(tau);
  if(K==1){weight <- 1;}
  else {
    btau <- quantile(r, prob=tau);
    ftau <- rep(0,K);
    ## empirical density estimate based on the empirical quantiles 
    for(jj in 1:K){ ftau[jj] <- density(r,from=btau[jj],to=btau[jj],n=1)$y;}   ## vector a in the paper 

    tmp2 <- matrix(tau,nrow=K, ncol=K);  ## repeat the quantile for K times; by column
    tmp1 <- t(tmp2);  
    tmp <- pmin(tmp1,tmp2)-tmp1*tmp2;
    tmpf <- matrix(ftau,ncol=1);
    m=tmp / (tmpf %*% t(tmpf));   

    d=dim(m)
    
    T=solve(m,diag(d[1]));
    #w1=T %*%ftau        ## optimal weight
    w1 <- solve.QP(Dmat=m,dvec=rep(0,K),Amat=cbind(matrix(1,ncol=1,nrow=K), diag(K)),bvec=c(1,rep(0,K)),meq=1)$solution;
    
    w2=w1/ftau
    weight <- round(w2/sum(w2), 8)
  }
  weight
}

quantile.weight.weighted.QP <- function(r, tau){
  require(quadprog)
  K <- length(tau);
  if(K==1){weight <- 1;}
  else {
    btau <- quantile(r, prob=tau) 
    ftau <- rep(0,K);
    ## empirical density estimate based on the empirical quantiles 
    for(jj in 1:K){ ftau[jj] <- density(r,from=btau[jj],to=btau[jj],n=1)$y;}   
    ## vector a in the paper 
    
    tmp2 <- matrix(tau,nrow=K, ncol=K);  ## repeat the quantile for K times; by column
    tmp1 <- t(tmp2);  
    tmp <- pmin(tmp1,tmp2)-tmp1*tmp2;
    tmpf <- matrix(ftau,ncol=1);
    m=tmp / (tmpf %*% t(tmpf));   
    
    d=dim(m)[1]   
    
    T=solve(m)  ## transpose of M in Bradic et al.
    w <- t(ftau) %*% T %*% (ftau)
    V <- diag(ftau)
    V.t <- solve(V)
    ## quadratic programming using the dual method 
    ## the target function in Koenker Cor. 5.1
    w1 <- solve.QP(Dmat=V.t %*% m %*% V.t, dvec=rep(0,K),
                   Amat=cbind(matrix(1,ncol=1,nrow=K),diag(K)), 
                   bvec=c(1,rep(0,K)), meq=1)$solution
    w2=w1/ftau
    weight <- round(w2/sum(w2), 8)
    
  }
  return(weight)
  
}
# function: scad
# first order derivative of SCAD
scad <- function(lambda,beta,a=3.7)
{
 # theta <- abs(beta);
  #lambda*(theta<=lambda)+(a*lambda-theta)*(a*lambda>theta)/(a-1)*(theta>lambda);
  #return(lambda)
  rep(lambda, length(beta))
}


