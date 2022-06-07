
#' this file is created for plotting Figure 1 in the draft

library(ggplot2)
library(reshape2)


quantile.loss <- function(x){
  
#  rho_cqr <- sapply(z, quantile.loss)
  K <- length(nu)
  
  left.end <- c(-Inf, btau[1 : K])
  right.end <- c(btau[1:K], Inf)
  
  left.in <- as.numeric(x >= left.end)
  right.in <- as.numeric(x < right.end)
  
  In <- which(apply(rbind(left.in, right.in), 2, sum) == 2)

  if(In == 1) return(sum(nu * (1 - tau) * abs(x - btau)))
  else if(In == (K + 1)) return(sum(nu*tau*abs(x - btau)))
  else{
     fx <- sum(c((nu * tau * abs(x - btau))[1:(In-1)], (nu * (1 - tau) * abs(x - btau))[In : K]))
     return(fx)
  }
  
  
}


z <- seq(-2.5, 2.5, by = 0.001)

#tau <- 1:5/6
#nu <- rep(1, 5)/5
# single tau
tau <- c(.25, .5, .75)
K <- length(tau)

nu <- 1/rep(K, K) # equal weights
#nu <- c(.15, .55, .3)
btau <- quantile(z, probs = tau)

rho_cqr <- sapply(z, quantile.loss)

quantileloss <- data.frame(x = z, rho_cqr = rho_cqr)



ggplot(quantileloss, aes(x = x, y = rho_cqr)) + geom_line() + ylab("rho(x) ") + 
 # ggtitle("Single quantile loss function with tau = 0.3") +
 # ggtitle(paste("Composite quantile loss functions with K = ", K, "and equal weights")) + 
 # ggtitle(paste("Composite quantile loss functions with K = ", K, "and unequal weights")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 16)) 


## 
tau <- .3
K <- length(tau)
nu <- 1/rep(K, K) # equal weights
#nu <- c(.15, .55, .3)
btau <- quantile(z, probs = tau)

rho_cqr_single <- sapply(z, quantile.loss)

tau <- c(.25, .5, .75)
K <- length(tau)
nu <- 1/rep(K, K) # equal weights
btau <- quantile(z, probs = tau)

rho_cqr_K3equal <- sapply(z, quantile.loss)

nu <- c(.15, .55, .3)
btau <- quantile(z, probs = tau)
rho_cqr_K3unequal <- sapply(z, quantile.loss)


quantileloss <- data.frame(x = rep(z, 3), rho_cqr = c(rho_cqr_single, rho_cqr_K3equal, rho_cqr_K3unequal),
                           loss = rep(c("Single", "Equal weights", "Unequal weights"), each  = length(z)))

quantileloss$loss <- factor(quantileloss$loss, levels = c("Single", "Equal weights", "Unequal weights"))

library(grid)
ggplot(quantileloss, aes(x = x, y = rho_cqr)) + geom_line() + ylab(expression(rho(x))) + facet_wrap(~ loss) + 
  # ggtitle("Single quantile loss function with tau = 0.3") +
   ggtitle(paste("Single versus composite quantile loss function with different weights")) + 
  # ggtitle(paste("Composite quantile loss functions with K = ", K, "and unequal weights")) + 
  ylim(c(0, 1)) + 
  theme(panel.grid.minor = element_blank(), 
       # legend.position = "bottom", legend.text = element_text(size = 16), 
      #  legend.title = element_text(size = 16), 
      plot.margin=unit(c(0,.1,00,.1),"cm"),
      axis.title = element_text(size = 16),
      aspect.ratio = 1.3,
      #legend.position = "none",
        axis.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 16)) 

