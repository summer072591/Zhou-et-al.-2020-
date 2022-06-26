## mixture normal 
mixnormal <- function(n, mean1, mean2, sd1, sd2, mix_prob){
  sample_prob <- runif(n)
  
  samples <- sapply(sample_prob, function(x){
                                 if(x <= mix_prob) return(rnorm(1, mean1, sd1))
                                 else return(rnorm(1, mean2, sd1))})
  return(samples)
}
