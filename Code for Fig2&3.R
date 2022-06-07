#'**********************************************************************************
#' This is an example R code for generating Figure 1 & 2
#'**********************************************************************************

library(ggplot2)
library(reshape2)

K <- 3
comb.est <- cbind(wave.ramp.recov, wave.opt.recov, wave.eq.recov, wave.lasso.recov)
reconst.melt <- data.frame(cbind(melt(comb.est), rep(wave.data, 4)))

colnames(reconst.melt) <- c("Sampling", "Estimator", "Amplitude", "Original_Signal")

levels(reconst.melt$Estimator) <- c("MAQR(w_1)", "MAQR(w_2)", "MAQR(w_eq)", "Lasso")


reconst.melt$Estimator <- as.factor(reconst.melt$Estimator)
plot <- ggplot(data = reconst.melt, aes(x = Sampling, y = Original_Signal)) + facet_wrap(~ Estimator) + 
  geom_path() + geom_path(aes(x = Sampling,  y = Amplitude, colour = Estimator, linetype = Estimator)) + 
  ylab("Amplitude ") + 
  ggtitle(paste("Reconstructed audio signal fractions using regularized MAQR with K = ", K, "and the Lasso")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 16)) 
