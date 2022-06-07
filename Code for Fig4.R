#'**********************************************************************************
#' This is the R code for generating Figure 4 the signal recovery plot
#' of the model-averaged estimator using Bates-Granger type weights
#'**********************************************************************************

library(ggplot2)
library(reshape2)
library(wavethresh)
K <- 3


## save the original audio signal 
beta.original <- as.numeric(read.table("wave_sound.txt", sep = " ", header = FALSE)[(6*1024+1):(8*1024)])
beta.wave <- wd(beta.original, filter.number = 8, family = "DaubLeAsymm")$D



#' read in data 
ramp.recov.BG.1 <- as.numeric(read.table("modavrg_beta_rampt3.txt", sep = " ", fill = TRUE, 
                                            header = FALSE))[-1]
ramp.recov.BG.2 <- as.numeric(read.table("modavrg_beta_rampmixnorm.txt", sep = " ", fill = TRUE, 
                                            header = FALSE))[-1]


wave.ramp.recov.BG.1 <- wr(ramp.recov.BG.1, filter.number = 8, family = "DaubLeAsymm")
wave.ramp.recov.BG.2 <- wr(ramp.recov.BG.2, filter.number = 8, family = "DaubLeAsymm")


wave.data <- wr(as.numeric(read.table("wave_sound_78.txt", sep = " ", header = FALSE)), 
                filter.number = 8, family = "DaubLeAsymm")

comb.est <- cbind(wave.ramp.recov.BG.1, wave.ramp.recov.BG.2)

reconst.melt <- data.frame(cbind(melt(comb.est), rep(wave.data, 2)))


colnames(reconst.melt) <- c("Sampling", "Error", "Amplitude", "Original_Signal")

levels(reconst.melt$Estimator) <- c("t3", "mixnormal")


reconst.melt$Estimator <- as.factor(reconst.melt$Error)
plot <- ggplot(data = reconst.melt, aes(x = Sampling, y = Original_Signal)) + facet_wrap(~ Error) + 
  geom_path() + geom_path(aes(x = Sampling,  y = Amplitude, colour = Error, linetype = Error)) + 
  ylab("Amplitude ") + 
  ggtitle(paste("Reconstructed audio signal fractions by regularized MAQR with K = ", K, "using Bates-Granger type weight")) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 16)) 



