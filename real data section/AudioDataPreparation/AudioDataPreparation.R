
#' ******************************************************************************************************
#'  This document is used to generate the wavelet coefficient used in the real data example section
#'  The generated coefficient is save under the name "wave_sound_78.txt"
#'  
#'  the original audio piece was stored under the name "wave_sound.txt"
#' ******************************************************************************************************


library(signal)
library(wavethresh)

## save the original audio signal 
data(wav)
write.table(t(wav$sound), "wave_sound.txt", sep = " ", row.names = FALSE, col.names = FALSE)



## obtain the wavelet coeffient
beta.original <- as.numeric(read.table("wave_sound.txt", sep = " ", header = FALSE)[(6*1024+1):(8*1024)])

write.table(t(wd(beta.original, filter.number = 8, family = "DaubLeAsymm")$D), 
            "wave_sound_78.txt", sep = " ", row.names = FALSE, col.names = FALSE)

## read in the wavelet coefficient
beta.wave <- wd(beta.original, filter.number = 8, family = "DaubLeAsymm")$D