## Helper functions

# remove cosmic rays for each spectrum
# spc: a intensities vector
removeCosmic <- function(spc) {
    medianFilt <- runmed(spc, 3)
    diff <- spc - medianFilt
    cutoff <- max(2000, 8*sd(diff))
    df <- data.frame(raw=spc, med=medianFilt, diff=diff)
    cosmic <- FALSE
    if (length(which(df$diff>cutoff)) > 0) {
      cosmic <- TRUE
      df[which(df$diff>cutoff),"raw"] <- df[which(df$diff>cutoff),"med"]
    }
    return(list("spc"=df$raw,"cosmic"=cosmic))
}


# search peaks for each spectrum
# size: need to be odd number
findPeaks <- function(spc, size=9, level=0.1) {
    half <- (size-1)/2
    maxI <- max(spc)
    peaks <- c()
    for (i in (half+1):(length(spc)-half)) {
      if ( (spc[i]>=level*maxI) && (spc[i] >= max(spc[(i-half):(i+half)])) ) {
        peaks <- c(peaks, i)
      }
    }
    peaks
}

