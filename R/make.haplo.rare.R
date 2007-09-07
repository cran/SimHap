`make.haplo.rare` <-
function(infer.object, min.freq){

  if(is.null(infer.object$initFreq)) initFreq <- infer.object$initFreq.controls
  else initFreq <- infer.object$initFreq
  
  if(is.null(infer.object$hap.freq)) hap.freq <- infer.object$hap.freq.controls
  else hap.freq <- infer.object$hap.freq
  freqs2 <- hap.freq[,1:2]
  rare.ind <- initFreq<min.freq #flag rare haplos 
  rare.haplos <- "no rare haplos"
  haplos.names <- names(initFreq)

  if(sum(rare.ind)>=1) { #then grouping of rare to be done *in design matrix only*
    rare.haplos<-haplos.names[rare.ind]

  }

  hapMat <- infer.object$hapMat
  hapMatRare <- hapMat

        # change rare haps to "rare"
        for(i in 1:nrow(hapMat)) {
          for (j in 1:length(rare.haplos)) {
            if (hapMat[i,2] == rare.haplos[j])
              hapMatRare[i,2] <- c("rare")
            if (hapMat[i,3] == rare.haplos[j])
              hapMatRare[i,3] <- c("rare")
          }
        }
   if(any(initFreq<min.freq)) {
   	freqs.rare <- rbind(freqs2[freqs2$Frequency>=min.freq,], c("rare", 1-sum(freqs2$Frequency[freqs2$Frequency>=min.freq])))
	}
	else freqs.rare <- freqs2[freqs2$Frequency>=min.freq,]
   
   #freqs.rare <- rbind(freqs2[freqs2$Frequency>=min.freq,], c("rare", 1-sum(freqs2$Frequency[freqs2$Frequency>=min.freq])))
   freqs.rare$Frequency <- as.numeric(as.vector(freqs.rare$Frequency))
   infer2 <- list(hapMat=infer.object$hapMat, hap.freq=infer.object$hap.freq, initFreq=infer.object$initFreq, final.freq=freqs.rare)
   haplist <- list(hapData=hapMatRare, hapObject=infer2)
   return(haplist)
}

