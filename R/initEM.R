`initEM` <-
function(haplo.object, freq=FALSE, maxit=100, tol=1/(2*sum(haplo.object$wt)*100)){

 hapMat <- haplo.object$hapMat

 ID <- haplo.object$ID

 N <- sum(haplo.object$wt)
 wts <- haplo.object$wt
 
 # Initial freq values, if no freq specified, calculate freq values based
 # on augmented dataset.

 if (freq==FALSE){
    allHaps<-c(hapMat[,1],hapMat[,2])
    allWts<-c(wts,wts)
    freq<-tapply(allWts,allHaps,sum)/(2*N)
 }

 freqdiff <- 1
 it <- 1
 num.prob <- vector(length=nrow(hapMat))

 # The EM loop

 while ( (it<maxit) && (freqdiff>tol) ){
   
 # multiplicative constant = 2 for heterozyg 1 for homozyg            
  haplo.probs <- rep(1,nrow(hapMat))+(as.numeric(haplo.object$hapMat[,1] != haplo.object$hapMat[,2]))
  haplo.probs <- haplo.probs*freq[hapMat[,1]]*freq[hapMat[,2]]
  num.prob <- haplo.probs

   # E step: Calculate the weights for everyone
   # Use the ID to determine the number of pseudo-individuals in the
   # denominator probability

        for (i in 1:nrow(hapMat)){               
            pseudo.index <- ID==ID[i]
            wts[i] <- num.prob[i]/sum(num.prob[pseudo.index])
        }

  # M step: Find new estimates using weighted haplotype counts
        allWts <- c(wts,wts)
        freqNew <- tapply(allWts,allHaps,sum)/(2*N)
        freqdiff <- max(abs(freq-freqNew), na.rm=TRUE)#maximum diff
        freq <- freqNew

        it <- it+1
 }

 if(freqdiff>tol) 
   warning(paste("no convergence in EM after ",maxit,"iterations\n"))
   
 freq <- freq[freq!= 0]
 results <- list(freq=freq, wts=wts)
 return(results)
}

