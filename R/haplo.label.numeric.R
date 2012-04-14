`haplo.label.numeric` <-
function(x,numSNPs=2) {
  len <- length(x)
  lab <- matrix(0,nrow=len,ncol=numSNPs)

    for(i in (numSNPs-1):0){
      digit<-floor(x/2^i)
      lab[,numSNPs-i]<-digit
      x <- x-digit*2^i
      }
    return(lab)
}

