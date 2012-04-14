`haplo.label` <-
function(x,numSNPs=2) {

  lab <- ""
  for(i in (numSNPs-1):0) {
    digit<-floor(x/2^i)
    lab <- paste(lab,as.character(digit),sep="")
    x <- x-digit*2^i
  }
  return(lab)
}

