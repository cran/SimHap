`codeDM` <-
function(haplos,haplo.labels){

  n <- length(haplos)
  numSNPs <- ncol(haplo.labels)

  ans1 <- t(haplos[1:(n/2)]==t(haplo.labels))
  ans2 <- t(haplos[(n/2+1):n]==t(haplo.labels))
  ans11 <- ans1[,1]
  ans22 <- ans2[,1]
  for(i in 2:numSNPs) {
     ans11 <- ans11&ans1[,i]
     ans22 <- ans22&ans2[,i]
  }
  ans=ans11+ans22

  return(ans)
}

