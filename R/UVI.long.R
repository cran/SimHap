'UVI.long' <- function(Betas,Variances,v0,ADJ=FALSE){
    m <- dim(Betas)[1]
    ses <- as.matrix(sqrt(Variances))
    coefs <- as.matrix(Betas)
    W <- colMeans(ses^2)
    Bhat <- colMeans(coefs)
    B.diff <- apply(coefs, 1, "-", Bhat)
    B <- apply(B.diff, 1, var)
    pooled <- list(coefficients = NULL, se = NULL,t.stat=NULL,
                   p.value=NULL, df=NULL, ratio=NULL)

    gamma <- (1 + 1/m)*B/(W + (1+1/m)*B)
    vobs <- (1-gamma)*v0*(v0+1)/(v0+3)
    vm <- (m-1)*(1+W/(B*(1+1/m)))^2

    pooled$coefficients <- Bhat
    pooled$se <- sqrt(W + (1 + 1/m) * B)
    pooled$ratio <- (1+1/m)*B/W
    if(ADJ){
      pooled$df <- (1/vm + 1/vobs)^(-1)
      } else{
        pooled$df <- vm
        if(any(vm > v0)){
          warning("The computed degrees of freedom is larger than the complete-data degrees of freedom. Consider using the adjusted degrees of freedom.")
        }
      }
    
    pooled$t.stat <- pooled$coefficients/pooled$se
    pooled$p.value <- 2*pt(-abs(pooled$t.stat),df=pooled$df)
    pooled
}
