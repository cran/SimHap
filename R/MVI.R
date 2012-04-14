'MVI' <- function(Betas, Variances,L){
    m <- length(Variances)
    num <- denom <- list()
    Wbar <- 0
    for(j in 1:m){
        num[[j]]<-L%*%Betas[[j]]
        denom[[j]]<-L%*%Variances[[j]]%*%t(L)
        Wbar <- Wbar+denom[[j]]

    }
    num <- t(do.call("cbind",num))
    p <- dim(num)[2]
    Qbar <- apply(num,2,mean)
    Wbar <- Wbar/m
    B <- 0
    for(j in 1:m){
        delta<- num[j,]-Qbar
        B <- B + delta%*%t(delta)
    }
    B <- B/(m-1)

#2 different methods of getting an F-test as in PROC MI documentation
    T0 <- Wbar + B*(1+1/m)
    Fstat <- c(t(Qbar)%*%solve(T0)%*%Qbar/p)
    r <- (1+1/m)*sum(diag(B%*%solve(Wbar)))/p
    nu <- (m-1)*(1+1/r)^2
    p1 <- 1-pf(Fstat, p, nu)

    T2 <- (1+r)*Wbar
    Fstat2 <- c(t(Qbar)%*%solve(T2)%*%Qbar/p)

    tt <- p*(m-1)
    nu1 <- 0.5*(p+1)*(m-1)*(1+1/r)^2
    nu2 <- 4+(tt-4)*(1+1/r*(1-2/tt))^2
    if(tt <=4)
        nu <- nu1
    else
        nu <- nu2
    p2 <- 1-pf(Fstat2, p, nu)
    c(FSmart=Fstat2, Naive=p1,Smart=p2,df1=p, df2=nu)
}
