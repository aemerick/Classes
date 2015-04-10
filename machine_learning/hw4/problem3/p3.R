# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 4
#
#  This file contains the clustering algorithm for problem 3
#  of homework 4
#
#
norm_vec <- function(v){

    if (dim(v)[1] > 1){
        result <- v / rowSums(v)
    } else{
        result <- v / sum(v)
    }

    return(result)
}

MultinomialEM <- function(H, K,tau){
    # multinomial EM with histogram matrix H
    # number of classes  K, and threshold tau
    n <- dim(H)[1] 

    # step 1: choose K at random and calculate centroids
    selection <- sample(1:K)
    H <- H  + 0.01    
    tk <- H[selection,] 
    tk <- tk / rowSums(tk) + 0.01

    # initialize phi and a 
    phi   <- matrix(rep(0,n*K),nrow=n,ncol=K)
    P     <- matrix(rep(0,n*K),nrow=n,ncol=K)
    prevP <- matrix(rep(0,n*K),nrow=n,ncol=K)
 
    # initialize ck to random values and normalize so
    # sum of ck = 1
    ck  <- runif(K,0,1)
    ck  <- ck / sum(ck)

    iter  <- 1
    keeplooping <- TRUE
    max_iter <- 500
    while (keeplooping){
        prevP <- P 
 
        # E-step
        phi <- exp(H %*% log(t(tk)))
        P   <- ck * phi
        P   <- P / rowSums(P)


        # M-step
        ck   <- colSums(P)/n
        bk   <- t(P) %*% H
        tk   <- bk / rowSums(bk)

        # Check for convergence
        delta <- norm( P - prevP,"O") #L1 matrix norm
        iter  <- iter + 1


        if ((delta < tau) || (iter > max_iter)){
            keeplooping <- FALSE
        }
        
    }

    print(sprintf('Converged in %i iterations',iter))
  
    # convert soft assignments to hard assignments
    m <- c(rep(0,n))
    for (i in 1:n){
        m[i] <- which.max(P[i,])
    }

    return(m)
}



