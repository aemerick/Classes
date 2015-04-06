# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 4
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
    print(dim(H))

    # step 1: choose K at random and calculate centroids
    selection <- sample(1:K)
    selection <- c(50,6000,10000)
    H <- H + 0.01
    H  <- H / rowSums(H)
    
    tk <- H[selection,] 
    tk <- tk / rowSums(tk)
    
    print('tk dim');    print(dim(tk))

    # initialize phi and a 
    phi   <- matrix(rep(0,n*K),nrow=n,ncol=K)
    P     <- matrix(rep(0,n*K),nrow=n,ncol=K)
    prevP <- matrix(rep(0,n*K),nrow=n,ncol=K)
 
    # initialize ck to random values and normalize so
    # sum of ck = 1
    ck  <- runif(K,0,1)
    ck  <- c(1,0.5,0.33333333)
    ck  <- ck / sum(ck)

    iter  <- 1
    keeplooping <- TRUE
    max_iter <- 500
    while (keeplooping){
        prevP <- P
        # Now do the E-step
        for (kk in 1:K) {
            phi[,kk] <- exp( H %*% (log(tk[kk,])))
           
            P[,kk]   <- ck[kk] * phi[,kk]
        }

        # normalize a 
        P <- P / rowSums(P)

        # Now do the M step
        ck   <- colSums(P)/n

        bk   <- t(P) %*% H

        tk   <- bk / rowSums(bk)

        # A is the matrix ak
        delta <- norm( P - prevP,"O") #L1 matrix norm
        iter  <- iter + 1


        if ((delta < tau) || (iter > max_iter)){
            keeplooping <- FALSE
        }
        print(tk)
    }

    print(iter)
    m <- c(rep(0,n))
    for (i in 1:n){
        m[i] <- which.max(P[i,])
    }
    print(P[sample(1:30),])
    return(m)
}


# loads in the Histogram matrix
# 40000 histograms, 16 bins per histogram
H <- matrix(readBin("histograms.bin", "double", 640000), 40000, 16)

# adjust for zero H's
#H <- H + 0.001
#H <- H[1:10,]

