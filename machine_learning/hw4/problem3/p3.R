# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 4
#

norm_vec <- function(v){

    if (dim(v)[1] > 1){
        result <- v / sqrt(rowSums(v^2))
    } else{
        result <- v / sqrt(sum(v^2))
    }

    return(v)
}

MultinomialEM <- function(H, K,tau){
    # multinomial EM with histogram matrix H
    # number of classes  K, and threshold tau
    n <- dim(H)[1]

    # step 1: choose K at random and calculate centroids
    selection <- sample(1:K)
    tk <- norm_vec(H[selection,])
#    H  <- norm_vec(H)
    # initialize phi and a 
    phi   <- matrix(rep(0,n*K),nrow=n,ncol=K)
    a     <- matrix(rep(0,n*K),nrow=n,ncol=K)
    preva <- matrix(rep(1,n*K),nrow=n,ncol=K)

    # initialize ck to random values and normalize so
    # sum of ck = 1
    ck  <- runif(K,0,1)
    ck  <- ck / sum(ck)

    delta <- 10*tau
    iter  <- 1
    keeplooping <- TRUE
    max_iter <- 100
    while (keeplooping){
        # Now do the E-step
 #       sum <- c(rep(0,K))
        for (k in 1:K) {
            phi[,k] <- exp( rowSums(H * log(tk[k,])))
            a[,k]   <- ck[k] * phi[,k]

#            sum <- sum + ck[k]*phi[,k]
        }
        #print(phi)
        # normalize a 
        a <- a / rowSums(a)

        # Now do the M step
        ck   <- colSums(a)/n
  #      print(dim(a))
 #       print(dim(H))
#        print(dim(t(a)))
        bk   <- t(a) %*% H
        sumb <- rowSums(bk)
        tk   <- bk / sumb

        # A is the matrix ak
        delta <- norm( a - preva,"O") #L1 matrix norm
        iter  <- iter + 1
        preva <- a
#        print(tau)
 #       print(delta)
  #      print(a[1,])
        if ((delta < tau) || (iter > max_iter)){
            keeplooping <- FALSE
        }
    }
    print(iter)
    m <- c(rep(0,n))
    for (i in 1:n){
        m[i] <- which.max(a[i,])
    }
    print(a[sample(1:10),])
    return(m)
}


# loads in the Histogram matrix
# 40000 histograms, 16 bins per histogram
H <- matrix(readBin("histograms.bin", "double", 640000), 40000, 16)

# adjust for zero H's
H <- H + 0.01
#H <- H[1:10,]

