# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2

# make sure classifier is loaded
#if (exists(classify) == false) {
source("classify.R")
#}

perceptrain <- function(S,y){
#
#
#    z[i+1] = z[i] - alpha[i] gradient(C_p(z^k))
# need to output the result (z)
# and Z_history, or z at every iteration

    max_iter  <- 10000    # maximum number of iterations
    tolerance <- 1.0E-12 # minimize cost until smaller than this
    d <- ncol(S) - 1

    k <- 1
    z <- c(runif(d + 1,0.0, 1.0)) # random number between 0 and 1
    prev_z <- c(rep(0, d+1)) # zeros of size n

    Z_history <- matrix( rep(0, (d+1) * (max_iter)), nrow=max_iter, ncol=d + 1)
    Z_history[1,] <- z

    cost <- 1000 # initialize cost to start the loop to > 0

    II <- c(rep(0, length(y)))

    while(cost > tolerance){
        k <- k + 1
        alpha_k <- 1.0 / k
        # do not need to sum elements where classification is correct
        # only those where it is wrong (i.e. where  f(x) =/= y),
        # since correct classifications give 0
        f_x    <- classify(S,z)
        II     <- II * 0
        II[f_x != y] <- 1   

         
   #     grad_c <- c(rep(0,1+d))
  #      for ( i in 1:length(y)){
 #           grad_c <- grad_c   -y[i]*S[i,] * II[i]
#        }
        grad_c = colSums(-y * S * II)
        #print(grad_c)
        z             <- prev_z  - alpha_k * grad_c
        prev_z        <- z
        Z_history[k,] <- z

        # evaluate the current cost
        f_x <- classify(S,z)
        II <- II* 0 ; II[f_x != y] <- 1

 #       for (i in 1:length(y)){
#            cost <- sum( abs(t(z) %*% S[i,]) )        
  #      }
        cost = sum(rowSums(abs(z * S)))
        print(cost)
    } # end while loop

    print(k)
    # return 
    Z_history = Z_history[1:k,]
    return (list(z=z, Z_history=Z_history))
}
