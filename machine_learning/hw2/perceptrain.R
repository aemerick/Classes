# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2

# make sure classifier is loaded
#if (exists(classify) == false) {
source("classify.R")
#}
norm_z <- function(z){
    # normalize z to the length of v_H
    d = length(z) - 1
    len <- sum(z[1:d]**2)**0.5
    
    z_norm <- z / len
    return(z_norm)
}

perceptrain <- function(S,y){
#
#
#    z[i+1] = z[i] - alpha[i] gradient(C_p(z^k))
# need to output the result (z)
# and Z_history, or z at every iteration

    max_iter  <- 1000    # maximum number of iterations
    tolerance <- 1.0E-12 # minimize cost until smaller than this
    d <- ncol(S) - 1

    z <- c(runif(d + 1,0.0, 1.0)) # random number between 0 and 1
    z <- norm_z(z)

    # allocate memory for Z_history matrix, deptermined by max allowed iterations
    Z_history <- matrix( rep(0, (d+1) * (max_iter)), nrow=max_iter, ncol=d + 1)
    Z_history[1,] <- z

    cost <- 1000 # initialize cost to start the loop to > 0

    II <- c(rep(0, length(y)))
    
    # initialize the first classification with the chosen z
    f_x    <- classify(S,z)
    II     <- II * 0
    II[f_x != y] <- 1  

    keep_looping <- TRUE; k <- 1
    while(keep_looping){
        k <- k + 1
        alpha_k <- 1.0 / k
        
        # gradient in the cost function
        grad_c = colSums(-y * S * II)

        # calculate new z and save 
        z             <- Z_history[k-1,]  - alpha_k * grad_c
        Z_history[k,] <- z

        # evaluate the cost of the new z
        f_x <- classify(S,z)
        II <- II* 0 ; II[f_x != y] <- 1

        cost = sum(II * abs(colSums( z * t(S))))

        # keep looking?
        if ((cost < tolerance) || (k > max_iter)) {
            keep_looping <- FALSE
        }
        
    } # end while loop

    print("number of iterations needed")
    print(k)
 
    # remove the unused Z_history rows 
    Z_history <- Z_history[1:k,]
    return (list(z=z, Z_history=Z_history))
}
