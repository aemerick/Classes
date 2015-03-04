# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 3
#
# 

# 

calc_error <- function(y, new_y, w=NULL){
    # returns the fraction of misclassified points
    # for the correct classifications "y" and the new "new_y"
    #
    # if a weighting is supplied, then the weighted "error" is returned 
    # instead. 

    n <- length(y)

    if (is.null(w)){
        val <- length(new_y[new_y != y]) / n
    } else {
        # do the weighted error
        if (length(w) == 1){
            w <- c(rep(w,n))
        }  

        # initialize 0 vector... fill with 1's for misclassifications
        II <- c(rep(0,n))
        II[new_y !=y] <- 1.0
   
        # weighted error
        val <- sum(w*II)/sum(w)
    }

    return(val)
}


train <- function(X, w, y){

    # get the classifier label ( y is +m or -m) 
    if (length(unique(y)) > 2){
        print(length(unique(y)))
        print("ERROR > 2 CLASS LABELS - BINARY CLASSIFIER ONLY")
        return()
    }

    n <- length(y)
    d <- dim(X)[2]

    # for each dimension, compute the best theta
    # choose theta's as a random number between each data point
    # ... or make it deterministic and pick halfway...
    possible_theta <- matrix(rep(0,((n + 1)*d)),nrow=n+1,ncol=d)
    sorted_X       <- apply(X,2,sort)
    
    #print(dim(possible_theta))
    #print(dim(X))


    # first and last values should be < min and max value in each axis
    # by how much is arbitrarily scaled by spacing between points
    possible_theta[1,  ] <- X[1,] - 0.5*(X[2,]-  X[1,])
    possible_theta[n+1,] <- X[n,] + 0.5*(X[n,]-X[n-1,]) 
    # rest of the values are 1/2 in between every point
    possible_theta[2:n,] <- 0.5*(X[2:n,] + X[1:n-1,])

    #print(possible_theta)

    #print(dim(possible_theta))
    #print(dim(X))
    # now for each value of theta in each axis, calculate error
    # assuming m = 1 or -1
    all_err_mup   <- matrix(rep(99999,((n+1)*d)),nrow=n+1,ncol=d)
    all_err_mdown <- 1.0 * all_err_mup

    #

    for (j in 1:d){
        for (i in 1:n+1){
        all_err_mup[i,j]   <- calc_error(y,classify(X,list(j=j,theta=possible_theta[i,j], m=1)),w=w)
        all_err_mdown[i,j] <- calc_error(y,classify(X,list(j=j,theta=possible_theta[i,j],m=-1)),w=w)
        
        }
    }

   # print(i)
    #print(j)
    # now find the argmin of each 
   # print("minumum errors for +1 m and -1 m")
    #print(max(all_err_mup))
    #print(min(all_err_mdown))
    argmin_mup   <- which.min(all_err_mup  )
    argmin_mdown <- which.min(all_err_mdown)
   
    

    min_m <- which.min(c(all_err_mup[argmin_mup],all_err_mdown[argmin_mdown]))

    # if min_m is 1, then m = +1
    # depending on sign of m, use above argmins to get best j and theta
    if (min_m == 1){
        m_star     <- 1
        j_star     <- col(all_err_mup)[argmin_mup]
        theta_star <- possible_theta[argmin_mup]
    } else { # then m is -1
        m_star     <- -1
        j_star     <- col(all_err_mdown)[argmin_mdown]
        theta_star <- possible_theta[argmin_mdown]
    }
    

    return(list(j=j_star,theta=theta_star,m=m_star))
}

classify <- function(X, pars){
    # Classify function is an implementation of the decision stump
    # classifier given in equation (1) of the homework.
    #
    # Inputs:
    # X : data to classify, with dimensions n x d, n is number of data
    #     points and d is number of axes in each data point
    # pars : list of classifier parameters with names j, theta, and m
    #        for the decision stump axis (j) split value (theta),
    #        and binary classifier value m
    #
    # Outputs:
    # y : array of classifier labels with dimensions n
    # --------------------------------------------------

    # initialize all of y to positive m
    y <- c(rep(pars$m,dim(test_data)[1]))

    # assign -m to y indeces where  x_j < theta
    y[X[,pars$j] <= pars$theta] <- -pars$m

    return(y)
} 

AdaBoost <- function(X, y, B){
    # boosting algortithm

    n <- length(y)
    d <- dim(X)[2]

    allPars <- matrix(rep(0,3*B), nrow=B, ncol=3)

    # initialize weights and alphas
    w     <- c(rep(1/n,n))
    alpha <- c(rep(0,B))
    II    <- c(rep(0,n))
    

    for (b in 1:B){
        c_pars  <- train(X,w,y)

        y_star  <- classify(X,c_pars)
        c_error <- calc_error(y, y_star, w=w)

        alpha[b]   <- log( (1.0 - c_error)/c_error)

        # new weights
        II <- II*0
        II[y_star != y] <- 1
        w  <- w * exp( alpha[b] * II)

        allPars[b,] <- c(c_pars$j,c_pars$theta,c_pars$m)
    }


    return(list(alpha=alpha, allPars=allPars))
}


agg_class <- function(X, alpha, allPars){
    # Aggrigated classifier from the AdaBoost algorithm using
    # weighted decision tree classifiers. 
    #    

    B <- length(alpha)
    n <- dim(X)[1]

    # initialize running total
    running_sum <- c(rep(0,n))
 
    # loop over all decision stumps to compute total weighted vote 
    print("Entering AdaBoost Loop --- printing progress")
    for (b in 1:B){
        print(b)
        #print(allPars[b,])
        running_sum <- running_sum + alpha[b]*classify(X,list(j=allPars[b,1],theta=allPars[b,2],m=allPars[b,3]))
    }

    # classification is the sign of the weighted vote
    return(sign(running_sum))
}


# load in the raw data
raw_data  <- as.matrix(read.table('uspsdata.txt'))
raw_class <- c(as.matrix(read.table('uspscl.txt')))
n         <- dim(raw_data)[1]

# draw a random sample of the data for testing
test_index <- sample(1:n, 0.2*n, replace=F)
test_data  <- raw_data[test_index,]
test_y     <- raw_class[test_index]

# remaining data to be used as training
train_data <- raw_data[- test_index,]
train_y    <- raw_class[- test_index]


# some tests
B = 11
print("adaboostin")
ada_results <- AdaBoost(train_data, train_y, B)
print(ada_results)
final_class <- agg_class(test_data, ada_results$alpha, ada_results$allPars)

print("error using n decision stumps")
print(B)
print(calc_error(test_y,final_class))
