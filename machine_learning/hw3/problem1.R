# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 3
#
# Problem 1: AdaBoost implementation with error estimation
#            via K-fold cross validation.
#
#    This code contains all functions asked for explicitly in
#    the assignment, as well as a few helper functions. Function
#    names should be descriptive, but code is commented.
#
#-------------------------------------------------------------------------

calc_error <- function(y, new_y, w=NULL){
    # If no weighting is given, returns the fraction of misclassified
    # points between two arrays of class labels. If array of weights is
    # supplied, calculates the "weighted error" instead:
    #    error = sum(w_i * II_i)/sum(w_i) where II is 1 if point is 
    #             misclassified, zero if classified correctly
    # Inputs
    # y     : correct class labels 
    # new_y : test class labels
    # w     : Optional weight array. Default is NULL (i.e. no weighting)
    # 
    # Outputs
    # error : Fraction of misclassified points OR weighted error 
    # -------------------------------------------------------------------
 
    n <- length(y)

    if (is.null(w)){ # no weights, mis. fraction
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
    #
    # Training function for a single decision stump, implemented
    # as per the assignment. Performs grid search over all axis,
    # a range of theta values for each, for m=+1 and m=-1, to 
    # choose best (j,theta,m) parameters. n+1 thetas are searched over
    # as determined by the data points in each axis. Test thetas
    # are placed 1/2 between each data point in each axis, + 2 additional
    # thetas at a point less than the min value and a point greater than max.
    #
    # Input
    # X    : matrix containing data points (rows) and their axis (columns)
    # w    : weights used to train the decision stump
    # y    : correct class labels for data
    #
    # Output
    # pars : list with names "j", "theta", and "m" for parameter set that
    #        minimized weighted classification error
    # --------------------------------------------------------------------

    # 
    if (length(unique(y)) > 2){
        print(length(unique(y)))
        print("ERROR > 2 CLASS LABELS - BINARY CLASSIFIER ONLY")
        return()
    }

    # dimensions of X
    n <- length(y)
    d <- dim(X)[2]

    # for each dimension, compute the best theta
    # choose theta's at 1/2 between each data point
    # 
    # first, sort all columns in ascending order
    sorted_X       <- apply(X,2,sort)
    # initialize theta matrix
    possible_theta <- matrix(rep(0,((n + 1)*d)),nrow=n+1,ncol=d)

    # first and last values should be < min and max value in each axis
    # by how much is arbitrarily scaled by spacing between points
    possible_theta[1,  ] <- X[1,] - 0.5*(X[2,]-  X[1,])
    possible_theta[n+1,] <- X[n,] + 0.5*(X[n,]-X[n-1,]) 
    # rest of the values are 1/2 in between every point
    possible_theta[2:n,] <- 0.5*(X[2:n,] + X[1:n-1,])

    # now for each value of theta in each axis, calculate error
    # assuming m = 1 or -1. Initialize error to large value
    all_err_mup   <- matrix(rep(99999,((n+1)*d)),nrow=n+1,ncol=d)
    all_err_mdown <- 1.0 * all_err_mup

    # loop over every j,theta combination and calc weighted error
    for (j in 1:d){
        for (i in 1:n+1){
        all_err_mup[i,j]   <-calc_error(y,
                              classify(X,list(j=j,theta=possible_theta[i,j], m=1)),
                              w=w)
        all_err_mdown[i,j] <- calc_error(y,
                              classify(X,list(j=j,theta=possible_theta[i,j],m=-1))
                              ,w=w)
        }
    }

    # find the argmin of the two
    argmin_mup   <- which.min(all_err_mup  )
    argmin_mdown <- which.min(all_err_mdown)
   
    # decide which choice of m was better
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
    y <- c(rep(pars$m,dim(X)[1]))

    # assign -m to y indeces where  x_j < theta
    y[X[,pars$j] <= pars$theta] <- -pars$m

    return(y)
} 

AdaBoost <- function(X, y, B){
    #
    # AdaBoost algorithm. Trains B decision stumps on data (X) with
    # correct class labels (y).
    #
    # Inputs
    # X     : data matrix with data points in rows and axes in columns
    # y     : correct class labels for data points
    # B     : number of decision stumps to train
    #
    # Outputs
    # results : Single list containing "alphas" voting weights for
    #           aggregated classifier, "allPars" for all decision stump
    #           parameters (rows) with colums of j, theta, and m. 
    #           "train_error" or misclassification rate as a function of b
    # ---------------------------------------------------------------------

    n <- length(y)
    d <- dim(X)[2]

    # init all parameters matrix
    allPars <- matrix(rep(0,3*B), nrow=B, ncol=3)

    # initialize weights and alphas
    w     <- c(rep(1/n,n))
    alpha <- c(rep(0,B))
    II    <- c(rep(0,n))
 
    # initialzie errors
    train_error <- c(rep(0,B))

    print("---- AdaBoost ----")
    for (b in 1:B){
        # train the decision stump
        c_pars  <- train(X,w,y)

        # classify with new DS, compute error and voting weights
        y_star  <- classify(X,c_pars)
        c_error <- calc_error(y, y_star, w=w) 
        alpha[b]   <- log( (1.0 - c_error)/c_error)

        # recompute new weights on data points (w)
        II <- II*0
        II[y_star != y] <- 1
        w  <- w * exp( alpha[b] * II)
    
        # save parameters
        allPars[b,] <- c(c_pars$j,c_pars$theta,c_pars$m)
        
        # running tally of training error as a function of b
        agg_y          <- agg_class(X,alpha[1:b],allPars[1:b,])
        train_error[b] <- calc_error(y, agg_y)

    } # end B

    return(list(alpha=alpha, allPars=allPars, train_error = train_error))
}

CV_AdaBoost <- function(X, y, B, K = 5){
    # Wrapper around the AdaBoost algorithm to perform K-fold 
    # cross validation purely for the purpose of estimating the 
    # error / misclassification rate as a function of b in order
    # to determine convergence properties.
    #
    # Inputs
    # X     : data matrix. n rows of points. d columns of axes
    # y     : correct class label
    # B     : Number of DS to train
    # K     : Optional. Number of K-folds. Default is 5.
    #
    # Outputs
    # List containing "train" and "test" errors as a function of b
    # ------------------------------------------------------------

    n <- dim(X)[1]
    d <- dim(X)[2]

    # initialize arrays for error
    train_error <- c(rep(0,B))
    test_error  <- c(rep(0,B))

    block_size <- n / K
    selection  <- sample(1:n) # random ordering for cross validation
    for (k in 1:K){ # loop over cross validation folds
        # select index range for test data in fold k
        low        <- (k-1)*block_size +1
        high       <- k * block_size
        test_index <- selection[low:high]
        test_x     <- X[test_index,]
        test_y     <- y[test_index]
        # rest of data is training
        train_x    <- X[- test_index,]
        train_y    <- y[- test_index]

        # run AdaBoost on training data...
        ada_results <- AdaBoost(train_x, train_y, B)

        # Classify test data and calculate errors
        # training data errors are already computed in AdaBoost
        for (b in 1:B){
            test_error[b] <- test_error[b] +
                             calc_error(
                               agg_class(test_x, ada_results$alpha[1:b], ada_results$allPars[1:b,]),
                               test_y)
        }
        print(calc_error(agg_class(test_x, ada_results$alpha[1:B], ada_results$allPars[1:B,]),test_y))
        print(ada_results$train_error)
        train_error <- train_error + ada_results$train_error
    }

    # now average the errors together over the K folds
    train_error <- train_error / K
    test_error  <- test_error / K

    return(list(train=train_error,test=test_error))
}


agg_class <- function(X, alpha, allPars){
    # Aggrigated classifier from the AdaBoost algorithm using
    # weighted decision tree classifiers. 
    #    
    # Inputs
    # X       : data matrix to classify
    # alpha   : array of voting weights for each decision stump
    # allPars : matrix containing decision stump parameters, with
    #           columns of (j,theta,m)
    #
    # Outputs
    # ystar   : final classification
    # ------------------------------------------------------------

    B <- length(alpha)
    n <- dim(X)[1]

    # initialize running total
    running_sum <- c(rep(0,n))
 
    # make sure allPars is matrix (if only 1 DS is used) 
    if ((length(allPars) == 3)){ 
        allPars <- t(matrix(allPars))
    }

    # loop over all decision stumps to compute total weighted vote
    for (b in 1:B){
        running_sum <- running_sum + 
                          alpha[b]*classify(X,list(j=allPars[b,1],theta=allPars[b,2],m=allPars[b,3]))
    }

    # classification is the sign of the weighted vote
    return(sign(running_sum))
}


