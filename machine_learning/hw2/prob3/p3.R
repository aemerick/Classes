# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2
#
# Problem 3: The linear and non-linear SVM classification
#            of USPS image data 
#

library(e1071) # load in svm library


# ----------------------------------------------
class_error <- function(y, new_y){
    # helper function to compute fraction of 
    # icorrectly classified data
    
    error <- length(y[new_y!=y]) / length(y)
    return(error)
}
# -----------------------------------------------



linear_svm <- function(x, y, cost, kblocks){
    # This uses a k cross validation scheme to compute the best
    # margin parameters to use in a linear SVM model.
    # Inputs:
    # x : data to use as training set
    # y : classifications of training data
    # cost : margin parameters to tune over
    # kblocks: number of cross validation blocks to use

    # Outputs:
    # model : svm best fit model trained on full training data set
    #         with best parameters
    # best_parameters : best margin parameters
    # errors : misclassification errors for each tested parameter
    
    # Find dimensions of the problem and initialize some things
    num_c <- length(cost)
    n     <- length(y) 
    error <- c(rep(0,length(cost)))
    block_size <- n / kblocks # size of the cross validation data blocks

    selection <- sample(1:n) # randomly order the data before cross validation
    for (i in 1:num_c){ # loop over cost parameters
    
        cval <- cost[i]
        sum_error <- 0 # initialize running sum error rate 
        for (k in 1:kblocks){
            # select and set aside the validation data 
            low  <- (k-1)*block_size + 1
            high <-  k   *block_size 
            val_index = selection[low:high]
            x_val <- x[val_index,]
            y_val <- y[val_index]
            
            # everything else is training data
            x_train <- x[- val_index,]
            y_train <- y[- val_index]

            # find the svm model
            model <- svm(type='C-classification', kernel = 'linear',
                         x=x_train, y=y_train, cost = cval)
                         
            # get the new classifcation and error rates of validation             
            new_y <- predict(model, x_val)
            sum_error <- sum_error + class_error(y_val,new_y)
        }# end the cross validation loop
        print(i) 
        sum_error <- sum_error / kblocks # normalize by number of blocks
        error[i]  <- sum_error         
    }# end cost loop

    # find best fit parameters and train data on them
    # if multiple bests, just choose one arbitrarily
    best_c <- cost[ min(error) == error ][1]
    model  <- svm(type='C-classification', kernel = 'linear',
                                            x=x, y=y, cost = best_c)

    return(list(model,best_c, error))
}

non_linear_svm <- function(x, y, costs, gammas, kblocks){
    # This uses a k cross validation scheme to compute the best
    # margin parameters to use in a non-linear SVM model. With 
    # RBF kernel. This is done by a grid search over every
    # possible combination of the margin and gamma parameters.
    #
    # Inputs:
    # x : data to use as training set
    # y : classifications of training data
    # costs : margin parameters to tune over
    # gammas : gamma parameters to tune over for RBF kernel.
    #          gamma = 1.0/(2.0*bandwith^2)
    # kblocks: number of cross validation blocks to use

    # Outputs:
    # model : svm best fit model trained on full training data set
    #         with best parameters
    # best_parameters : best margin parameters
    # errors : misclassification errors for each tested parameter
    # model_params: The entire list of margin and gamma pairs tested


    # initialize error array
    n          <- length(y) 
    num_costs  <- length(costs)
    num_gammas <- length(gammas)

    error        <- c(rep(0, num_costs*num_gammas))
    model_params <- matrix(rep(0,num_costs*num_gammas*2),nrow=num_costs*num_gammas,ncol=2)
    block_size   <- n / kblocks # size of the cross validation data blocks
    selection    <- sample(1:n)

    l<-1 # counter
    for (i in 1:num_costs){  # loop over margin params
        cval <- costs[i]
  
        for (j in 1:num_gammas){ # loop over gammas
            gamma <- gammas[j]      
        
            model_params[l,] <- c(cval,gamma)
            sum_error <- 0 # initialize running sum error rate 
            for (k in 1:kblocks){
                # select and set aside the validation data 
                low  <- (k-1)*block_size + 1
                high <-  k   *block_size 
                val_index = selection[low:high]
                x_val <- x[val_index,]
                y_val <- y[val_index]
            
                # everything else is training data
                x_train <- x[- val_index,]
                y_train <- y[- val_index]

                model <- svm(type='C-classification', kernel = 'radial',
                             x=x_train, y=y_train, cost = cval, gamma=gamma)
                         
                # get the new classifcation and error rates of validation             
#                new_y <- sign(predict(model, x_val))
                new_y <- predict(model,x_val)
                sum_error <- sum_error + class_error(y_val,new_y)
            }# end the cross validation loop

            sum_error <- sum_error/kblocks
            error[l]  <- sum_error
            l <- l + 1
        } # end gamma loop
        print(i)
        sum_error <- sum_error / kblocks # normalize by number of blocks
        error[i] <- sum_error
    }# end cost loop

    # choose the best parameters and fit the model
    best_params <- model_params[error==min(error),]
    best_params <- matrix(best_params,ncol=2, nrow=length(best_params)/2)
    model  <- svm(type='C-classification', kernel = 'radial',
                         x=x, y=y, cost = best_params[1,1], gamma=best_params[1,2])

    return(list(model,best_params,error,model_params))
}

plot_linear_svm <- function(costs, errors, outfile){
    # Plots the results of the linear SVM
    #
    # Input:
    # costs: Margin parameters used
    # errors: classification errors for each margin parameter
    # outfile: filename to save
    pdf(outfile)
    plot(costs,errors, xlab='Margin Parameter (cost)', ylab ='Misclassification Rate',asp=1,xlim=c(2**(-10),1),ylim=c(0,0.5),cex=2)
    min_err <- min(errors)
    min_costs <- costs[min_err == errors]
    min_err <- min_err * c(rep(1,length(min_costs)))
    points(min_costs, min_err, col='red', bg='red',cex=2,xlim=c(2**(-10),1),ylim=c(0,0.25))
    legend(x='topleft',c('All points', 'Best Param'),col=c('black','red'),pch=c(1,1),cex=1.5)
    dev.off()
} # end plot linear SVM

plot_non_linear_svm <- function(parameters, gammas, errors, best_param, outfile){
    # Plots the results of the non-linear SVM
    # Inputs:
    # parameters : the full list of model parameters tested
    # gammas     : the gammas tested
    # error      : the classification errors for each model parameter pair
    # outfile    : output filename (must be .png!)
    #
    pdf(outfile)
    plot(-1,-1,xlab='Margin Parameter (cost)',
            ylab ='Misclassification Rate',xlim=c(2**(-10),1),ylim=c(0,1),asp=1,cex=1.5)

    colors <- c("black","red","blue","green","black")
    pch    <- c(0,1,2,3,18)
    for (i in 1:length(gammas)){
        
        select <- parameters[,2] == gammas[i]
        points(parameters[select,][,1], errors[select], col=colors[i], bg=colors[i],pch=pch[i],cex=1.5)
    }
    points(best_param[1], min(errors), col='black', bg='black',pch=18,cex=2)

    legend(x='topright',c("log(Gamma) = -4","log(Gamma) = -3", "log(Gamma) = -2", "log(Gamma) = -1","Best Param"),
           cex=1.5,pch=pch,col=colors)
    dev.off()	
} # end plot non-linear SVM


# ---------------------------------------------------------------
#
# The below selects the data, finds the SVM's, tests, and plots 
#
# ---------------------------------------------------------------


# load in the raw data
raw_data  <- as.matrix(read.table('uspsdata.txt'))
raw_class <- c(as.matrix(read.table('uspscl.txt')))
n    <- dim(raw_data)[1]

# draw a random sample of the data for testing
test_index <- sample(1:n, 0.2*n, replace=F)
test_data  <- raw_data[test_index,]
test_y     <- raw_class[test_index]

# remaining data to be used as training
train_data <- raw_data[- test_index,]
train_y    <- raw_class[- test_index]

# set the costs and gamma ranges
costs <-  2**seq(-10,0,0.5)
gammas <- 10**c(-4,-3,-2,-1)

#
#
# --------------- Linear SVM ----------------------------
#
#

# find the best linear SVM using k=5 cross validation
print("Finding the best linear svm") 
results <- linear_svm(train_data, train_y, costs, 5)
model   <- results[1][[1]]
error   <- results[3][[1]]

# find the test error
plot_linear_svm(costs, error, "linear_SVM.pdf")
pred_test <- predict(model, test_data)
final_error <- class_error(test_y, pred_test)
print("Final error on test data set for linear SVM")
print(final_error)

#
#
# ------------------- Non - linear (RBF kernel) SVM -----------
#
#
# find the best non-linear SVM using k=5 corss validation
print("finding the best non-linear svm")
non_lin <-non_linear_svm(train_data, train_y, costs, gammas, 5)
nl_model <- non_lin[1][[1]]
best_param <- non_lin[2][[1]]
errors   <- non_lin[3][[1]]
params   <- non_lin[4][[1]]

# calculate the test error
pred_test <- predict(nl_model, test_data)
final_error <- class_error(test_y, pred_test)
print("Final error on test data set for non-linear SVM")
print(final_error)
plot_non_linear_svm(params, gammas, errors, best_param, "non_linear_svm.pdf")

