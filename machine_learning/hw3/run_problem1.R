# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 3
#
#  This script contains the functions used to run the AdaBoost 
#  problem. i.e. this is the run script. AdaBoost test is run
#  with do_adaboost to generate errors as a function of B.



source("problem1.R")

do_adaboost <- function(B){
    #
    # Given B, function to load in data and do the AdaBoost
    # algorithm with errors computed via cross validation
    #
    raw_data  <- as.matrix(read.table('uspsdata.txt'))
    raw_class <- c(as.matrix(read.table('uspscl.txt')))
    n         <- dim(raw_data)[1]
    cv_results <- CV_AdaBoost(raw_data,raw_class,B,K=5)
    train_error <- cv_results$train
    test_error <- cv_results$test
    print("Training and test errors")
    print(train_error)
    print(test_error)
    return(list(train=train_error,test=test_error))
}

load_results <- function(filename){
    # load previously computed adaboost results from file
    #
    #
    results <- as.matrix(read.table(filename))
    return(list(b=results[,1],train=results[,2],test=results[,3]))       
}

plot_results <- function(train, test, filename='error.pdf'){
    #
    # Given the arrays for training and test error, plot the
    # results
    #
    pdf(filename)
    B <- length(train)
    b <- c(1:B)
    plot(b, train, xlab='b', ylab='Error', col='black', type ='l',cex.lab=1.2,lwd=3,ylim=c(0.0,0.2))
    lines(b, test,col='blue',cex=2.5,lwd=3)
    
    # averaged final error
    llim <- floor(0.7*B)
    avg_test <- ave(test[llim:B])[1]
    avg_train <- ave(train[llim:B])[1]
    lines(c(1,B),c(avg_train,avg_train),lwd=2.5,col='black',lty=2)
    lines(c(1,B),c(avg_test,avg_test),lwd=2.5,col='blue',lty=2)
   
    legend(x='topright', c(sprintf("Training - Final = %0.3f",avg_train), sprintf("Test - Final = %0.3f",avg_test)), col=c('black','blue'),pch=c('-','-'),cex=1.5)
    dev.off()
    graphics.off()
}

B <- 1

results <- do_adaboost(B)
#results <- load_results("final_output.dat")
#final <- matrix(c(c(1:B), results$train, results$test),nrow=B,ncol=3)
#write(t(final), "final_output.dat", ncolumns=3)
#plot_results(results$train,results$test)

