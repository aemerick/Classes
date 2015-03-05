source("problem1.R")

do_adaboost <- function(B){
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


plot_results <- function(train, test, filename='error.png'){
    png("error.png")
    b <- c(1:B)
    plot(b, train_error, xlab='b', ylab='Error', col='black', type ='l',cex=2.5,lwd=3,ylim=c(0.0,0.2))
    lines(b, test_error,col='blue',cex=2.5,lwd=3)
    legend(x='topright', c("Training Error", "Test Error"), col=c('black','blue'),pch=c('-','-'),cex=1.5)
    dev.off()
    graphics.off()
}

#results <- do_adaboost(14)
#plot_results(results$train,results$test)

