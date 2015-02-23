# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2

# NOTES: This script produces the training & test data, trains the 
#        classifier, tests, and plots the final results.

error <- function(y, new_y){
    # returns the fraction of misclassified points
    # for the correct classifications "y" and the new "new_y"
    n <- length(y) * 1.0
    val <- length(new_y[new_y != y]) / n
    return(val)
}


# num dimensions and data points
d = 2
n = 100

# load in the functions
source("classify.R")
source("fakedata.R")
source("perceptrain.R")

# create the data set
z <- c(runif(d + 1, 0.01, 1.0))
z <- norm_z(z)
training_data <- fakedata(z, n)

# now train using perceptron
results <- perceptrain(training_data$S, training_data$y)
prev_c = results$z[3]
results$z = norm_z(results$z)

# make sure perceptrain worked:
new_y <- classify(training_data$S,results$z)
print('fraction of points misclassified with training data')
print(error(training_data$y,new_y))

# now, classify a new batch of data:
test_data <- fakedata(z,n)
new_y <- classify(test_data$S,results$z)
print('fraction of points misclassified with test data')
print(error(test_data$y,new_y))

# now plot everything
# y = -1 is blue
# y = +1 is red
S <- test_data$S
y <- test_data$y
minval <- min( min(c(S[,1],S[,2])))
maxval <- max( max(c(S[,1],S[,2])))
xlim <- c(minval,maxval)
ylim <- c(minval, maxval)
png("classify.png",res=300,width=3,height=3,units='in')

plot(S[,1][y==1],S[,2][y==1], xlab='S_1', ylab='S_2', col='red',bg='red',pch=0,xlim=xlim,ylim=ylim)
points(S[,1][y==-1],S[,2][y==-1],col='blue',bg='blue',pch=2)

z_plot <- z
theta <- atan2(z_plot[2],z_plot[1])
line_length <- 10.0
p1 <- -prev_c * c(-sin(theta), cos(theta))
p2 <- -prev_c * c(sin(theta), -cos(theta))
#p1 <- -results$z[3] * c(-sin(theta), cos(theta))
#p2 <- -results$z[3] * c(sin(theta), -cos(theta))
#f  <- line_length / sum((p2-p1)**2)**0.5
#p1 <- p1*f; p2<- p2*f
#lines(c(p1[1],p2[1]),c(p1[2],p2[2]))

p1 <- -c(z[3]/z[1],0)
p2 <- -c(0, z[3]/z[2])
vdir <- (p2-p1)/(sum((p2-p1)**2)**0.5)
f  <- line_length / sum((p2-p1)**2)**0.5 
#p1 <- p1*f; p2<- p2
p2 <- p2+vdir*f
p1 <- p1-vdir*f
lines(c(p1[1],p2[1]),c(p1[2],p2[2]))

dev.off()
