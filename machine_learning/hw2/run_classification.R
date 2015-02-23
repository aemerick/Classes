# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2

# NOTES: This script generates a training set and test set
#        of data from fakedata, uses perceptrain to develop
#        a classifier using the training data, and tests
#        it on the test data. Finally, plot the results.


# load in the functions
source("classify.R")
source("fakedata.R")
source("perceptrain.R")

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

# choose a z randomly, normalize, generate training data
z <- c(runif(d + 1, 0.01, 1.0))
z <- norm_z(z)
ztrue <- z
training_data <- fakedata(z, n)

# now train using perceptron... normalize final z
results <- perceptrain(training_data$S, training_data$y)
prev_c    <- results$z[3]
results$z <- norm_z(results$z)

# make sure perceptrain worked:
new_y <- classify(training_data$S,results$z)
print('fraction of points misclassified with training data')
print(error(training_data$y,new_y))

# now, generate and classify a new batch of data:
test_data <- fakedata(z,n)
new_y <- classify(test_data$S,results$z)
print('fraction of points misclassified with test data')
print(error(test_data$y,new_y))

################################################################
#
# THE BELOW IS MOSTLY GROSS CODE TO PLOT THE RESULTS....
#
################################################################

# now plot everything
# y = -1 is blue
# y = +1 is red
S <- training_data$S
y <- training_data$y

# find the x1 and x2 limits to plot
minval <- min( min(c(S[,1],S[,2])))-1 
maxval <- max( max(c(S[,1],S[,2])))+1
xlim <- c(minval,maxval)
ylim <- c(minval, maxval)

pdf("classify.pdf") # open a file
plot(S[,1][y==1],S[,2][y==1], xlab=expression('S'[1]), ylab=expression('S'[2]), col='red',bg='red',pch=0,xlim=xlim,ylim=ylim,asp=1)
points(S[,1][y==-1],S[,2][y==-1],col='blue',bg='blue',pch=2)


line_length <- 50.0 # how long to draw the hyperplanes
zhist <- results$Z_history
nz    <- length(zhist)/3
select <- ceiling(seq(1,nz,length=6)) # select some z_hists to plot
colors <- c('red', 'orange', 'green', 'blue', 'purple','black','blue','red')
legend_string <- c(rep('',8)) # init legend format string
pch <- c(rep('-',6))          # points used

# this loops over chosen z histories and plots the hyperplanes
for (i in 1:length(select)){
    l  <- select[i]
    zk <- zhist[l,]
    p1 <- -c(zk[3]/zk[1],0)    # start point of line
    p2 <- -c(0, zk[3]/zk[2])   # end point of line
    vdir <- (p2-p1)/(sum((p2-p1)**2)**0.5)  # the line vector
    f  <- line_length / sum((p2-p1)**2)**0.5 # scaling the line
    p2 <- p2+vdir*f            # final end point
    p1 <- p1-vdir*f            # final start point
    lines(c(p1[1],p2[1]),c(p1[2],p2[2]),col=colors[i], lwd=2)
    legend_string[i] <- sprintf("z[%i]",l)
    
}

legend_string[6] <- sprintf("z[%i] - final",l)
legend_string[7] <- "Training Data: -1"
legend_string[8] <- "Training Data: +1"

# legend for scatter plot points and the lines
legend(x="bottomleft",legend_string[1:6],cex=1.5,col=colors[1:6],pch=pch[1:6])
legend(x="topright",legend_string[7:8],cex=1.5,col=colors[7:8],pch=c(2,0))

# done!
dev.off()
graphics.off()
# ------------------------------------------------------------------------
# make a new plot with test data 
# now plot everything
# y = -1 is blue
# y = +1 is red
S <- test_data$S
y <- test_data$y
zf <- results$z

# find the x1 and x2 limits to plot
minval <- min( min(c(S[,1],S[,2])))-1
maxval <- max( max(c(S[,1],S[,2])))+1
xlim <- c(minval,maxval)
ylim <- c(minval, maxval)

pdf("classify_final.pdf") # open a file
plot(S[,1][y==1],S[,2][y==1], xlab=expression('S'[1]), ylab=expression('S'[2]), col='red',bg='red',pch=0,xlim=xlim,ylim=ylim,asp=1)
points(S[,1][y==-1],S[,2][y==-1],col='blue',bg='blue',pch=2)

line_length <- 50.0 # how long to draw the hyperplanes
p1 <- -c(ztrue[3]/ztrue[1],0)    # start point of line
p2 <- -c(0, ztrue[3]/ztrue[2])   # end point of line
vdir <- (p2-p1)/(sum((p2-p1)**2)**0.5)  # the line vector
f  <- line_length / sum((p2-p1)**2)**0.5 # scaling the line
p2 <- p2+vdir*f            # final end point
p1 <- p1-vdir*f            # final start point
lines(c(p1[1],p2[1]),c(p1[2],p2[2]),col='green', lwd=2)

line_length <- 50.0 # how long to draw the hyperplanes
p1 <- -c(zf[3]/zf[1],0)    # start point of line
p2 <- -c(0, zf[3]/zf[2])   # end point of line
vdir <- (p2-p1)/(sum((p2-p1)**2)**0.5)  # the line vector
f  <- line_length / sum((p2-p1)**2)**0.5 # scaling the line
p2 <- p2+vdir*f            # final end point
p1 <- p1-vdir*f            # final start point
lines(c(p1[1],p2[1]),c(p1[2],p2[2]),col='black', lwd=2)


legend_string <- c('z used to generate data', 'Perceptron final z')
legend(x="bottomleft",legend_string,cex=1.5,col=c('green','black'),pch=c('-','-'))
legend(x="topright",c("Test Data: -1","Test Data: +1"),cex=1.5,col=c('blue','red'),pch=c(2,0))
dev.off()

