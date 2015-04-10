# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 4
#
#  This file contains the code used to run and plot the results
#  from the clustering algorithm in problem 3 of homework 4
#

source("p3.R")

# choose a K and a tau
K <- 2
tau <- 1.0E-10

# load in H
H <- matrix(readBin("histograms.bin", "double", 640000), 40000, 16)

# run the algorithm
m <- MultinomialEM(H,K,tau)

n <- length(m)
frac <- c(rep(0,K))

for (k in 1:K){
    print(noquote(sprintf("k = %.1f - %.3f",k,length(m[m==k])/n)))
    frac[k] = length(m[m==k])/n
}
write(frac, noquote(sprintf("k%i.out",K)),ncolumns=1)
matrix_m <- matrix(m,nrow=200,ncol=200)

# function to rotate image matrix by 90 degrees
rotate_function <- function(x) t(x)[,nrow(x):1]

# plot and save
pdf(noquote(sprintf("k%i.pdf",K)))
image(rotate_function(t(matrix_m)),col=gray( seq(0,1,1/(K-1) )))
dev.off()
graphics.off()
