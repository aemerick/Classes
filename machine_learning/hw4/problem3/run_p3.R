#

K = 3
tau = 10.0

m <- MultinomialEM(H,K,tau)
hist(m, breaks=c(0,1,2,3,4,5),freq=FALSE)

