# Author: Andrew Emerick - aje2123@columbia.edu
# UNI   : aje2123
# Course: STATW4400
# HW    : HW 2

classify <- function(S,z){
# 
# f(x) = sgn( <(1,x),z>) as defined on slide 35 in the lecture
#
   y = c(sign(S %*% z)) 

   return(y)
} # end of function classify

