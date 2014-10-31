import numpy as np

def crrate(ratio, x, nh, T):

    T2 = T / 100.0
    xe = x
    
    
    k169 = 4.1E-8*T2**(-0.52)
    k1610 = 7.7E-8*T2**(-0.52)
    
    phi = (1.0 - xe/1.2)*((0.67)/(1.0+ xe/0.05))
    
    CR = (k169 + k1610)*nh*xe*0.5*ratio/(1.0+phi)
    
    return CR
    
    
rmin = 1.0E-9
rmax = 1.0E-7

x = 1.5E-4
nh = 200.0
T = 70.0

crmin, crmax = crrate(rmin,x,nh,T), crrate(rmax,x,nh,T)

print "crmin - crmax"
print "%5.4e %5.4e"%(crmin,crmax)
