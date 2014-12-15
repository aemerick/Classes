import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
plt.rc('font',size=16)

kb = 1.3806E-16
c  = 2.9979E10
me = 9.109E-28
Q  = 13.6 *  1.602E-12  # ev -> erg

def neutral_fraction(T,eta):


    consts = 3.84 * eta * (kb*T/(me*c*c))**(3./2.) * np.exp(Q/(kb*T))

    x = (-1.0 + (1.0+4.0*consts)**0.5) / (2.0*consts)       
        
    return x
    
    
T = np.linspace(2200.0, 5300.0, 1000.0)
etas = [1.0E-9, 1.0E-10, 5.5E-10]
lw = 1.5

for eta in etas:
    x = neutral_fraction(T,eta)
    plt.plot(T, neutral_fraction(T,eta), label = r'$\eta$ = %3.2e'%(eta), lw=lw)
    

plt.plot([3740.0,3740.0],plt.ylim(),label = r'T$_{rec}$', color='black', ls = '-')
plt.plot([3000.0,3000.0],plt.ylim(),label = r'T$_{dec}$ = T$_{ls}$', color='black', ls = '--')
plt.minorticks_on()
plt.legend(loc='best',fancybox=True)
plt.savefig('x_T.png')
plt.close()



print neutral_fraction(3000.0,5.5E-10)




