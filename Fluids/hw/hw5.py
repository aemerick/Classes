import matplotlib.pyplot as plt
import numpy as np

r = 20.0 # kpc
L = np.linspace(-4.0,4.0,10000.0)
L = 10**(L)
L = 10**(1.0)
a = 3.8 # unitless
v = 3500.0 # km/s


gcm = 0.15*(1000.0*100.0)*v*v/(r/3.2E-17) /a 
# convert to Myr kpcl
v = v * 3.24E-17
v = v * 60 * 60 * 24 * 365 * 1.0E6

g = 0.15 * v*v/r /a 

tau = L/(2.0*np.pi) * (abs(2.0*np.pi * (1.0+a)/(1.0-a) * (1.0/(g*L))))**(0.5)

tau = np.log10(abs(tau))
L = np.log10(L)
plt.plot(L,tau)
plt.show()
