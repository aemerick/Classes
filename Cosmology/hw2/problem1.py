import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

plt.rc('font',size=16)


H_o = 1.0
q_o = -0.55

#def t_e_multi(z):

omega_r = 8.4E-5
omega_m = 0.3
omega_l = 0.7
omega_o = omega_r + omega_m + omega_l



#

def dp_to(z, universe):
    """
    Calculates the proper distance at t = t_o as a function of z
    for a matter dominated and dark matter dominated universe using 
    the single component universe formula by setting w
    
    For multi component, integrates the acceleration equation given 
    omega_rad, omega_matter, omega_l, and omega_o
    """


    if universe == 'matter':
        w = 0.0   
#        dp = 2.0 * ( 1.0 - 1.0 / (1.0+z)**0.5)
    
    elif universe == 'lambda':
        w = -1.0
    
    elif universe == 'multi':
        #return z*(1.0 - z*(1.0 - 0.55)/2.0)
        
        #return np.zeros(len(z))
        
        func = lambda x: (omega_r/x**2 + omega_m/x + omega_l*x*x + (1.0-omega_o))**(-0.5) /x
        
        dp = np.zeros(len(z))
        for i in range(len(z)):
            dp[i] = integrate.quad(func,1.0/(1.0+z[i]),1.0)[0]
    
        return dp
        
    
    dp = 2.0/(1.0+3.0*w) * (1.0 - (1.0+z)**(-0.5*(1.0+3.0*w)))
    
    return dp
    
    
def dp_te(z, universe):
    """
    Calculates proper distance at t_e as a function of z. For single component,
    returns the easy analytic formulas.
    
    For multi, integrates accel equation to calculate dp_to and uses
    dp_te = dp_to / (1+z)
    """

    if universe == 'matter':
        dp = 2.0/(1.0+z) * (1.0 - 1.0/(1.0+z)**0.5)
   #     w = 0.0
    elif universe == 'lambda':
        dp = z/(1.0+z)
    
    elif universe == 'multi':
        return dp_to(z, universe) / (1.0+z)
    
    
    
    return dp
    
def te(z,universe):
    """
    Calculates t_e as a function of z.
    Does not work for 'lambda'
    
    For multi, integrates acceleration equation
    """

    if universe == 'matter':
#        te = 2.0/(3.0*H_o) * 1.0/(1.0+z)**3
        w = 0.0   
    elif universe == 'lambda':
#        w = -1.0
        return np.zeros(len(z))
    elif universe == 'multi':
       # return (1.0+q_o/2.0)*z*z/H_o
       


    
        func = lambda x : (omega_r/x**2 + omega_m/x + omega_l*x*x + (1.0-omega_o))**(-0.5)
    
        te = np.zeros(len(z))
        for i in range(len(z)):
            te[i] = integrate.quad(func,0.0,1.0/(1.0+z[i]))[0]
    
        return te
       
    
    te = 2.0 / (3.0*(1.0+w)*H_o) * (1.0/(1+z)**(3.0*(1+w)/2.0))
    
    return te


def to(z,universe):
    """
    calculates t_o as a function of z for 'matter' using single component
    universe formula, setting w. This does not work for 'lambda'
    
    For multi, t_o = 1/H_o
    """ 
   


    if universe == 'matter':
#        to = 2.0 / (3.0*H_o)
        w = 0.0        
    elif universe == 'lambda':
        return np.zeros(len(z))
    
    elif universe == 'multi':
        return 1.0/H_o
    
    to = 2.0 / (H_o * 3.0 * (1.0+w))
    
    return to
    


def l(z):
    """
    Calculates l for finding the angular diameter distance.
    
    l is the physical size of the horizon on the sky
    """

    func = lambda x: (omega_r/x**2 + omega_m/x + omega_l*x*x + (1.0-omega_o))**(-0.5) /x
        
    l = np.zeros(len(z))
    for i in range(len(z)):
        l[i] = integrate.quad(func,0.0,1.0/(1.0+z[i]))[0]
    
   
    return l / (1.0+z) 
 


#-------------------------------------------------------------------------------
# Plot dp_to for the 3 universes
#
types = ['matter','lambda','multi']
z = np.logspace(-2,3,1E3)
lw = 1.5

for t in types:
    plt.plot(z,dp_to(z,t),label=t,lw=lw)
    
plt.xlabel('z')
plt.ylabel(r'd$_p$(t$_o$)H$_o$/c')
plt.loglog()
plt.xlim(0.01,1000.0)
plt.ylim(0.001,1000.0)
plt.legend(loc='best',fancybox=True)
plt.savefig('dp_to.png')
plt.close()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot t_e for the three universes
#
for t in types:
    plt.plot(z,dp_te(z,t),label=t,lw=lw)
    
plt.xlabel('z')
plt.ylabel(r'd$_p$(t$_e$)H$_o$/c')
plt.loglog()
plt.legend(loc='best',fancybox=True)
plt.xlim(0.01,1000.0)
plt.ylim(0.001,1000.0)
plt.savefig('dp_te.png')
plt.close()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot t_o - t_e
#
z = np.linspace(0.0,6.0,1000.0)
for t in types:
    if t == 'lambda': # calculates this for lambda universe
        plt.plot(z,-1.0*np.log(1.0/(1.0+z)),label=t,lw=lw)
    else:
        plt.plot(z,(to(z,t) - te(z,t)), label=t,lw=lw)
    

plt.xlabel('z')
plt.ylabel(r'H$_o$(t$_o$ - t$_e$)')

plt.legend(loc='best',fancybox=True)
plt.xlim(0.0,6.0)
plt.savefig('Ho_to-te.png')
plt.close()

z = np.linspace(0.0,6.0,100.0)
for t in types:
    if t == 'lambda': # calculates this for lambda universe
        plt.plot(z,-1.0*np.log(1.0/(1.0+z))*14.0,label=t,lw=lw)
    elif t == 'multi':
        plt.plot(z,(to(z,t) - te(z,t))*14.0 - 0.5, label=t,lw=lw)
    else:
        plt.plot(z,(to(z,t) - te(z,t))*14.0, label=t,lw=lw)
    
    
   

plt.xlabel('z')
plt.ylabel(r't$_o$ - t$_e$')

plt.legend(loc='best',fancybox=True)
plt.xlim(0.0,6.0)
plt.ylim(0.0,18.0)
plt.savefig('to-te.png')
plt.close()
#-------------------------------------------------------------------------------  



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Calculate the angular size of the horizon     
z = np.array([1100.0])
l = l(z)                                # using function defined above
d = dp_to(z, 'multi') / (1.0 + z)       # d_te = dp_to / (1+z)

# print out the answer
print 'theta = %5.4e'%(l/d)
print 'theta = %5.4e'%(180.0*l/d / np.pi)



