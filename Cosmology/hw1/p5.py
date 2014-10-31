import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def RHS(sign_of_derivative):
    """
    This function is the right hand side of the differential equation da/dt = 
    RHS. However, this function returns a function corresponding to whether or
    not I am working in dtau, where tau = H_o*t, or ds, where s = t_o - t. The
    latter just is the normal equation times -1.0.
    
    Parameters
    ------------
    sign_of_derivative : integer
        Supply > 0 if working in dtau or < 0 if working in ds.
        
    Returns
    -------
    function :
        returns the right hand side of the differential equation with/without
        a negative sign depending on the dependent variable used
    """
    # definitions of the omega's
    omega_r = 8.4E-5
    omega_m = 0.3
    omega_l = 0.7
    omega_o = omega_r + omega_m + omega_l 
    
    # the "normal" function of da / dtau = 
    def positive(a,t):
        return np.sqrt( (omega_r /(a*a) + omega_m  / a + omega_l*a*a -(1.0-omega_o)) )  
    
    # the rewritten function of da / ds = 
    def negative(a,t):
        return -1.0 * np.sqrt( (omega_r /(a*a) + omega_m  / a + omega_l*a*a - (1.0-omega_o)) )

    # depending on dtau or ds, return one of the above functions
    if sign_of_derivative < 0:
        return negative
    else:
        return positive
      
          
          
def solve_for_a(t0, a0, final_time, 
                tf_shift = 10.0, number_of_samples = 1E3, ode_function = RHS(1.0),
                loop_max = 200):
    """
    This function takes in an initial time and value for a and solves the
    differential equation beginning at t0 and ending at final_time. The function
    scales the time step size by computing the solution in multiplicative 
    intervales of tf_shift using number_of_samples points in each interval. This
    is essentially equivalent to uniform dt in logbase tf_shift space (with
    the defaults, logbase 10). 
    
    Parameters
    -----------
    t0: float
        Initial time 
        
    a0: float
        initial value for the scale factor
        
    final_time: float    
        Final time to compute to
        
    tf_shift: optional, float
        Computes solution from t0 to final_time by over a time range in factors
        of tf_shift. Default = 10.0
    
    number_of_samples: optional, float
        Number of sample points within the tf_shift range. Default 1000.0
        
    ode_function: optional, function
        RHS function. Default = RHS(1.0), or da / dt
        
    loop_max: optional, float
        Maximum number of loops to get from t0 to final time. here to prevent
        locking up the computer when testing. Default = 20.0
    
    Returns
    --------
    tsolved : array
        Array of all of the t values computed in solution
    
    asolved : array
        Array of all of the scale factor values computed
    """

    # if statements on setting initial final time as scaled by tf_shift
    # depending on whether or not moving forwards or backwards in time
    if t0 < final_time and t0 > 0:
        increasing = True
        tf = t0 * tf_shift
    elif t0 < final_time and t0 == 0:
        increasing = True
        tf = t0 + 0.000001
    elif t0 < final_time and t0 < 0:
        increasing = True
        tf = t0 / tf_shift
        sign_of_t = -1.0
    else:
        increasing = False
        tf = t0 / tf_shift

    tout = []
    aout = []
    j = 0
    
    current_time = t0
    
    # loop until current time is equal to final time. 
    # boolean clauses are there because may be solving towards larger 
    # or smaller t (t0 < final_time or t0 > final_time)
    while ((current_time < final_time and increasing) or\
          ( current_time > final_time and not increasing)) and\
          j < loop_max:

        # t values to sample
        t = np.linspace(t0, tf, number_of_samples)
        
        # solve!
        a = odeint(ode_function, a0, t)

        # set the next final time down by a factor of 10
        # set to final_time if this moves tf past final_time
        if increasing and tf > 0:
            tf = np.min([tf*tf_shift, final_time])
        elif tf > 0:
            tf = np.max([tf/tf_shift, final_time])
        else:
            tf = np.max([tf / tf_shift, final_time])
    
        # take the last a and t as the boundary conditions for
        # the next loop. Round them to some number of sig figs
        # solver got upset when bc's were near machine precision
        a0 = float('%.8g' % a[-1]) 
        t0 = float('%.8g' % t[-1])
       
        current_time = t0   
        
        # append the found a values and sample points to a list
        aout = aout + list(a)
        tout = tout + list(t)
        
        j = j + 1 
        
    return np.array(tout), np.array(aout)

#-------------------------------------------------------------------------------
# try this first using the defined boundary conditions

#a0, t0 = 1.0, 1.0
#tf     = 1.0E-10
a0, t0 = 1.0E-8, 1.0E-12
tf = 10.0

tsolved, asolved = solve_for_a(t0, a0, tf, tf_shift = 2.0,
                               ode_function= RHS(1.0))
plt.plot(tsolved, asolved, color = 'blue', lw = 2.0,ls = '-')

#-------------------------------------------------------------------------------
# try this again using the defined boundary conditions
# but recast axis to s = t0 - t so we are moving up in 's'
#a0, t0 = 1.0, 1.0
#tf     = 1E-10

# recast in s = t0 - t
#sf = t0 - tf
#s0 = 0.0

# using RHS(-1.0) which changes the sign in the equation since ds = - dt
#tsolved, asolved = solve_for_a(s0, a0, sf, ode_function= RHS(-1.0))
#plt.plot(t0 - tsolved, asolved, color = 'black', lw = 2.0, ls = '-',
#          label='forwards in: s = t0 - t')

#-------------------------------------------------------------------------------
# try going the other direction from the start points
#a0, t0 = 1.0, 1.0
#tf     = 10.0

#tsolved, asolved = solve_for_a(t0, a0, tf, ode_function= RHS(1.0))
#plt.plot(tsolved, asolved, color = 'blue', lw = 2.0, ls = '-.',
#          label='forwards from a0,t0=1.0,1.0')

#-------------------------------------------------------------------------------
# try going the other direction from a lower point
#a0, t0 = 2.8E-4, 4.7E4 / 14.0E9
#tf     = 10.0

#tsolved, asolved = solve_for_a(t0, a0, tf)
#plt.plot(tsolved, asolved, color = 'red', lw = 2.0, ls = ":", label='forwards')

#-------------------------------------------------------------------------------
# try going the other direction from a lower point
#a0, t0 = 2.8E-4, 4.7E4 / 14.0E9
#tf     = 1.0E-10

#tsolved, asolved = solve_for_a(t0, a0, tf)
#plt.plot(tsolved, asolved, color = 'red', lw = 2.0, ls = "-", label='backwards')

#-------------------------------------------------------------------------------
# formatting
plt.loglog()
ymin, ymax = 1.0E-6, 1000
xmin, xmax = 1.0E-10, 10
plt.plot([1.0,1.0],[ymin,ymax],lw=1.5,ls=':',color='black')
plt.plot([4.7E4/14.0E9,4.74E4/14.0E9],[ymin,ymax],lw=1.5,ls=":",color='black')
plt.plot([9.8E9/14.0E9,9.8E9/14.0E9],[ymin,ymax],lw=1.5,ls=":",color = 'green')
plt.plot([xmin,xmax],[1.0,1.0],lw=1.5,ls=":",color='black')
plt.plot([xmin,xmax],[2.8E-4,2.8E-4],lw=1.5,ls=":",color='black')
plt.plot([xmin,xmax],[0.75,0.75],lw=1.5,ls=":",color='green')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel(r"scale factor, a")
plt.xlabel(r"H_o t")
plt.legend(loc='best',fancybox=True)
plt.savefig('p5.png')
