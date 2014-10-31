"""
    1D Fluid Equations evolved with the Roe type, linearized Riemann solver...
    Author: Andrew Emerick
"""
import numpy as np
import copy

# number of grid points and xrange
N    =  100
xmin, xmax = -10.0, 10.0

tstepmax = 50 # maximum time step... not time

outarray = np.arange(0,1000,1) # array list of time steps to write out a dump

# some globals used as switches for debugging purposes
FIXED_DT = [True, 0.05]
EPSILON_FIX = [False, 1.0E-6]
FLUX_LIM = 'Lax-Wendroff'
INIT_COND = 'sod1'
ADD_VISC = False
BOUNDARY_CONDITION = 'constant' # uses initial values at boundaries
# Slope limiter options:
# minmod
# Lax-Wendroff
# beam-warming
# superbee
# fromm  

class simGrid:
    """
    Object defines the grid. The object contains knowledge of the boundary
    conditions being used.

    Initialized with the number of bins (i.e. cell centers), not the number 
    of interfaces.... this is a centered cell implementation...

    """

    def __init__(self, N):
        self.N  = N
                
    def setLim(self,xmin,xmax):
        """
        set the xmin and xmax of the grid
        """
        self.xmin = xmin
        self.xmax = xmax
        
    def setBoundaryType(self,boundary):
        """
        Set the boundary type. Options are:
        
        constant -> fixed at starting values
        periodic -> self explanatory
        """
        self.bc = boundary
        
    def makeGrid(self):
        """
        Makes the grid given xmin, xmax, and N
        """
        try:
            xmin = self.xmin
            xmax = self.xmax
        except:
            print "must set xmin and xmax first via  grid.setLim(xmin,xmax)"
        
        self.grid = 1.0*np.linspace(xmin,xmax,self.N + 1)
        self.dx   = self.grid[1] - self.grid[0]
        
    def center(self, index=None):
        """
        returns array with coordinates of bin centers
        
        only return designated bin center if index is provided
        
        returns array of size of self.grid, or a single value
        """
                      
        if index == None:
            center = 0.5*(self.grid[1:] + self.grid[:-1])
        else:
            center = 0.5*(self.grid[index + 1] + self.grid[index])
            
        return center



def set_initial_conditions(state, name, inputfile=None):
    """
    Takes in simstate class and a initial condition type
       
    Specify name of initial conditions to be set. Alternatively, provide a 
    previous dump, and the simState is set to that dump.
    
    gamma is currently fixed at 7.0/5.0 ..... 
    """ 
    N = state.grid.N
    xmin = state.grid.xmin
    xmax = state.grid.xmax
    c = state.grid.center()
    cgrid = 0.5*(xmax+xmin)
        
    rho = np.zeros(N)
    u   = np.zeros(N)
    P   = np.zeros(N)
        
        
    if not inputfile == None:    
        data = np.genfromtxt(inputfile,names=True)
        
        rho = data['Density']
        P   = data['Pressure']
        u   = data['Velocity']
        gamma = 7.0/5.0
        
    elif name == 'sod1':
        # test case from computational gasdynamics book..... using whatever units
        # given in there.... 
        L = c < cgrid
        R = (c > cgrid) + (c==cgrid)
            
        sizeL, sizeR = np.size(L[L]), np.size(R[R])

        
        rho[L], rho[R] = 1.0E5 * np.ones(sizeL), 1.25E4 * np.ones(sizeR)  
        u[L]  ,   u[R] = 0.0   * np.ones(sizeL), 0.0    * np.ones(sizeR)
        P[L]  ,   P[R] = 1.0   * np.ones(sizeL), 0.1    * np.ones(sizeR)
        gamma = 7.0 / 5.0
            
            
    elif name == 'linear':
        rho[0:0.25*N] = np.ones(0.25*N)*1.0E5
        rho[0.25*N:0.75*N] = np.linspace(1.0E5,1.0E4,0.5*N)
        rho[0.75*N:] = np.ones(0.25*N)*1.0E4
        u   = np.zeros(N)
        P[0:0.25*N] = np.ones(0.25*N)*1.0
        P[0.25*N:0.75*N] = np.linspace(1.0,0.1,0.5*N)
        P[0.75*N:] = np.ones(0.25*N)*0.1
        gamma = 7.0 / 5.0
    
    elif name == 'curve':
        # sinusoidal wave taking up 0.5 of the grid space.. centered at grid 
        # center... rest of box is a constant... rho, P, and u are continious
    
        l = np.pi*(2.0/(N))
        
        rho[0:0.25*N ] = np.ones(0.25*N)*3.0E4
        rho[0.75*N:]   = np.ones(0.25*N)*3.0E4
        rho[0.25*N:0.75*N] = (2.0E4-1.0E4)* np.sin( l*np.arange(0.5*N)) + 3.0E4
        
        P[0:0.25*N] = np.ones(0.25*N)*1.5
        P[0.75*N:]  = np.ones(0.25*N)*1.5
        P[0.25*N:0.75*N] = (1.0-0.5)    * np.sin( l*np.arange(0.5*N)) + 1.5
        
        u[0:0.25*N] = np.ones(0.25*N)*0.01
        u[0.75*N:]  = np.ones(0.25*N)*0.01
        u[0.25*N:0.75*N] = (0.1-0.01)   * np.sin( l*np.arange(0.5*N)) + 0.01
        
        l = np.pi*(4.0/(1.0*N))
        u = np.ones(N)*0.00
        gamma = 0.75
    

    # set the simstate time and initial conditions....    
    state.setStateFromPrim(rho,u,P,gamma)
    state.time(0.0)
        
        
def arrayList(N, ndim=1):
    """
    Returns a 3 x N array... set to zero....
    I made a function to do this since it comes up a lot in order to generate
    a vector quantity that is defined at every grid cell. The number of 
    dimensions is default to 1 since this is a 1D code.
    """

    return np.array([np.zeros(N) for i in range(ndim*3)])


class simState:
    """
    This object defines the simulation at a given point in time....
    Nothing is requested at instantiation... 
    
    do simState.initializeSimulation(grid) where grid is a grid object in order
    to initialize things to zero
    
    """

    def __init__(self):
        print "initalized"
        
    def initializeSimulation(self, grid):
        """
        simState object stores the state vector and flux vector at all cell 
        centers, as well as the grid object.
        """
        self.q = arrayList(grid.N)
        self.f = arrayList(grid.N)
        self.grid = grid
        


    def time(self, t=None):
        """
        Sets or returns the current time depending on whether or not t 
        is passed to the current simState
        """
    
        if t == None:
            try:
                return self.t
            except:
                print "NO TIME ASSOCIATED WITH THIS SIMULATION STATE"
        else:
            self.t = t
    
    
    def add_visc(self):
        """
        Adds von Neumann-Richtmyer artificial viscosity. This adds a bulk
        viscosity that acts as an additional pressure term.... i.e. 
        P -> P + PI, where PI is the 'visc' pressure term calculated below.
        
        xi is a freely chosen value, usually 2 or 3 which is a tuning 
        parameter that determines the spread of the shock
        """
    
        xi = 3.0
    
        rho = self.getPrimitive('Density')
        u = self.getPrimitive('Velocity')
        P = self.getPrimitive('Pressure')
        
        visc = np.zeros(np.size(P))
        
        i = 1
        while i < np.size(P) - 1:
            if u[i+1] < u[i-1]:
                visc[i] = 0.25*xi*xi*(u[i+1] - u[i-1])**2 * rho[i]
            else:
                visc[i]=0.0
            
            i = i + 1
        
        visc[0] = 0.0
        visc[-1] = 0.0
            
        
        # modify the state and flux vectors with the new pressure
        self.setStateFromPrim(rho, u, P + visc, self.gamma)
    
    
    
    def setStateFromPrim(self, rho, u, P, gamma):
        """
        Given arrays of density, velocity, and pressure at all points,
        sets the state vector and flux. ... also gamma
        """
               
        # set gamma
        self.gamma = gamma
     
        # energy and total energy
        e = P/((self.gamma - 1.0)*rho)
        etot = e + 0.5 * u * u
              
        # total enthalpy
        htot = etot + P/rho
        
        # set the state vector q        
        self.q[0] = rho
        self.q[1] = rho * u
        self.q[2] = rho * etot

        # set the flux vector
        self.setFlux()

        # make sure boundary conditions are enforced
        try:
            t = self.t
            
            self.setBoundaryCondition()
        except:
            self.t = 0
            self.setBoundaryCondition()

    def setFlux(self):
        """
        Set the flux vector using the state vector q. 
        """   
        self.f[0] = self.q[1]
        
        self.f[1] = (self.gamma-1.0)*self.q[2] +\
                    0.5*(3.0-self.gamma)*(self.q[1]**2)/self.q[0]

                    
        self.f[2] = self.gamma*self.q[2]*self.q[1]/self.q[0] +\
                    0.5*(1.0-self.gamma)*((self.q[1])**3)/(self.q[0])**2

    def getPrimitive(self, name, q = None):
        """
        Calculate a primitive variable from the state vector
        
        can provide a separate q vector to calculate primitives from... this
        is used really only for the boundaries
        """

        if q == None:
            q = self.q

        
        # define functions to calculate a given primitive variable 
        def _density():
                 
            return q[0]
            
        def _velocity():
            return q[1] / q[0]
        
        def _pressure():
            P = (q[2] - 0.5*(q[1]*q[1])/q[0])*(self.gamma - 1.0)

            # would be really nice to not have to do this........
            P[P < 0] = np.zeros(np.size( P[P<0] ) )
            return P
            
        def _etot():
            return q[2] / q[0]
        
        def _htot():
            return (q[2]/q[0] + _pressure()/q[0])   
        
        def _cs():
            return (self.gamma * _pressure() / _density())**0.5
        
        
        primDict = {'Density': _density,
                    'Pressure': _pressure,
                    'Velocity': _velocity,
                    'etot': _etot,
                    'htot': _htot,
                    'cs': _cs
                    }
                    
       
          
       

        return primDict[name]()
        


    def writeOut(self, outname='',include_state=False):
        """
        Function writes out to file... only doing primitive variables for now.
        rho, u, p, maybe tack on e and h ....
        
        This needs to be nice, but for debugging purposes, only doing to 
        write out to ascii for now... just ot be quick and easy to focus on
        coding rather than fancy outputting.....
        """
        x   = self.grid.center()
        rho = self.getPrimitive('Density')
        u   = self.getPrimitive('Velocity')
        P   = self.getPrimitive('Pressure')
        
        if include_state:
            data = np.column_stack((x,rho,u,P,self.q[0],self.q[1],self.q[2]))
            header = '# x Density Velocity Pressure q0 q1 q2'
            
        else:
            data = np.column_stack((x,rho,u,P))
            header = '# x Density Velocity Pressure'
    
        
        np.savetxt(outname + '_simstate.txt', data, header=header, fmt='%1.4e')
        

        
    def getInterfaceValue(self, name, direction):
        """
        Get the value at the i - 1/2 interface... direction is supplied to 
        either obtain the right or left value to the interface...
        
        A primitive variable can be supplied for 'name' or a kth element of the
        state or flux vector. (i.e. q0,q1,q2,f0,f1,f2).
        
        This function also implements the specified boundary conditions and 
        reads from 'ghost cells' when needed
        """

        try: # see if it is a primitive variable... if fail, go to except
            variable = self.getPrimitive(name)
        except:

            # if not primitive, check if one of the state vectors or flux 
            if 'q' in name:

                for i in range(3):
                    if str(i) in name:
                        index = i


                variable = self.q[index]

            elif 'f' in name:

                for i in range(3):
                    if str(i) in name:
                        index = i

                variable = self.f[index]
        
        iVar = np.zeros( np.size(variable) + 1 ) # N + 1 interfaces
        
        # apply the boundary conditions
        if direction == "R":
            iVar[0:-1] = variable
            iVar[-1]   = self.applyBC(name,'right')
            
        elif direction == "L":
            iVar[1:] = variable
            iVar[0]  = self.applyBC(name, 'left')
            
            

        return iVar
                
    def applyBC(self, name, side):
        """
        retrieves the value at the boundary... i.e. the "ghost" cell...
        """
        if self.grid.bc == 'constant':

            if 'q' in name:
                                  
                for i in range(3):
                    if str(i) in name:
                        index = i
                
                if side == 'left':
                    bVal = self.__qL__[index]
                elif side == 'right':
                    bVal = self.__qR__[index] 

            elif 'f' in name:
                for i in range(3):
                    if str(i) in name:
                        index = i
                        
                if side == 'left':
                    bVal = self.__fL__[index]
                elif side == 'right':
                    bVal = self.__fR__[index]                  

            else:         
                if side == 'left':
                    bVal = self.getPrimitive(name, q = self.__qL__)
                elif side == 'right':
                    bVal = self.getPrimitive(name, q = self.__qR__)

  
  
            return bVal
            
    

    def setBoundaryCondition(self):
        """
        Enacts the specified boundary condition when called. If 'constant' is
        the specified conditon, then this function does nothing unless the time
        of the simulation object is 0.0 .... otherwise it will overwrite the 
        fixed ghost cells
        """
        
    
        if self.grid.bc == 'constant' and self.t == 0.0:
            # conditions are fixed to their starting values at edges
            self.__qR__ = np.array([[self.q[0][-1]],[self.q[1][-1]],[self.q[2][-1]]])
            self.__qL__ = np.array([[self.q[0][0]] ,[self.q[1][0]] ,[self.q[2][0]]])
            
            self.__fR__ = np.array([[self.f[0][-1]],[self.f[1][-1]],[self.f[2][-1]]])
            self.__fL__ = np.array([[self.f[0][0]] ,[self.f[1][0]] ,[self.f[2][0]]])
            
                                 
        elif self.grid.bc == 'periodic':
            self.__qR__ = np.array([[self.q[0][0]],[self.q[1][0]],[self.q[2][0]]])
            self.__qL__ = np.array([[self.q[0][-1]],[self.q[1][-1]],[self.q[2][-1]]])

            self.__fR__ = np.array([[self.f[0][0]],[self.f[1][0]],[self.f[2][0]]])
            self.__fL__ = np.array([[self.f[0][-1]],[self.f[1][-1]],[self.f[2][-1]]])
            
        elif not self.grid.bc == 'constant':
            print "nothing set with boundary conditions... check bc settings"       
                   
            
            

def roe_evolve(simPrev):
    """
    This is the meat and potatoes of the program. This evolves the simulation
    state according to the roe algorimth... welll.... it should....
    """
    C = 0.5 # fixed CFL condition

    N = simPrev.grid.N
    gamma = simPrev.gamma
    dx = np.abs(simPrev.grid.dx)

    sim = copy.deepcopy(simPrev)

    # add artificial viscosity (i.e. modify pressure term accordingly)
    if ADD_VISC:
        simPrev.add_visc()

    # Calculate the roe averages at the cell interfaces
    rhoL = simPrev.getInterfaceValue('Density','L')
    rhoR = simPrev.getInterfaceValue('Density','R')
    

    rho_roe = (rhoL * rhoR)**0.5  # Density

    denom = (rhoL)**0.5 + (rhoR)**0.5  # denom. for normalization
    
    # velocity
    u_roe = ( (rhoL**0.5)*simPrev.getInterfaceValue('Velocity','L') +\
              (rhoR**0.5)*simPrev.getInterfaceValue('Velocity','R') ) / denom
    
    # total enthalpy
    htot_roe = ((rhoL**0.5)*simPrev.getInterfaceValue('htot','L') +\
                (rhoR**0.5)*simPrev.getInterfaceValue('htot','R') ) / denom
   
    # sound speed.... calc from roe velocity and enthalpy
    Cs_roe = ((gamma-1.0) * (htot_roe - 0.5*u_roe*u_roe) )**0.5
   
    # make sure dt satisfies cfl... adjust if not
    cfldt = C * simPrev.grid.dx / (np.max(u_roe) + np.max(Cs_roe))
    
    dt = FIXED_DT[1]
   
    if dt > cfldt:
        dt = cfldt

    # set up eigenvalues and eigenvectors
    eigenval = arrayList(N + 1)
    e1 = arrayList(N + 1) 
    e2 = arrayList(N + 1)
    e3 = arrayList(N + 1)
    
    # calculate the eigenvalues at the interface
    eigenval[0] =  u_roe - Cs_roe
    eigenval[1] =  u_roe
    eigenval[2] =  u_roe + Cs_roe
    
    # turn on some fix
    if EPSILON_FIX[0]:
        epsilon = EPSILON_FIX[1]
        for i in [0,2]:
            selection = np.abs(eigenval[i]) < epsilon
            eigenval[i][selection] = 0.5*(eigenval[i][selection]**2 / epsilon + epsilon)
    
    # calculate the eigenvectors
    e1[0] = 1.0*np.ones(N+1)
    e1[1] = u_roe - Cs_roe
    e1[2] = htot_roe - Cs_roe * u_roe
    
    e2[0] = 1.0*np.ones(N+1)
    e2[1] = u_roe
    e2[2] = 0.5*u_roe*u_roe
    
    e3[0] = 1.0*np.ones(N+1)
    e3[1] = u_roe + Cs_roe
    e3[2] = htot_roe + Cs_roe * u_roe

    # ----------------------------

    
    dq_roe = arrayList(N+1)
    dq     = arrayList(N+1)
    
    # calculate the jump of the state vector across the interface
    for i in range(3):
        dq[i]  = simPrev.getInterfaceValue('q'+str(i), 'R') -\
                                               simPrev.getInterfaceValue('q'+str(i), 'L')
 

   
    # some things used to calculate the modified dq's
    A = (gamma - 1.0) / (2.0 * Cs_roe * Cs_roe)
    B = (dq[1] - u_roe * dq[0]) / (2.0 * Cs_roe)
    ksi = u_roe * dq[1] - dq[2]
    ekin = 0.5 * u_roe * u_roe
    
    # calculate the modified dq's    
    dq_roe[0] = A * (ekin * dq[0] - ksi) - B
    dq_roe[1] = A * ((htot_roe - 2.0*ekin)*dq[0] + ksi)
    dq_roe[2] = A * (ekin * dq[0] - ksi) + B
    

    dq_roe_tilde = arrayList(N+1)
    for i in range(3):
        dq_roe_tilde[0] += dq_roe[i] * e1[i]
        dq_roe_tilde[1] += dq_roe[i] * e2[i]
        dq_roe_tilde[2] += dq_roe[i] * e3[i]
    
    

    # make theta, the sign of the eigenvalues
    theta = arrayList(N+1)
    for i in range(3):
        theta[i] = 1.0*np.sign(eigenval[i])


    # calculate epsilon, the eigenval * (dt/dx)
    eps = arrayList(N+1)
    eps = eigenval * dt / (1.0*dx)
    
    # calculate phi, as determined by the specified slope limiter method
    phi = slope_limiter(simPrev, eigenval, name=FLUX_LIM)
     
    # calculate the averaged flux across the interface 
    favg = arrayList(N+1)
    favg[0] = (simPrev.getInterfaceValue('f0','R') +\
                                          simPrev.getInterfaceValue('f0', 'L'))
    favg[1] = (simPrev.getInterfaceValue('f1','R') +\
                                          simPrev.getInterfaceValue('f1', 'L'))
    favg[2] = (simPrev.getInterfaceValue('f2','R') +\
                                          simPrev.getInterfaceValue('f2', 'L'))
    
    favg = favg * 0.5
    
    # loop to calculate the correction term
    sigma = arrayList(N+1)

    
    for i in range(3):
        sigma[0] += eigenval[i] * dq_roe[i] * e1[i] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
        sigma[1] += eigenval[i] * dq_roe[i] * e2[i] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
        sigma[2] += eigenval[i] * dq_roe[i] * e3[i] * (theta[i] + phi[i] * (eps[i] - theta[i]) )

        
        
    
    # apply correction to get    f^n+1/2
    fnew = favg - 0.5*sigma
    
    # update the state vector with the right and left interface fluxes
    for i in range(3):
        sim.q[i] = simPrev.q[i] - (dt/(1.0*dx))*(fnew[i][1:] - fnew[i][0:-1]) 

    # update flux within each cell and set the new time  
    sim.time( simPrev.t + dt )
    sim.setFlux()

    sim.setBoundaryCondition()
    

    return sim


def slope_limiter(sim, eigenval, name):
    """
    Calculates phi according to the set slope limiter.. for those that switch
    between values, this requires the eigenvalues of all characteristics in
    order to calculate the 'r' value...
    """
    
    N   = sim.grid.N
    phi = arrayList(N+1)

    def _calcr():
            
        r = arrayList(N+1)
        
        for i in range(3):
            
            j = 2
            while j < N - 2:
                
                denom = sim.q[i][j] - sim.q[i][j-1]
            
                if denom == 0:
                    r[i][j] = 0
                elif eigenval[i][j] > 0:
                    r[i][j] = (sim.q[i][j-1] - sim.q[i][j-2]) / denom
                    print sim.q[i][j-1] - sim.q[i][j-2]
                    print denom
                elif eigenval[i][j] < 0:
                    r[i][j] = (sim.q[i][j+1] - sim.q[i][j]) / denom
                else:
                    r[i][j] = 0.0
                    
                
                j = j + 1
                
            r[i][0] = 0.0 
            r[i][1] = 0.0
            r[i][-2] = 0.0
            r[i][-1] = 0.0
        r[ r == 0.0 ] = 0.0
        
        return r
    
    def _donor_cell():
        return phi # phi is fixed at zero in donor cell
    
    def _lax_wendroff():
        return phi + 1.0 # phi is fixed at one
    
    def _beam_warming():
        return _calcr()
    
    def _fromm():
        return 0.5*(1.0+_calcr())
    
    def _superbee():
        
        r = _calcr()
        
        for i in range(3):
        
            j = 0
            for j in range(N+1):
                b = np.min([1.0, 2.0*r[i][j]])
                c = np.min([2.0,r[i][j]])
            
                phi[i][j] = np.max(np.array([0.0,b,c]))
                
               
    
        return phi
    
    def _minmod():
        a = 1.0

        r = _calcr()
        
        for i in range(3):
        
            j = 0
            for j in range(N+1):

                if a*r[i][j] <= 0.0:
                    phi[i][j] = 0.0

                elif np.abs(a) > np.abs(r[i][j]):
                    phi[i][j] = r[i][j]
                    
                elif np.abs(a) < np.abs(r[i][j]):
                    phi[i][j] = a

        return phi
         
    phiDict = {'minmod': _minmod,
               'donor-cell': _donor_cell,
               'Lax-Wendroff': _lax_wendroff,
               'beam-warming': _beam_warming,
               'superbee': _superbee,
               'fromm': _fromm
               
              }   
        
        
    return phiDict[name]()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #




# Below runs the code and sets up a few things





init_cond = INIT_COND
##############################################################

# make grid with boundary
grid = simGrid(N)
grid.setLim(xmin,xmax)
grid.makeGrid()
grid.setBoundaryType(BOUNDARY_CONDITION)


# make simulation volume
sim = simState()
sim.initializeSimulation(grid) # make variables on grid   
set_initial_conditions(sim, init_cond,inputfile=None) # set initial conditions
   
   
tstep = 0    
sim.writeOut(outname="data/%.4d_"%(tstep)+init_cond,include_state=True) # write out initial conditions

# loop over the roe solver... output at specified times
while tstep < tstepmax:
    tstep   = tstep + 1

    simNext = roe_evolve(sim)  # evolve according to roe scheme
   
    if np.mod(tstep,1) ==0 or simNext.t == tmax:
        print "=======================",tstep,"=================="
        print simNext.t - sim.t
   
    if tstep in outarray or tstep == tstepmax-1:
        simNext.writeOut(outname="data/%.4d_"%(tstep)+init_cond,include_state=True)   # write out
    

    sim    = simNext


   
   
   
   
