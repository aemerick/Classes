import numpy as np
import copy

FIXED_DT = [True, 0.1]
FLUX_LIM = 'Lax-Wendroff'
INIT_COND = 'sod1'
ADD_VISC = True
#    phiDict = {'minmod': _minmod,
 #              'donor-cell': _donor_cell,
  #             'Lax-Wendroff': _lax_wendroff,
   #            'beam-warming': _beam_warming,
    #           'superbee': _superbee,
     #          'fromm': _fromm

class simGrid:
    """
    Class defined grid. 

    """

    def __init__(self, N):
        self.N  = N
                
    def setLim(self,xmin,xmax):
        """
        """
        self.xmin = xmin
        self.xmax = xmax
        
    def setBoundaryType(self,boundary):
        """
            zero : boundaries are set to zero
            constant: boundary set to initial condition value
        """
        self.bc = boundary
        
    def makeGrid(self):
        """
        """
        self.grid = 1.0*np.linspace(self.xmin,self.xmax,N+1)
       # self.centers = 0.5*(self.grid[1:] - self.grid[:-1])
        self.dx   = self.grid[1] - self.grid[0]
        
    def center(self, index=None):
        """
        returns array of bin centers
        
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
       
    Specify name of initial conditions. To generalize this,
    need to take primitive variables as input paramaters... but
    too lazy for now. Current options are Sod 1 and Sod 2
    """ 
    N = state.grid.N
    xmin = state.grid.xmin
    xmax = state.grid.xmax
    c = state.grid.center()
    cgrid = 0.5*(xmax+xmin)
        
    rho = np.zeros(N)
    u = np.zeros(N)
    P = np.zeros(N)
        
        
    if not inputfile == None:    
        data = np.genfromtxt(inputfile,names=True)
        
        rho = data['Density']
        P   = data['Pressure']
        u   = data['Velocity']
        gamma = 7.0/5.0
        
    elif name == 'sod1':
        # below is as defined as 'test 1' on page 352 of comp. gas dynamics
        # nvm... using chapter 6 of the Max Plank notes.... Fig 6.5
        L = c < cgrid
        R = (c > cgrid) + (c==cgrid)
            
        sizeL, sizeR = np.size(L[L]), np.size(R[R])
         
             #    rho[L], rho[R] =   1.0, 0.125  # kg/m^3
       #       u[L],   u[R] =   0.0,   0.0  # m/s
       #       p[L],   p[R] = 1.0E5, 1.0E4  # N/m^2  
         #            gamma = 7.0 / 5.0     # not too sure about this......
        
        rho[L], rho[R] = 1.0E5 * np.ones(sizeL), 1.25E4 * np.ones(sizeR) # 
        u[L]  ,   u[R] = 0.0   * np.ones(sizeL), 0.0    * np.ones(sizeR)
        P[L]  ,   P[R] = 1.0   * np.ones(sizeL), 0.1    * np.ones(sizeR)
        gamma = 7.0 / 5.0
            
            
    elif name == 'linear':
        rho[0:0.25*N] = np.ones(0.25*N)*1.0E5
        rho[0.25*N:0.75*N] = np.linspace(1.0E5,1.0E4,0.5*N)
        rho[0.75*N:] = np.ones(0.25*N)*1.0E4
#        rho = np.linspace(1.0E5,1.0E4,N)
        u   = np.zeros(N)
        P[0:0.25*N] = np.ones(0.25*N)*1.0
        P[0.25*N:0.75*N] = np.linspace(1.0,0.1,0.5*N)
        P[0.75*N:] = np.ones(0.25*N)*0.1
        gamma = 7.0 / 5.0
    
    elif name == 'curve':
        l = np.pi*(2.0/(N))
        
        rho[0:0.25*N ] = np.ones(0.25*N)*1.0E4
        rho[0.75*N:]   = np.ones(0.25*N)*1.0E4
        rho[0.25*N:0.75*N] = (1.0E5-1.0E4)*np.sin(l*np.arange(0.5*N)) + 1.0E4
        
        P[0:0.25*N] = np.ones(0.25*N)*0.1
        P[0.75*N:]  = np.ones(0.25*N)*0.1
        P[0.25*N:0.75*N] = (1.0-0.1) * np.sin( l*np.arange(0.5*N)) + 0.1
        
        u[0:0.25*N] = np.ones(0.25*N)*0.01
        u[0.75*N:]  = np.ones(0.25*N)*0.01
        u[0.25*N:0.75*N] = (0.1-0.01) * np.sin( l*np.arange(0.5*N)) + 0.01
        
        l = np.pi*(1.0/(1.0*N))
        rho = (2.0E4-1.0E4)* np.sin( l*np.arange(N)) + 1.0E4
        P   = (1.0-0.5)    * np.sin( l*np.arange(N)) + 0.5
        u   = (0.1-0.01)   * np.sin( l*np.arange(N)) + 0.01
        u = np.ones(N)*0.00
        gamma = 0.75
    

    # set the simstate time and initial conditions....    
    state.setStateFromPrim(rho,u,P,gamma)
    state.time(0.0)
        
        
def arrayList(N, ndim=1):
    """
    Returns zeroed array with ndim*3 rows and N
    columns
    """
#    return [np.zeros(ndim*3) for i in range(N)]        


    return np.array([np.zeros(N) for i in range(ndim*3)])


class simState:
    """
        This is a object that completely defins the system
        at all times
    
    """

    def __init__(self):
        print "initalized"
        
    def initializeSimulation(self, grid):
        self.q = arrayList(grid.N)
        self.f = arrayList(grid.N)
        self.eigenVal = arrayList(grid.N)
        self.e1 = arrayList(grid.N)
        self.e2 = arrayList(grid.N)
        self.e3 = arrayList(grid.N)
        self.grid = grid
        
        self.setBoundaryCondition()

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
        xi = 3.0
    
        rho = self.getPrimitive('Density')
        u = self.getPrimitive('Velocity')
        P = self.getPrimitive('Pressure')
        
        visc = np.zeros(np.size(P))
        
        i = 1
        while i < np.size(P) - 1:
            if u[i+1] <= u[i-1]:
                visc[i] = 0.25*xi*xi*(u[i+1] - u[i-1])**2 * rho[i]
            else:
                visc[i]=0.0
            
            i = i + 1
        
        visc[0] = 0.0
        visc[-1] = 0.0
            
        
        
        self.setStateFromPrim(rho, u, P + visc, self.gamma)
    
    def setStateFromPrim(self, rho, u, P, gamma):
        """
        Given arrays of density, velocity, and pressure at all points,
        sets the state vector and flux
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

        # set the flux vector f_x
        #self.f[0] = rho*u
        #self.f[1] = rho * u * u
        #self.f[2] = rho*htot*u

        self.setFlux()

        try:
            t = self.t
        except:
            self.t = 0
            self.setBoundaryCondition()

    def setFlux(self):   
        self.f[0] = self.q[1]
  #      self.f[1] = self.q[1]**2 / self.q[0]
   #     self.f[2] = self.q[1] * self.getPrimitive('htot')
        self.f[1] = (self.gamma-1.0)*self.q[2] +\
                    0.5*(3.0-self.gamma)*(self.q[1]**2)/self.q[0]
                    
        self.f[2] = self.gamma*self.q[2]*self.q[1]/self.q[0] +\
                    0.5*(1.0-self.gamma)*(self.q[1])**3/(self.q[0])**2

    def getPrimitive(self, name, q = None):
        """
        Calculate a primitive variable from the state vector
        
        can provide a separate q vector to calculate primitives from
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



            #K = 
            #P = q[0]**(self.gamma) * K
            P[P < 0] = np.zeros(np.size( P[P<0] ) )
            return P
            
        def _etot():
            return q[2] / q[0]
        
        def _htot():
            return (q[2] + _pressure()) / q[0]    
        
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
    
        
       # np.savetxt(outname + '_simstate_%3.3f_.txt'%(self.t), data,
        np.savetxt(outname + '_simstate.txt', data, header=header, fmt='%1.4e')
        
        # np output
        
    def getInterfaceValue(self, name, direction):
        """
        direction is L or R
        """
        try:
            variable = self.getPrimitive(name)
        except:
            if 'q' in name:
                problem = True
                for i in range(3):
                    if str(i) in name:
                        index = i
                        problem = False

                variable = self.q[index]
                
                if problem:
                    print 'HOUSTON WE HAVE A PROBLEM'
            elif 'f' in name:
                problem = True
                for i in range(3):
                    if str(i) in name:
                        index = i
                        problem = False
                 
                if problem:
                    print 'HOUSTON WE HAVE A PROBLEM'
                variable = self.f[index]
        
        iVar = np.zeros( np.size(variable) + 1 ) # N + 1 interfaces
        
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
        # more bc's here
            
    
    #####
    def setBoundaryCondition(self):
        if self.grid.bc == 'constant':
            # conditions are fixed to their starting values at edges
            self.__qR__ = np.array([[self.q[0][-1]],[self.q[1][-1]],[self.q[2][-1]]])
            self.__qL__ = np.array([[self.q[0][0]] ,[self.q[1][0]] ,[self.q[2][0]]])
            
            self.__fR__ = np.array([[self.f[0][-1]],[self.f[1][-1]],[self.f[2][-1]]])
            self.__fL__ = np.array([[self.f[0][0]] ,[self.f[1][0]] ,[self.f[2][0]]])
            
            
            
       # print self.__qR__,self.__qL__
            
            
            

def roe_evolve(simPrev):
    C = 0.5 # fixed CFL condition

    N = simPrev.grid.N
    gamma = simPrev.gamma
    dx = simPrev.grid.dx
    # make n + 1 simulation state (simPrev + 1)
    #sim = simState()
    #sim.initializeSimulation(simPrev.grid)
    sim = copy.deepcopy(simPrev)
  #  print sim.grid is simPrev.grid
   # print simPrev.getPrimitive('Velocity')

    if ADD_VISC:
        simPrev.add_visc()
    # calculate dt with CFL condition
    max_vel = np.max(np.abs( simPrev.getPrimitive('Velocity') ) )
    

    


  #  print simPrev.getPrimitive('Density')
   # print simPrev.getPrimitive('Pressure')

    # below is all calculates for roe averages needed to calculate the eigs
    rhoL = simPrev.getInterfaceValue('Density','L')
    rhoR = simPrev.getInterfaceValue('Density','R')
    



   # print "time = " ,simPrev.time()
    denom = (rhoL)**0.5 + (rhoR)**0.5
    
    u_roe = ( (rhoL**0.5)*simPrev.getInterfaceValue('Velocity','L') +\
                        (rhoR**0.5)*simPrev.getInterfaceValue('Velocity','R') ) / denom

    #rho_roe = (( rhoL**0.5)*rhoL + (rhoR**0.5)*rhoR) / denom
    
    rho_roe = (rhoL * rhoR)**0.5  ### "standard practice' to use this ??
    
    htot_roe = ((rhoL**0.5)*simPrev.getInterfaceValue('htot','L') +\
                        (rhoR**0.5)*simPrev.getInterfaceValue('htot','R') ) / denom
   # print simPrev.getPrimitive('htot')
   # print htot_roe
  #  P_roe    = ((rhoL**0.5)*simPrev.getInterfaceValue('Pressure','L') +\
    #                    (rhoR**0.5)*simPrev.getInterfaceValue('Pressure','R') ) / denom
                        
   # Cs_roe   = (gamma * P_roe / rho_roe)**0.5
   
    Cs_roe = ((gamma-1.0) * (htot_roe - 0.5*u_roe*u_roe) )**0.5
    cfldt = simPrev.grid.dx / (np.max(u_roe) + np.max(Cs_roe))
    
    dt = FIXED_DT[1]

    
    if dt > cfldt:
        #dt = C * simPrev.grid.dx / (np.max(u_roe) + np.max(Cs_roe))
        dt = cfldt
  
  
   #d Cs_roe = (gamma*(gamma-1.0)/(gamma-2.0) *(htot_roe - 0.5*u_roe*u_roe))**0.5
   # print "--- after roe ---", dt
    # set up eigenvalues and eigenvectors
    eigenval = arrayList(N + 1)
    e1 = arrayList(N + 1) 
    e2 = arrayList(N + 1)
    e3 = arrayList(N + 1)
    
    # eigenvalues defined below
    eigenval[0] =  u_roe - Cs_roe
    eigenval[1] =  u_roe
    eigenval[2] =  u_roe + Cs_roe
    
    # eigenvectors defined below
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
    ekin = 0.5 * u_roe * u_roe
    
    dq_roe = arrayList(N+1)
    dq     = arrayList(N+1)
    
    dq[0]  = simPrev.getInterfaceValue('q0', 'R') -\
                                           simPrev.getInterfaceValue('q0', 'L')
    dq[1]  = simPrev.getInterfaceValue('q1', 'R') -\
                                           simPrev.getInterfaceValue('q1', 'L')
    dq[2]  = simPrev.getInterfaceValue('q2', 'R') -\
                                           simPrev.getInterfaceValue('q2', 'L')
  
#    dq[0] = ( (rhoL)**0.5 * simPrev.q[0] 
    
    A = (gamma - 1.0) / (2.0 * Cs_roe * Cs_roe)
    B = (dq[1] - u_roe * dq[0]) / (2.0 * Cs_roe)
    ksi = u_roe * dq[1] - dq[2]
    
    dq_roe[0] = A * (ekin * dq[0] - ksi) - B
    dq_roe[1] = A * ((htot_roe - 2.0*ekin)*dq[0] + ksi)
    dq_roe[2] = A * (ekin * dq[0] - ksi) + B
    
    # ---------------------------
    
    
    dq_roe_tilde = arrayList(N+1)
    for i in range(3):
        dq_roe_tilde[0] += dq_roe[i] * e1[i]
        dq_roe_tilde[1] += dq_roe[i] * e2[i]
        dq_roe_tilde[2] += dq_roe[i] * e3[i]
    
    # flux limiter
    theta = arrayList(N+1)
    theta = 1.0*np.sign(eigenval)
    
    eps = arrayList(N+1)
    eps = eigenval * dt / (1.0*dx)
    
    # calculate R
    #phi = arrayList(N+1) # upwind is just zeroes... leave as that
    
    # delta needs to be screwed around with
    # this is the Harten correction to the classic Roe scheme, which should be
    # able to handle both expansion and compression shocks. 
    # ---- phi is only not the eigenvalue if there is an expansion shock
    #delta = u_roe
    #for i in range(3):
    #    condition = eigenval[i] - delta
    #    ltarray =  condition < 0.0
    #    gtarray =  condition >= 0.0
    #    
    #    ltarray[delta == 0] = False
    #    gtarray[delta == 0] = True

    #    phi[i][ltarray] = (eigenval[i][ltarray]**2 + delta[ltarray]**2)/(2.0*delta[ltarray])
    #    phi[i][gtarray] = np.abs(eigenval[i][gtarray])
   
    # calculate r
   # for i in range(3):
   #     q_roe[i] = ((rhoL**0.5)*sim.getInterfaceValue('q%.1d'%d,'L') +\
   #                (rhoR**0.5)*sim.getInterfaceValue('1%.1d'%d,'R') ) / denom
    
    
    phi = slope_limiter(simPrev, u_roe, name=FLUX_LIM)
    

    
    #phi = 1.0+arrayList(N+1) # upwind is just zeroes... leave as that
    ######################
    favg = arrayList(N+1)
    favg[0] = (simPrev.getInterfaceValue('f0','R') +\
                                          simPrev.getInterfaceValue('f0', 'L'))
    favg[1] = (simPrev.getInterfaceValue('f1','R') +\
                                          simPrev.getInterfaceValue('f1', 'L'))
    favg[2] = (simPrev.getInterfaceValue('f2','R') +\
                                          simPrev.getInterfaceValue('f2', 'L'))
    
    favg = favg * 0.5
    
#    sigma = np.zeros(N+1)
    sigma = arrayList(N+1)
    for i in range(3):
        #sigma = sigma + eigenval[i] * dq_roe[i] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
        sigma[0] += eigenval[i] * dq_roe_tilde[0] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
        sigma[1] += eigenval[i] * dq_roe_tilde[1] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
        sigma[2] += eigenval[i] * dq_roe_tilde[2] * (theta[i] + phi[i] * (eps[i] - theta[i]) )
    
    #print np.max(favg), np.max(sigma)
    fnew = favg - 0.5*sigma

    
    # update sim

    
    for i in range(3):
        #sim.q[i] = simPrev.q[i] - (dt/dx)*(fnew[i][1:] - fnew[i][:-1])
        sim.q[i] = simPrev.q[i] + (dt/(1.0*dx))*(fnew[i][0:-1] - fnew[i][1:]) 
       # i = i + 1
    
    sim.setFlux()
    
    
    # set the time of the yet to be evolved simulation state
    sim.time( simPrev.t + dt )
    
    return sim
#    u_roe = 
 #   rho_roe = 
  #  htot_roe = 

def slope_limiter(sim, u_roe, name):
    """
    Calculates the slope limiter phi
    """
    
    N   = sim.grid.N
    phi = arrayList(N+1)

              

    

    
    
        
    
    def _calcr():
            
        r = arrayList(N+1)
        
        for i in range(3):
            
            j = 2
            for u in u_roe[2:-2]:
                
                denom = sim.q[i][j] - sim.q[i][j-1]
            
                if denom == 0:
                    r[i][j] = 0
                elif u > 0:
                    r[i][j] = (sim.q[i][j-1] - sim.q[i][j-2]) / denom
                elif u < 0:
                    r[i][j] = (sim.q[i][j+1] - sim.q[i][j]) / denom
                else:
                    r[i][j] = 0.0
                    
                
                j = j + 1
                
            r[i][0] = 0.0 
            r[i][1] = 0.0
            r[i][-2] = 0.0
            r[i][-1] = 0.0
    
    
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
        a = 1.0
        
        r = _calcr()
        
        for i in range(3):
        
            j = 0
            for j in range(N+1):
                b = np.min([a, 2.0*r[i][j]])
                c = np.min([2.0*a,r[i][j]])
            
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




# Below are "input" that can be set by external file later
N    =  255# NUMBER OF BINS, NOT NUMBER OF POINTS!!!!!

xmin, xmax = 0.0, 100.0



boundary = 'constant'
init_cond = INIT_COND
##############################################################

# make grid with boundary
grid = simGrid(N)
grid.setLim(xmin,xmax)
grid.makeGrid()
grid.setBoundaryType("constant")
  
#inputfile = "./data/130000_sod1_simstate.txt"
# make simulation volume
sim = simState()
sim.initializeSimulation(grid) # make variables on grid   
set_initial_conditions(sim, init_cond,inputfile=None) # set initial conditions
   
   
tstep = 0    
sim.writeOut(outname="data/%.4d_"%(tstep)+init_cond,include_state=True) # write out initial conditions


  #  while tstep < tstepmax:


#outarray = [0,10000,20000,25000,50000,100000,200000,300000,400000,500000,600000,700000,800000,900000]
#outarray = [100000,110000,120000,130000,140000]
outarray = np.arange(0,50000,100)
#tstep = 130000
#tstepmax = 1000001
tstepmax = 50005
while tstep < tstepmax:
    
    tstep   = tstep + 1

    simNext = roe_evolve(sim)  # evolve according to roe scheme
   
    if np.mod(tstep,100) ==0:
        print "=======================",tstep,"=================="
        print simNext.t - sim.t
   
    if tstep in outarray or tstep == tstepmax-1:
        simNext.writeOut(outname="data/%.4d_"%(tstep)+init_cond,include_state=True)   # write out
    

    sim    = simNext
#        simOld = sim # copy over sim...   
   
   
   
   
