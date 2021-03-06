import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib import rc

# constants from WMAP and as given
#fH = 0.75
fH = 0.93
no = 2.541E-7
T4 = 2.0

#fHe/fH = 0.07

def heating(z,species):

    h = 6.6260755E-27 # planck constant in cgs

    if species=='HI':
        # number density of neutral hydrogen
        xHI = calcxHI(z,T4)
        n = fH * xHI * no*(1.0+z)**3
        
        Z = 1.0
        I = 13.59 * 1.602E-12 # ionization energy in erg
        vi = I / h
        
    if species == 'HeII':
        xHeII = calcxHeII(z,T4)
        # right part gives fraction of total He
        # multiplaction by HeII frac gives number density of HeII
        n   = xHeII * (1.0-fH)*no*(1.0+z)**3
        
        Z = 2.0
        I = 54.416 * 1.602E-12 # ionization energy in erg
        vi = I / h
        
    I = 13.59 * 1.602E-12 # i think the i in the equation below is that of H
    sigma_o = 6.304E-18 * (Z**(-2.0))
    
    # de = crap dv... function to integrate
    de=  lambda v:  4.0*np.pi*v**(-1.2)/(h*v) * h*(v-vi)*\
          sigma_o*( (Z*Z*I/(h*v))**4.0) *\
          np.exp(4.0-4.0*np.arctan((h*v/(Z*Z*I)-1)**0.5)/(h*v/(Z*Z*I) -1)**0.5)\
          /(1.0 - np.exp(-2.0*np.pi / (h*v/(Z*Z*I) -1)**0.5 )) 

   # print "here"
    e, error = integrate.quad(de, vi*1.00000000000000012, np.inf)
    
    de2=  lambda v:  4.0*np.pi*v**(-1.2)/(h*v) * h*(v-vi)*sigma_o
    e2, error = integrate.quad(de, vi,vi*1.00000000000000012)
    e = e +e2 # for the v = v_i bound
   # print "made it past integrating"
  #  e2, error2 = integrate.fixed_quad(de, vi*1.0001,np.inf)
  #  print e1
   # print n
    H = e*n
    print "eeee"
    print e,error
    
    return H
    

def aB(T4,species):
    # taken from 14.5 in Draine
    if species == 'HI':
        Z = 1
    elif species == 'HeII':
        Z = 2 
    
    ab = 4.13E-13 * Z * (T4/(Z*Z))**(-0.7131-0.0115*np.log(T4/(Z*Z)))
    return ab
    
def Gamma(z,species):
    # as given in the problem
    gamma = 6.7E-13 * (1.0+z)**(0.73)
    gamma = gamma * np.exp(- (z-2.30)**2 / 1.90 )
    
    if species == 'HeII':
        gamma = gamma * 0.01
    
    return gamma
    
    
def calcxHI(z,T4):
    species = 'HI'

    a = fH*no*(1.0+z)**(3.0)*aB(T4,species)/Gamma(z,species)
    
    return 1.0 / (1.0 + a**(-1.0))
    
def calcxHeII(z,T4):
    species = 'HeII'
    
  #  a = Gamma(z,species) / (aB(T4,species) * fH * no *(1.0+z)**3.0 *(1.0-calcxHI(z,T4)))
    
    A = fH * no * (1.0+z)**3.0 * (1.0 - calcxHI(z,T4)) * aB(T4,species) / Gamma(z,species)
    
    a = (1.0 + A**(-1.0))**(-1.0)
    
    return a
#    a = (fH)*no*(1.0 + z)**(3.0)*aB(T4,species)/Gamma(z,species)
    
  #  return 1.0 / (1.0 + a**(-1.0))

    
z = np.linspace(8.0,0.0,1000.0)

xHI = calcxHI(z,T4)
xHeII = calcxHeII(z,T4)


# find where xHI and xHeII are at 0.01 %
diff_HI = np.abs(xHI - 0.01)
diff_HeII = np.abs(xHeII - 0.01)

min_HI = diff_HI.argsort()[0]
min_HeII = diff_HeII.argsort()[0]

z_HI = z[min_HI]
z_HeII = z[min_HeII]


# plot xhi
plt.plot(z,xHI,label =r'x$_{HI}$',linewidth=1.5,color="black")
plt.xlabel("z")
plt.ylabel("Fraction")
plt.xlim([np.max(z),np.min(z)])
plt.semilogy()
xlim = plt.xlim()
ylim = plt.ylim()

plt.plot(z,xHeII,label=r'x$_{HeII}$',linewidth=1.5,color="black",linestyle='--')

# plot 1 % line
plt.plot([xlim[0],xlim[1]],[0.01,0.01], label='1% line')

plt.plot([z_HI,z_HI],ylim,label='z = %0.3f' % z_HI,linewidth=1.5,linestyle=':',color='red')
plt.plot([z_HeII, z_HeII],ylim,label='z = %0.3f' % z_HeII,linewidth=1.5,linestyle=':',color='red')

print ylim
print z_HI
print z_HeII


plt.legend(loc="best",fancybox=True)
plt.savefig("./xHI_xHe.png")
plt.close()

#################################################################################

H_HI = heating(z,'HI')
H_HeII = heating(z,'HeII')

z_trans = z[np.where(np.abs(H_HeII - H_HI)<1E-5)][0]

plt.plot(z,H_HI,label="HI Heating", linewidth=1.5, linestyle='-')
plt.plot(z,H_HeII,label="HeII Heating", linewidth=1.5,linestyle='--')
plt.xlabel("z")
plt.ylabel(r"log(H) K s$^{-1}$ cm$^{-3}$")

plt.xlim(xlim)
plt.semilogy()
ylim = plt.ylim()
plt.plot([z_trans,z_trans],ylim,label="z = %0.4f" % z_trans,linewidth=1.5,linestyle=":",color="red")




plt.legend(loc="best",fancybox=True)
plt.savefig("./HI_HeII_heating.png")
plt.close()

