import numpy as np

def define_mu():

    mu_e = {'H': 1.17,
            'He': 2.,
            'C' : 2.,
            'O' : 2.,
            'Si' : 2.}
    
    mu = {'H': 0.61,
          'He': 1.0/(1./2. + 1./4.),
          'C' : 1.0/(1./2. + 1./12.),
          'O' : 1.0/(1./2. + 1./16.),
          'Si': 1.0/(1./2. + 1./28.) }

    return mu, mu_e

def minimum_mass(fuel,T):

    mu, mu_e = define_mu()

    K = 1.E13
    R = 1.3806E-16 / 1.6733E-24
    G = 6.67E-8
    print fuel, mu[fuel], mu_e[fuel]


    return ( (15./np.pi)**(2./3.) * K * R * T\
           /(G**2 * mu[fuel]**2 * mu_e[fuel]**(2./3.) ) ) **(3./4.) / 1.99E33

fuel_list = ['H','H','He','C','O','Si']
core_temp = np.array([4.,15.,100.,600.,1000.,3000.]) * 1.0E6


for (fuel,T) in zip(fuel_list,core_temp):


    print fuel, "%5.4e %5.4e" % (T, minimum_mass(fuel,T) )
