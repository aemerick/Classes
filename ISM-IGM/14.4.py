import numpy as np

def _define_fit_dict():
    
    fitParamsDict = {"HII": [   12.25, 8.074E-6,     1.378, 5.087E2,  #C0-3
                             1.586E-2,   0.4723, 1.102E-5],     #C4-6
                     "CII": [   45.58, 6.089E-3,     1.128, 4.331E2,  #C0-3
                             4.845E-2,   0.8120, 1.333E-4]\
                    }                      
    
    
    return fitParamsDict


def alpha_gr(G,ne,T,species):
    CDict = _define_fit_dict()
    C = CDict[species]
    
    phi = phiCalc(G,T,ne)
    
    agr = 1.0E-14 *C[0]/\
          (1.0+C[1]*phi**(C[2])*\
          (1.0+C[3]*T**(C[4]) * phi**(-C[5]-C[6]*np.log(T))))
          
    return agr


def phiCalc(G,T,ne):
    return G*(T)**0.5/ne

# Results from first parts
PH = 7.01E-14
PC = 8.63E-14

# Givens from the problem
ne = 0.01  # cms
npr = 0.005 # cms
nH = 30.0  # cms
nC = 0.005 # cms
G  = 1.0 # 
T  = 100.0 # K

# Do this for Hydrogen
aHgr = alpha_gr(G,ne,T,species="HII")
aCgr = alpha_gr(G,ne,T,species="CII")

PHgr = aHgr*ne
PCgr = aCgr*ne

fgrH = (PHgr)/(PHgr + PH)
fgrC = (PCgr)/(PCgr + PC)

fH = (PH)/((PH) + (PHgr))
fC = (PC)/((PC) + (PCgr))

print "----- PART C AND D -----"
print "Alpha gr:"
print "HII : %.4E" %(aHgr)
print "CII : %.4E" %(aCgr)

print "P gr:"
print "HII : %.4E" %(PHgr)
print "CII : %.4E" %(PCgr)

print "----- PART E ------"
print "Fraction of recombinations by grains:"
print "HII : %.4f" %(fgrH)
print "CII : %.4f" %(fgrC)
print "Fraction of recombinations by electrons:"
print "HII : %.4f" %(fH)
print "CII : %.4f" %(fC)


