import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def Tss(a,U,name):

    if name == 'graphite':
        return  22.3*(a/0.1)**(-1.0/40.0)*U**(1.0/6.0)
        
    elif name =='silicon':
        return  16.4*(a/0.1)**(-1.0/15.0)*U**(1.0/6.0)
        
a = np.linspace(-8,2,10000)
a = 10**(a)       
#U = np.linspace(0.1,1E7,10000)

for name in ['graphite','silicon']:
    fig = plt.figure(figsize = [6,6])
    ax = fig.add_subplot(111)
    plt.tight_layout() 
    for U in np.array([1,10,90,100,110,1000]):
        
        
        T = Tss(a,U,name)
        
        ax.plot(a,T,label='U = ' + str(U))
        
        
    ax.plot([np.min(a),np.max(a)],[50.0,50.0],label='T = 50 K',linestyle='-.',color='black',
             linewidth=1.5)
    ax.loglog()
    ax.set_ylabel('T (K)')
    ax.set_xlabel(r'a ($\mu$m)')
    ax.legend(loc='best',fancybox=True)
    plt.tight_layout() 
    fig.savefig(name+"_Tss.png")

    #plt.show()
    #plt.close()
    plt.close(fig)
