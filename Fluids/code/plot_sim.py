import numpy as np
import matplotlib.pyplot as plt
import os
plot_q = False
#time_dumps = [0,10000,20000,25000,50000,100000,200000,300000,400000,500000,600000,700000,800000,900000]
time_dumps = np.arange(0,50000,500)
time_dumps = np.arange(0,26,1)
name='sod1'
for i in time_dumps:
    infile = "data/%.4d_"%(i)+name+"_simstate.txt"
    #if 1:
    
    
    if os.path.isfile(infile):
        data = np.genfromtxt(infile,names=True)

        x   = data['x']
        rho = data["Density"]
        u   = data['Velocity']
        P   = data['Pressure']

        
        ax1 = plt.subplot(2,2,1)
        ax1.plot(x,rho,color='black')


        ax1.set_ylabel('Density')
        ax1.set_xlabel('x')
        
        ax2 = plt.subplot(2,2,2)
        ax2.plot(x,P,color='black')
        ax2.set_ylabel('Pressure')
        ax2.set_xlabel('x')
        
        ax3 = plt.subplot(2,2,3)
        ax3.plot(x,u,color='black')
        ax3.set_ylabel('Velocity')
        ax3.set_xlabel('x')
        
        
        e = P/((7.0/5.0 - 1.0)*rho)
        etot = e + 0.5 * u * u
                  
            # total enthalpy
        htot = etot + P/rho
            
        
        
        ax4 = plt.subplot(2,2,4)
        ax4.plot(x,htot,color='black')
        ax4.set_ylabel('H$_{tot}$')
        ax4.set_xlabel('x')
        
        plt.tight_layout()
        print "plotted", infile
        plt.savefig("plots/all_%.4d.png"%(i))
        plt.close()
        
        if plot_q:
        #################
            x   = data['x']
            q0 = data['q0']
            q1 = data['q1']
            q2 = data['q2']

            
            ax1 = plt.subplot(2,2,1)
            ax1.plot(x,q0,color='black')

            ax1.set_ylabel('q0')
            ax1.set_xlabel('x')
            
            ax2 = plt.subplot(2,2,2)
            ax2.plot(x,q1,color='black')
            ax2.set_ylabel('q1')
            ax2.set_xlabel('x')
            
            ax3 = plt.subplot(2,2,3)
            ax3.plot(x,q2,color='black')
            ax3.set_ylabel('q2')
            ax3.set_xlabel('x')
            
            plt.tight_layout()
            print "plotted", infile
            plt.savefig("plots/state/all_state_%.4d.png"%(i))
            plt.close()

