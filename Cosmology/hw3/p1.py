import matplotlib.pyplot as plt
import numpy as np
plt.rc('font',size=16)

print "Omega baryon *h*h", 0.05*0.7*0.7
print "cdm", 0.25*0.7*0.7, 0.5*0.7*0.7
print "lam", 0.7*0.7*0.7

#use_real = True
use_rad = False
# model values of omega's
omega_baryon = [0.05, 0.05, 0.05]
omega_cdm    = [0.25, 0.25, 0.50]
omega_lambda = [0.7, 0.0, 0.7]
omega_rad    = 0.000

i = 0
h = 0.7
for i in range(len(omega_baryon)):
    print '-------------- model', i + 1, ' ----------------------'
    print 'baryon', omega_baryon[i] *h *h
    print 'cdm'   , omega_cdm[i] * h*h
    print 'lambda', omega_lambda[i]
    print 'curvature', 1.0 - omega_baryon[i] - omega_cdm[i] - omega_lambda[i] - omega_rad/(h*h)


lw = 1.5

planck1 = np.genfromtxt('planck_data.txt',names=True,delimiter=",")
planck2 = np.genfromtxt('planck_data_2.txt',names=True,delimiter=",")

#if use_real:
#    model1 = np.genfromtxt('model1_real.dat', names=True)
#    model2 = np.genfromtxt('model2_real.dat', names=True)
#    model3 = np.genfromtxt('model3_real.dat', names=True)
#    radstr = ''

if not use_rad:
    model1 = np.genfromtxt('model1.dat', names=True)
    model2 = np.genfromtxt('model2.dat', names=True)
    model3 = np.genfromtxt('model3.dat', names=True)
    radstr = ''
else:
    model1 = np.genfromtxt('model1_r.dat', names=True)
    model2 = np.genfromtxt('model2_r.dat', names=True)
    model3 = np.genfromtxt('model3_r.dat', names=True)
    radstr = 'rad'

models = {'model1': model1, 'model2': model2, 'model3': model3}
color  = {'model1': 'red', 'model2': 'blue', 'model3':'green'}


#plt.scatter(planck1['ELL'], planck1['D_ELL'], label='Planck Data', color='black')
plt.errorbar(planck1['ELL'], planck1['D_ELL'], yerr=[planck1['ERRUP'],planck1['ERRDOWN']],
               color='black', fmt='--o')
               
#plt.scatter(planck2['ELL'], planck2['D_ELL'], color='black')
plt.errorbar(planck2['ELL'], planck2['D_ELL'], yerr=planck2['ERR'], color='black', fmt='o')


i = 0
for model in models:
    plt.plot(models[model]['l'], models[model]['yval'],
              label = model, color=color[model] , lw=lw)
    i = i + 1


plt.xlim(0,1500)
plt.ylim(0,6500)
plt.xlabel(r'l $\sim$ 180$\deg$/$\theta$')
plt.ylabel(r'l(l+1)C$_{l}$/2$\pi$ [$\mu$K$^{2}$]')
plt.legend(loc='best',fancybox=True)

plt.savefig('p1_' + radstr + '.png')

plt.close()
