import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')

# import data
data = np.genfromtxt('hipparcos.txt', names=True)

# convert apparent V mag to absolute
selection = data['parallax'] > 0 # non-zero, non-negative parallax measurements

data = data[selection]
d    = 1000.0 / data['parallax'] # mas -> as => parsecs
Vmag = data['Vmag'] - 5.0*( np.log10(d) - 1.0)

# make a sub sample selection based upon quality of measurements
selection = ((data['err_parallax'] / data['parallax']) <  0.1)

# make sure the distances look OK
hist,bins = np.histogram(d[selection])
plt.plot((bins[:-1]+bins[1:])*0.5,hist, lw = 1.75)
plt.semilogy()
plt.ylim(1,1E5) ; plt.xlabel('distance (pc)'); plt.ylabel('counts')


# In[35]:

# select out the binaries
selection = selection * (data['multiplicity_flag'] == 0)

plt.scatter(data['BminusV'][selection], Vmag[selection], s=0.5, alpha=0.75)
plt.xlabel(r'B-V')
plt.ylabel(r'M$_V$')
plt.xlim(-0.5,2.0)
plt.ylim(15,-5) 

x1,y1 = 0.0, 2.0
x2,y2 = 1.5, 12.0
m1 = (y2-y1)/(x2-x1)
b1 = y2 - m1*x2

x1,y1 = 1.5,8.0
x2,y2 = 0.0, -1.0

m2 = (y2-y1)/(x2-x1)
b2 = y2 - m2*x2


plt.plot(data['BminusV'][selection],data['BminusV'][selection]*m1+b1,color = 'red')
plt.plot(data['BminusV'][selection],data['BminusV'][selection]*m2+b2,color='blue')





# In[36]:

# lets make this into a color histogram though
# I have a python script to do this already
# sorry... i realize this prevents you from running the script....
from plotting import plotTools as myplot

nbins = 75
xbins = np.linspace(-0.5,2.0,nbins)
ybins = np.linspace(-5,15,nbins)

BVmesh, Vmagmesh, H2d = myplot.my_histogram2d(data['BminusV'][selection],Vmag[selection],logscale=True, xbins=xbins, ybins=ybins)

mappable = plt.pcolormesh(BVmesh,Vmagmesh, H2d, cmap=plt.cm.gnuplot)
#plt.cm.gnuplot.set_bad('black',1.)

plt.ylim(15,-5)
plt.xlim(-0.5,2.0)
plt.xlabel(r'B-V')
plt.ylabel(r'M$_V$')


## Problem #2

# See paper for the derivations of A and B in terms of circular velocity and of proper motion in terms of A, d, and l

# We have $$\frac{v_{t}}{d}= Acos(2l) + B = \mu_{l}$$. 

# In[37]:

import emcee
import triangle


# In[38]:



def model(pars, x):
    A,B = pars
    return A*x + B

def ln_likelihood(pars, x, y, sigma_y):
    sum_term = 2.0*np.log(sigma_y) + ((y - model(pars, x))/sigma_y)**2
    return -0.5*(len(x)*np.log(2.0*np.pi) + np.sum(sum_term))

def ln_prior(pars, x, y, sigma_y):
    A,B = pars
    
    p = 0.
    
    A_range = (0.0,40.0)
    if A < A_range[0] or A > A_range[1]:
        return -np.inf
    else:
        p += -np.log(A_range[1]-A_range[0])
        
    B_range = (-40.,0.0)
    if B < B_range[0] or B > B_range[1]:
        return -np.inf
    else:
        p += -np.log(B_range[1] - B_range[0])
        
    return 0.

def ln_posterior(pars, *args):
    return ln_prior(pars, *args) + ln_likelihood(pars, *args)

def make_MS_cut(dummy, x, y, selection):
    
    # made two lines in the HR diagram and used them to select MS stars only
    # see the red and blue lines in the HR diagram above...
    if dummy:

        x1,y1 = 0.0, 2.0
        x2,y2 = 1.5, 12.0
        m1 = (y2-y1)/(x2-x1)
        b1 = y2 - m1*x2

        x1,y1 = 1.5,8.0
        x2,y2 = 0.0, -1.0

        m2 = (y2-y1)/(x2-x1)
        b2 = y2 - m2*x2
        
        selection = selection * (y > m2*x+b2)*(y < m1*x+b1)
        

    return selection
    


# For this analysis, I take $x = cos(2l)$, $y = \mu_l$, and $\sigma^2_{\rm{total}} = \sigma^2_{\mu_l} + \sigma^2_{v}$, where $\sigma_{\mu_l}$ is the uncertainty on the proper motion measurement, and $\sigma_{v}$ is the assumed isotropic velocity dispersion of the stars' random motion, fixed at 25 km/s. It is bad to assumed a fixed dispersion here, because the dispersion varies greatly (by a factor of a few) depending on stellar type / age, with a positive correlation in B-V. Assuming a fixed $\sigma_v$ will place less weight than there should be on stars with a low $\sigma_v$ (and low B-V), and place too much weight on stars with larger $\sigma_v$ (and high B-V). This will skew the resulting Oort constants, and likely result in larger estimated uncertainties on A and B.
# 
# In the below, I just show the analysis selecting the MS stars only. I make this selection using a separate function selecting only stars that lie within the red and blue lines shown in the scatter plot HR diagram above. I give the results of this analysis and the analsysis using the full sample at the end of the notebook.

# In[39]:

nwalkers = 40
nsteps = 2000

p0 = np.zeros((nwalkers,2))
p0[:,0] = np.random.uniform(4.,24,nwalkers)
p0[:,1] = np.random.uniform(-22,2,nwalkers)

# define x, y, and sigma_y from data
# lets do a cut though to remove some velocity outliers
v_t = 4.74*(d/1000.0)*data['mu_l']
selection = selection * (np.abs(v_t) < 125)

# make a selection for main sequence stars only
selection = make_MS_cut(True, data['BminusV'], Vmag, selection)

# cos(2l)
x = np.cos(2.0*data['l'][selection] * np.pi / 180.)

# mu_l
#y = 4.74 * data['mu_l'][selection] # km/s kpc^-1
y = data['mu_l'][selection] / 4.74
# uncertainty in the measured proper motion
#sigma_mu_l = 4.74 * data['err_mu_l'][selection]
sigma_mu_l = data['err_mu_l'][selection] / 4.74
# assuming isotropic velocity dispersion for each star, translated to dispersion in mu_l units
sigma_v = 25.0/(d[selection]/1000.0)

#sigma_v = (15.0 + data['BminusV'][selection]*15.0)/(d[selection]/1000.0)

# total uncertainty .... added in quadrature
sigma_y =  (sigma_v**2 + sigma_mu_l**2)**0.5 # in units of km/s kpc^-1


# In[40]:

# now plot the selected data to be fit to in the HR diagram 
plt.scatter(data['BminusV'][selection], Vmag[selection], s=0.5, alpha=0.75)
plt.xlabel(r'B-V')
plt.ylabel(r'M$_V$')
plt.xlim(-0.5,2.0)
plt.ylim(15,-5) 

x1,y1 = 0.0, 2.0
x2,y2 = 1.5, 12.0
m1 = (y2-y1)/(x2-x1)
b1 = y2 - m1*x2

x1,y1 = 1.5,8.0
x2,y2 = 0.0, -1.0

m2 = (y2-y1)/(x2-x1)
b2 = y2 - m2*x2


plt.plot(data['BminusV'][selection],data['BminusV'][selection]*m1+b1,color = 'red')
plt.plot(data['BminusV'][selection],data['BminusV'][selection]*m2+b2,color='blue')


# In[41]:


sampler = emcee.EnsembleSampler(nwalkers, dim=2, lnpostfn=ln_posterior,
                                args = (x,y, sigma_y)) 
pos, prob, state = sampler.run_mcmc(p0,nsteps//10) #???
sampler.reset()
pos, prob, state = sampler.run_mcmc(pos,nsteps)


# In[42]:

names = ['A','B']
for i in range(2):
    plt.figure()
    for walker in sampler.chain[...,i]:
        plt.plot(walker, drawstyle='steps', marker=None, linestyle='-', lw=1.5, alpha=0.5)
    plt.xlabel('Step Number')
    plt.ylabel(names[i], rotation='horizontal', labelpad=12)


# In[43]:

true_pars = [14.82,-12.37]

fig = triangle.corner(sampler.flatchain, truths=true_pars, labels=names, 
                      extents=[(0,40),(-40,0)])


# In[44]:

# plotting the data and fit line
plt.figure(figsize=(8,6))
plt.errorbar(x,y,sigma_y,marker='o', ls='none', ecolor='black',
             elinewidth=2, markersize=2,alpha=0.5)
x_l = np.linspace(np.min(x),np.max(x),100)

# draw 100 samples from the inferred
for i in np.random.randint(nsteps*nwalkers,size=100):
    plt.plot(x_l, model(sampler.flatchain[i], x_l), marker=None, 
             linestyle='-', color='#2166AC', alpha=0.2)
    
    
plt.plot(x_l, model(true_pars, x_l), marker=None, label='true',
         linestyle='-', color='green', linewidth=3.)
plt.xlabel("cos(2l)")
plt.ylabel(r"$\mu_l$ (km s$^{-1}$ kpc$^{-1}$)")
plt.legend(loc='best',fancybox=True)

# take the best fit parameters
histA, binA = np.histogram(sampler.flatchain[0])
histB, binB = np.histogram(sampler.flatchain[1])

A, err_A = np.median(sampler.flatchain[:,0]), np.std(sampler.flatchain[:,0])
B, err_B = np.median(sampler.flatchain[:,1]), np.std(sampler.flatchain[:,1])



print "%4.3f %4.3f %4.3f"%(A,np.average(sampler.flatchain[:,0]),err_A)
print "%4.3f %4.3f %4.3f"%(B,np.average(sampler.flatchain[:,1]),err_B)


# In the above computation of A and B using the Bayesian analysis approach, I only show the results using the MS stars. Using the full sample, however, I get: $$A = 15.23 \pm 3.57~ \rm{km}~\rm{s}^{-1}~\rm{kpc}^{-1}\\B = -16.07 \pm 2.61~ \rm{km}~ \rm{s}^{-1}~ \rm{kpc}^{-1}.$$ Selecting for only MS stars, I get $$A = 18.09 \pm 4.10~ \rm{km}~\rm{s}^{-1}~\rm{kpc}^{-1}\\B = -18.93 \pm 3.03~ \rm{km}~ \rm{s}^{-1}~ \rm{kpc}^{-1}.$$
# 
# The MS star only values are surprisingly inconsistent with the Oort constants identified from the Hipparcos team. I am not sure why this is the case. Yes the assumed isotropic velocity dispersion will play a significant role here, but I would not have expected it to have this dramatic of an effect. The constants obtained using the full sample aren't terribly great either, but are better. The A in the full sample is near Hipparcos' results, but my uncertainty is significantly larger. My B value in the full sample is inconsistent as well

# In[44]:



