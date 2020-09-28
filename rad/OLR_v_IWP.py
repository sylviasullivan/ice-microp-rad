import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
#from scipy.stats import kde
import scipy.stats
from matplotlib import cm
import sys

basedir = '/work/bb1131/b380873/tropic_run2_output/'
fi1 = basedir + 'LWFLXALL_0085_remapdis_global0.025.nc'
lwfi = xr.open_dataset(fi1).lwflxall[0,0,::10,::10]
lwfi = -1*np.reshape(np.array(lwfi),(180*460,))
fi2 = basedir + 'CLCONV_2D_icon_tropic_0085_remapdis_global0.025.nc'
tqi = 1000*xr.open_dataset(fi2).tqi[0,::10,::10]
tqi = np.reshape(np.array(tqi),(180*460,))
del fi1
del fi2

def mean_confidence_interval(data,confidence=0.99):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.median(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h

tqibins = np.logspace(-3,2.8,50)
tqicenters = np.zeros((49,))
olr = np.zeros((49,2))
for i in np.arange(len(tqibins)-1):
    j = np.argwhere((tqi > tqibins[i]) & (tqi <= tqibins[i+1]))
    tqicenters[i] = (tqibins[i] + tqibins[i+1])/2
    olr[i,0] = np.nanmean(lwfi[j[:,0]])
    olr[i,1] = mean_confidence_interval(lwfi[j[:,0]])

fs = 13
fig = plt.figure(figsize=(5,5))
#plt.hist2d(tqi,lwfi,bins=[np.linspace(100,375,50),np.logspace(-3,2.5,50)],cmap=cm.viridis)
plt.plot(tqicenters,olr[:,0],color='darkblue',linewidth=1.25)
plt.fill_between(tqicenters,olr[:,0]-olr[:,1],olr[:,0]+olr[:,1],color='lightblue')

plt.gca().set_xscale('log')
plt.xlabel(r'Ice water path [g m$^{-2}$]',fontsize=fs)
plt.ylabel(r'Outgoing longwave radiation [W m$^{-2}$]',fontsize=fs)
#fig.savefig('iwp-olr.pdf',bbox_inches='tight')
plt.show()
