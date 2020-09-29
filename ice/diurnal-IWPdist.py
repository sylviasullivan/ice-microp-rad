import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

basedir = '/work/bb1131/b380873/tropic_run2_output/'
tqi_fi1 = xr.open_dataset(basedir + 'tqi_2017080800-2017080806.nc')
tqi_fi2 = xr.open_dataset(basedir + 'tqi_2017080806-2017080812.nc')
tqi_fi3 = xr.open_dataset(basedir + 'tqi_2017080812-2017080818.nc')
tqi_fi4 = xr.open_dataset(basedir + 'tqi_2017080718-2017080800.nc')


tqi_vals_unwrapped = np.zeros((4,12*1800*4600))

tqi_vals = tqi_fi1.tqi.values
tqi_vals_unwrapped[0] = np.reshape(tqi_vals,(12*1800*4600))
tqi_vals = tqi_fi2.tqi.values
tqi_vals_unwrapped[1] = np.reshape(tqi_vals,(12*1800*4600))
tqi_vals = tqi_fi3.tqi.values
tqi_vals_unwrapped[2] = np.reshape(tqi_vals,(12*1800*4600))
tqi_vals = tqi_fi4.tqi.values
tqi_vals_unwrapped[3] = np.reshape(tqi_vals,(12*1800*4600))


fs = 11
fig, ax = plt.subplots(2, 2, figsize=(9,9))
farbe = ['red','gold','green','blue']
titre = ['(a) 2017080800-2017080806 \n early morning',\
         '(b) 2017080806-2017080812 \n late morning',\
         '(c) 2017080812-2017080818 \n early afternoon',\
         '(d) 2017080718-2017080800 \n late afternoon']

c = 0
n = 40   # number of bins
h = np.zeros((n-1,4))
for i in np.arange(2):
    for j in np.arange(2):
        wgts = np.ones_like(tqi_vals_unwrapped[c])/len(tqi_vals_unwrapped[c])*100
        h1 = ax[i,j].hist(tqi_vals_unwrapped[c]*1000,weights=wgts,edgecolor='k',\
                 color=farbe[c],bins=np.logspace(-3.5,3.5,n))
        h[:,c] = h1[0]
        ax[i,j].set_title(titre[c])
        ax[i,j].set_ylim([0,4])
        if j == 0:
           ax[i,j].set_ylabel('Probability [%]',fontsize=fs)
        ax[i,j].set_xlim([0.0003,1500])
        if i == 1:
           ax[i,j].set_xlabel(r'IWP [kg m$^{-2}$]',fontsize=fs)
        ax[i,j].set_xscale('log')
        c += 1

print(h)
#fig.savefig('IWPdist-log-diurnal-abs.pdf',bbox_inches='tight')
plt.show()
