import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'

# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2
#H = np.load('../output/H_1mom' + suffix2 + '.npy')
H = np.load('../output/H_2mom' + suffix2 + '.npy')
#H = np.load('../output/H_no2mom' + suffix2 + '.npy')
#H = np.load('../output/H_novgrid' + suffix2 + '.npy')
colors = plt.cm.jet(np.linspace(0,1,24))

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,5.5))
for i in np.arange(24):
    ax[0].plot(H[i,0],pl/100,color=colors[i],label='ICON-1mom')
    ax[1].plot(H[i,1],pl/100,color=colors[i],label='ICON-1mom')

ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'LW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].set_ylim([50,800])
ax[0].set_xlim([-2,2])
ax[0].set_yscale('log')
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)

ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[1].set_xlabel(r'SW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)
ax[1].set_ylim([50,800])
ax[1].set_xlim([-2,2])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)
#plt.colorbar()
#ax[1].legend(fontsize=fs-3,loc='center right')
fig.suptitle('Hourly, domain-mean heating profiles',fontsize=fs)

#fig.savefig('../output/diurnal-heating-rates-2mom.pdf')
plt.show()


