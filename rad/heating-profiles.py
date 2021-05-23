import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import sys
sys.path.append('/work/bb1018/b380873/tropic_vis/')
from plotting_utilities import sim_colors, sim_ls, file_prefix, sexy_axes2

# Which set of pressure levels to look at?
suffix = '_PL2' # '_PL'

# Which simulation acronyms to calculate heating rates for?
acronym = []
others = ['0V1M0A0R', '1V1M0A0R', '0V2M0A0R', '1V2M0A0R', '0V2M0A1R', '1V2M0A1R',
          '0V2M1A0R', '1V2M1A0R', '1V2M1A1R'] # '0V2M1A1R',

# Heat capacity [J kg-1 K-1]
cv = 0.718*10**(3)   # constant volume
cp = 1.08*10**(3)    # constant pressure

# Gravity [m s-2]
g = 9.8

# basedir is the base directory where the nc files are found.
# acronym is how to label the output npy.
# f is the fraction of high cloud coveraged required.
# startindx and endinx are the files over which to iterate.
def meanProfile(basedir, acronym, f, startindx, endindx, fileprefix):
    # How many vertical levels depends on which set we look at
    if suffix == '_PL2':
       c = 120
    elif suffix == '_PL':
       c = 18
    H = np.zeros((24,2,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        flx = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix + '.nc')
        clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
        # The factor of -1 is because pressure decreases upward
        temp = -1.*(flx.lwflxall.isel(time=0)-flx.lwflxclr.isel(time=0)).where(clch > f)
        temp = temp.assign_coords(plev_2=temp.plev_2.values)
        temp = temp.differentiate('plev_2')
        # Factor of 86400 to convert K s-1 to K day-1.
        H[i-startindx,0] = g/cp*86400*temp.mean(dim={'ncells'})
        temp = -1.*(flx.swflxall.isel(time=0)-flx.swflxclr.isel(time=0)).where(clch > f)
        temp = temp.assign_coords(plev_2=temp.plev_2.values)
        temp = temp.differentiate('plev_2')
        H[i-startindx,1] = g/cp*86400*temp.mean(dim={'ncells'})
    print('Saving mean heating rates from ' + basedir + '...')
    np.save('../output/H_' + acronym + suffix + '.npy',H)
    return H

for a in acronym:
   H = meanProfile('/scratch/b/b380873/' + a + '/', a, 0, 1, 25, 'FLX')

ll = len(others + acronym)
H_all = np.zeros((ll, 24, 2, 120))
for indx, o in enumerate(others + acronym):
    H_all[indx] = cv/cp * np.load('../output/H_' + o + suffix + '.npy')


# Load the ERA5 heating profiles. Subtract all-sky and clear-sky values.
# Multiply by 86400 to convert K s-1 to K day-1.
basedir2 = '/work/bb1018/b380873/tropic_vis/obs/ERA5/'
vals_ERA5 = xr.open_dataset(basedir2 + 'ERA5_ddtmean_20170807-20170808_55e115e5s40n_PL.nc')
H_ERA5 = np.zeros((24,2,120))
vals = vals_ERA5.mttlwr.mean(dim={'lat','lon'}) - vals_ERA5.mttlwrcs.mean(dim={'lat','lon'})
H_ERA5[:,0] = vals.values[6:30]*86400
vals = vals_ERA5.mttswr.mean(dim={'lat','lon'}) - vals_ERA5.mttswrcs.mean(dim={'lat','lon'})
H_ERA5[:,1] = vals.values[6:30]*86400

# Retrieve the pressure levels
pl = np.loadtxt('../remapping/PMEAN_48-72.txt')
farbe = sim_colors()
stil = sim_ls()

fs = 13
fig, ax = plt.subplots( nrows=1, ncols=2, figsize=(9,5.5) )
ax[0].plot( np.nanmean(H_ERA5[:,0],axis=0), pl/100, color='black', ls='-', label='ERA5' )
ax[0].plot([0,0], [50,800], lw=0.75, linestyle='--', color='k')
for h, o in zip(H_all, others + acronym):
    ax[0].plot( np.nanmedian(h[:,0],axis=0), pl/100, color=farbe[o[1:]], ls=stil[o[:2]], label=o )
ax[0].set_ylabel('Pressure [hPa]', fontsize=fs)
ax[0].set_xlabel(r'LW cloud radiative heating [K day$^{-1}$]', fontsize=fs)
ax[0].text(0.05, 0.92, '(a)', weight='bold', fontsize=fs+4, transform=ax[0].transAxes)
ax[0].set_xlim([-1,1])
ax[0].legend(fontsize=fs-3,loc='upper right')
sexy_axes2(ax[0], fs, True)

ax[1].plot(np.nanmean(H_ERA5[:,1],axis=0),pl/100,color='black',linestyle='-',label='ERA5')
for h, o in zip(H_all, others + acronym):
    ax[1].plot(np.nanmedian(h[:,1],axis=0), pl/100, color=farbe[o[1:]], ls=stil[o[:2]], label=o)
ax[1].plot([0,0], [50,800], lw=0.75, linestyle='--', color='k')
ax[1].set_xlabel(r'SW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[1].text(0.05, 0.92, '(b)', weight='bold', fontsize=fs+4, transform=ax[1].transAxes)
ax[1].set_xlim([-1,1])
sexy_axes2(ax[1], fs, True)
fig.suptitle('Daily-median domain-mean heating profiles',fontsize=fs)
#fig.savefig('../output/lw-sw-heating-profiles_all.pdf')
plt.show()


fig = plt.figure(figsize=(5.5,6.5))
plt.plot(np.nanmedian(H_ERA5[:,0]+H_ERA5[:,1],axis=0), pl/100, color='black', ls='-', label='ERA5')
for h, o in zip(H_all, others + acronym):
    plt.plot(np.nanmedian(h[:,0]+h[:,1],axis=0), pl/100, color=farbe[o[1:]], ls=stil[o[:2]], label=o)
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative heating [K day$^{-1}$]',fontsize=fs)
#plt.title('Daily-mean domain-mean heating profiles',fontsize=fs)
plt.xlim([-0.75,1.25])
sexy_axes2(plt.gca(), fs, True)
plt.gca().legend(fontsize=fs-1)
#fig.savefig('../output/heating-profiles_all.pdf')

# Heating difference between 0V2M1A0R and 1V2M1A0R simulations
print(np.nanmedian(H_all[6,:,0]+H_all[6,:,1],axis=0) - np.nanmedian(H_all[7,:,0]+H_all[7,:,1],axis=0))
#print(np.nanmedian(H_PDA[:,0]+H_PDA[:,1],axis=0)/np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0))
#jj = np.argwhere((np.abs(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0))<0.05))
#kk = np.argwhere((np.abs(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0))<0.05))
#print('Altitudes of lowest radiative heating (< 0.05 K per day): ')
#print(pl.values[jj[:,0]])
#print(pl.values[kk[:,0]])
#jj = np.argwhere((np.abs(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0))>0.425))
#kk = np.argwhere((np.abs(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0))>0.75))
#ll = np.argmax(np.nanmedian(H_rad2mom[:,0]+H_rad2mom[:,1],axis=0))
#print('Altitudes of largest radiative heating (> 0.425 / 0.75 K per day): ')
#print(pl[jj[:,0]])
#print(pl[kk[:,0]])
#print(pl[ll[:,0]])
#print('Largest radiative heatings at those altitudes: ')
#print(cv/cp*np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)[jj[:,0]])
#print(cv/cp*np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)[kk[:,0]])
#print(cv/cp*np.nanmax(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)))
#print(cv/cp*np.nanmax(np.nanmedian(H_rad2mom[:,0]+H_rad2mom[:,1],axis=0)))
#3print(cv/cp*np.nanmax(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)))
#jj = np.argwhere((np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)<-0.065))
#kk = np.argwhere((np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)<-0.5))
##print('Altitudes of largest radiative cooling (< -0.065 / -0.5 K per day): ')
##print(pl[jj[:,0]])
##print(pl[kk[:,0]])
#print('Largest radiative coolings at those altitudes: ')
#print(cv/cp*np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)[jj[:,0]])
#print(cv/cp*np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)[kk[:,0]])

plt.show()
