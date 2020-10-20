import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid, 2mom
arr = [False, False, False, False]
#arr = [True, True, True, True]

def file_prefix(j):
    if len(str(j)) == 1:
       return '000'
    elif len(str(j)) == 2:
       return '00'
    elif len(str(j)) == 3:
       return '0'
    else:
       return 'Inappropriate length of input to file_prefix'

# How many vertical levels depends on which set we look at
if suffix1 == '_PL2':
   c = 120
elif suffix1 == '_PL':
   c = 18

WT_1mom = np.zeros((24,2,c))
WT_2mom = np.zeros((24,2,c))
WT_no2mom = np.zeros((12,2,c))
WT_novgrid = np.zeros((12,2,c))

# Which fraction of high cloud coverage are you requiring?
f = 0

if arr[0] == True:
   for i in np.arange(1,24):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run2/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'WT_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       WT_1mom[i-1,0] = flx.omega.isel(time=0).where(clch > f).mean(dim={'ncells'})
       WT_1mom[i-1,1] = flx.temp.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 1mom omega and T profiles...')
   np.save('../output/WT_1mom' + suffix2 + '.npy',WT_1mom)
else:
   WT_1mom = np.load('../output/WT_1mom' + suffix2 + '.npy')


if arr[1] == True:
   for i in np.arange(1,13):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_no2mom/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'WT_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       WT_no2mom[i-1,0] = flx.omega.isel(time=0).where(clch > f).mean(dim={'ncells'})
       WT_no2mom[i-1,1] = flx.temp.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving no2mom fluxes...')
   np.save('../output/WT_no2mom' + suffix2 + '.npy',WT_no2mom)
else:
   WT_no2mom = np.load('../output/WT_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   for i in np.arange(1,13):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_novgrid/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'WT_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       WT_novgrid[i-1,0] = flx.omega.isel(time=0).where(clch > f).mean(dim={'ncells'})
       WT_novgrid[i-1,1] = flx.temp.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving novgrid fluxes...')
   np.save('../output/WT_novgrid' + suffix2 + '.npy',WT_novgrid)
else:
   WT_novgrid = np.load('../output/WT_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   for i in np.arange(48,72):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'WT_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       WT_2mom[i-48,0] = flx.omega.isel(time=0).where(clch > f).mean(dim={'ncells'})
       WT_2mom[i-48,1] = flx.temp.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 2mom fluxes...')
   np.save('../output/WT_2mom' + suffix2 + '.npy',WT_2mom)
else:
   WT_2mom = np.load('../output/WT_2mom' + suffix2 + '.npy')


# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,6.5))
p = 0.1
#ax[0].plot(np.percentile(WT_1mom[:,0],1,axis=0),pl/100,color='red',label='ICON-1mom')
ax[0].plot((np.percentile(WT_2mom[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='green',label='ICON-2mom')
ax[0].plot((np.percentile(WT_no2mom[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='blue',label='ICON-no2mom')
ax[0].plot((np.percentile(WT_novgrid[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='gold',label='ICON-novgrid')
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'Domain-mean $\Delta\omega_{99.9}$ '
                  '\n'
                 'from 1-mom [hPa s$^{-1}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].legend(loc='upper right')

#ax1 = fig.add_axes([0.33,0.55,0.07,0.25])
#ax1.plot(np.percentile(WT_1mom[:,0],p,axis=0)*100,pl/100,color='red',label='_nolegend_')
#ax1.set_xlabel('$\omega_{99.9}$',fontsize=fs-3)
#ax1.set_yscale('log')
#ax1.set_ylim([50,800])
#ax1.set_yticks([800,500,300,100])
#ax1.set_yticklabels(['800','500','300','100'])
#ax1.invert_yaxis()

ax[0].set_ylim([50,800])
ax[0].set_yscale('log')
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)

#ax[1].plot(np.nanmean(WT_1mom[:,0],axis=0),pl/100,color='red',label='ICON-1mom')
ax[1].plot((np.nanmean(WT_2mom[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='green',label='ICON-2mom')
ax[1].plot((np.nanmean(WT_no2mom[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='blue',label='ICON-no2mom')
ax[1].plot((np.nanmean(WT_novgrid[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='gold',label='ICON-novgrid')
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
#ax[1].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[1].set_xlabel(r'Domain-mean daily-mean $\Delta\omega$'
                 '\n'
                 'from 1-mom [hPa s$^{-1}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)

ax2 = fig.add_axes([0.8,0.55,0.07,0.25])
ax2.plot(np.nanmean(WT_1mom[:,0],axis=0)*100,pl/100,color='red',label='_nolegend_')
ax2.set_xlabel('$<\omega>$',fontsize=fs-3)
ax2.set_yscale('log')
ax2.set_ylim([50,800])
ax2.set_yticks([800,500,300,100])
ax2.set_yticklabels(['800','500','300','100'])
ax2.invert_yaxis()

ax[1].set_ylim([50,800])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)

fig.savefig('../output/wT-profiles2.pdf')
plt.show()
