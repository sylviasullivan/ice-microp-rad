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
# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8

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

FLX_1mom = np.zeros((24,4,c))
FLX_2mom = np.zeros((24,4,c))
FLX_no2mom = np.zeros((24,4,c))
FLX_novgrid = np.zeros((24,4,c))

# Which fraction of high cloud coverage are you requiring?
f = 0

if arr[0] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run2/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       FLX_1mom[i-1,0] = flx.lwflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_1mom[i-1,1] = flx.lwflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_1mom[i-1,2] = flx.swflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_1mom[i-1,3] = flx.swflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 1mom fluxes...')
   np.save('../output/FLX_1mom' + suffix2 + '.npy',FLX_1mom)
else:
   FLX_1mom = np.load('../output/FLX_1mom' + suffix2 + '.npy')


if arr[1] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_no2mom/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       FLX_no2mom[i-1,0] = flx.lwflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_no2mom[i-1,1] = flx.lwflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_no2mom[i-1,2] = flx.swflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_no2mom[i-1,3] = flx.swflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving no2mom fluxes...')
   np.save('../output/FLX_no2mom' + suffix2 + '.npy',FLX_no2mom)
else:
   FLX_no2mom = np.load('../output/FLX_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_novgrid/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       FLX_novgrid[i-1,0] = flx.lwflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_novgrid[i-1,1] = flx.lwflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_novgrid[i-1,2] = flx.swflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_novgrid[i-1,3] = flx.swflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving novgrid fluxes...')
   np.save('../output/FLX_novgrid' + suffix2 + '.npy',FLX_novgrid)
else:
   FLX_novgrid = np.load('../output/FLX_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   for i in np.arange(49,73):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       FLX_2mom[i-49,0] = flx.lwflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_2mom[i-49,1] = flx.lwflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_2mom[i-49,2] = flx.swflxall.isel(time=0).where(clch > f).mean(dim={'ncells'})
       FLX_2mom[i-49,3] = flx.swflxclr.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 2mom fluxes...')
   np.save('../output/FLX_2mom' + suffix2 + '.npy',FLX_2mom)
else:
   FLX_2mom = np.load('../output/FLX_2mom' + suffix2 + '.npy')


# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,5.5))
ax[0].plot(np.nanmedian(FLX_1mom[:,0],axis=0),pl/100,color='red',label='ICON-1mom')
print(np.nanmedian(FLX_1mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(FLX_2mom[:,0],axis=0),pl/100,color='green',label='ICON-2mom')
print(np.nanmedian(FLX_2mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(FLX_no2mom[:,0],axis=0),pl/100,color='blue',label='ICON-no2mom')
print(np.nanmedian(FLX_no2mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(FLX_novgrid[:,0],axis=0),pl/100,color='gold',label='ICON-novgrid')
print(np.nanmedian(FLX_novgrid[:,0],axis=0).max())
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'LW cloud flux [W m$^{-2}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].set_ylim([50,800])
#ax[0].set_xlim([-1,1])
ax[0].set_yscale('log')
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)

ax[1].plot(np.nanmedian(FLX_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
ax[1].plot(np.nanmedian(FLX_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
ax[1].plot(np.nanmedian(FLX_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
ax[1].plot(np.nanmedian(FLX_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[1].set_xlabel(r'SW cloud flux [W m$^{-2}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)
ax[1].set_ylim([50,800])
#ax[1].set_xlim([-1,1])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)
ax[1].legend(fontsize=fs-3,loc='center right')
fig.suptitle('Daily-median domain-mean flux profiles',fontsize=fs)

#fig.savefig('../output/lw-sw-flux-profiles2.pdf')
plt.show()


fig = plt.figure(figsize=(8,4.5))
plt.plot(np.nanmedian(FLX_1mom[:,0]+FLX_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
plt.plot(np.nanmedian(FLX_2mom[:,0]+FLX_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
plt.plot(np.nanmedian(FLX_no2mom[:,0]+FLX_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
plt.plot(np.nanmedian(FLX_novgrid[:,0]+FLX_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative fluxes [W m$^{-2}$]',fontsize=fs)
plt.title('Daily-mean domain-mean flux profiles',fontsize=fs)
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.gca().tick_params(labelsize=fs)
plt.gca().legend(fontsize=fs-1)
#fig.savefig('../output/flux-profiles2.pdf')
plt.show()


fig = plt.figure(figsize=(8,4.5))
colors = plt.cm.jet(np.linspace(0,1,24))
for i in np.arange(24):
    plt.plot(FLX_1mom[i,0],pl/100,color=colors[i],label='ICON-1mom')
#plt.plot(np.nanmedian(FLX_2mom[:,0]+FLX_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
#plt.plot(np.nanmedian(FLX_no2mom[:,0]+FLX_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
#plt.plot(np.nanmedian(FLX_novgrid[:,0]+FLX_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative fluxes [W m$^{-2}$]',fontsize=fs)
plt.title('Hourly domain-mean flux profiles',fontsize=fs)
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.gca().tick_params(labelsize=fs)
plt.gca().legend(fontsize=fs-1)
#fig.savefig('../output/diurnal-flux-profiles.pdf')
plt.show()

