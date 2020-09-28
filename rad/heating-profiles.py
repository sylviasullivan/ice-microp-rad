import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid, 2mom
arr = [False, False, False, False]
#arr = [False, True, True, True]
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

H_1mom = np.zeros((24,2,c))
H_2mom = np.zeros((24,2,c))
H_no2mom = np.zeros((24,2,c))
H_novgrid = np.zeros((24,2,c))

# Which fraction of high cloud coverage are you requiring?
f = 0.5

if arr[0] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run2/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       temp = (flx.lwflxall.isel(time=0)-flx.lwflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       # Factor of 86400 to convert K s-1 to K day-1.
       H_1mom[i-1,0] = g/cp*86400*temp.mean(dim={'ncells'})
       temp = (flx.swflxall.isel(time=0)-flx.swflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       H_1mom[i-1,1] = g/cp*86400*temp.mean(dim={'ncells'})
   print('Saving 1mom heating rates...')
   np.save('../output/H_1mom' + suffix2 + '.npy',H_1mom)
else:
   H_1mom = np.load('../output/H_1mom' + suffix2 + '.npy')


if arr[1] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_no2mom/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       temp = (flx.lwflxall.isel(time=0)-flx.lwflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       # Factor of 86400 to convert K s-1 to K day-1.
       H_no2mom[i-1,0] = g/cp*86400*temp.mean(dim={'ncells'})
       temp = (flx.swflxall.isel(time=0)-flx.swflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       H_no2mom[i-1,1] = g/cp*86400*temp.mean(dim={'ncells'})
   print('Saving no2mom heating rates...')
   np.save('../output/H_no2mom' + suffix2 + '.npy',H_no2mom)
else:
   H_no2mom = np.load('../output/H_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   for i in np.arange(1,25):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_novgrid/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       print(clch)
       temp = (flx.lwflxall.isel(time=0)-flx.lwflxclr.isel(time=0)).where(clch > f)
       print(temp)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       # Factor of 86400 to convert K s-1 to K day-1.
       H_novgrid[i-1,0] = g/cp*86400*temp.mean(dim={'ncells'})
       temp = (flx.swflxall.isel(time=0)-flx.swflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       H_novgrid[i-1,1] = g/cp*86400*temp.mean(dim={'ncells'})
   print('Saving novgrid heating rates...')
   np.save('../output/H_novgrid' + suffix2 + '.npy',H_novgrid)
else:
   H_novgrid = np.load('../output/H_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   for i in np.arange(49,73):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'FLX_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       temp = (flx.lwflxall.isel(time=0)-flx.lwflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       # Factor of 86400 to convert K s-1 to K day-1.
       H_2mom[i-49,0] = g/cp*86400*temp.mean(dim={'ncells'})
       temp = (flx.swflxall.isel(time=0)-flx.swflxclr.isel(time=0)).where(clch > f)
       temp = temp.assign_coords(plev_2=temp.plev_2.values)
       temp = temp.differentiate('plev_2')
       H_2mom[i-49,1] = g/cp*86400*temp.mean(dim={'ncells'})
   print('Saving 2mom heating rates...')
   np.save('../output/H_2mom' + suffix2 + '.npy',H_2mom)
else:
   H_2mom = np.load('../output/H_2mom' + suffix2 + '.npy')


# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,5.5))
ax[0].plot(np.nanmedian(H_1mom[:,0],axis=0),pl/100,color='red',label='ICON-1mom')
print(np.nanmedian(H_1mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(H_2mom[:,0],axis=0),pl/100,color='green',label='ICON-2mom')
print(np.nanmedian(H_2mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(H_no2mom[:,0],axis=0),pl/100,color='blue',label='ICON-no2mom')
print(np.nanmedian(H_no2mom[:,0],axis=0).max())
ax[0].plot(np.nanmedian(H_novgrid[:,0],axis=0),pl/100,color='gold',label='ICON-novgrid')
print(np.nanmedian(H_novgrid[:,0],axis=0).max())
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'LW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].set_xlim([-1,1])
ax[0].set_yscale('log')
ax[0].set_ylim([50,800])
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)

ax[1].plot(np.nanmedian(H_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
ax[1].plot(np.nanmedian(H_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
ax[1].plot(np.nanmedian(H_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
ax[1].plot(np.nanmedian(H_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[1].set_xlabel(r'SW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)
ax[1].set_ylim([50,800])
ax[1].set_xlim([-1,1])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)
ax[1].legend(fontsize=fs-3,loc='center left')
fig.suptitle('Daily-median domain-mean heating profiles',fontsize=fs)

fig.savefig('../output/lw-sw-heating-profiles3.pdf')
plt.show()

sys.exit()

fig = plt.figure(figsize=(8,4.5))
plt.plot(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
plt.plot(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
plt.plot(np.nanmedian(H_no2mom[:,0]+H_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
plt.plot(np.nanmedian(H_novgrid[:,0]+H_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative heating [K day$^{-1}$]',fontsize=fs)
plt.title('Daily-mean domain-mean heating profiles',fontsize=fs)
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.gca().tick_params(labelsize=fs)
plt.gca().legend(fontsize=fs-1)
#fig.savefig('../output/heating-profiles2.pdf')
plt.show()


