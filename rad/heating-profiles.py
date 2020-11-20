import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to recalculate heating rates for? Order is 1mom, no2mom, novgrid, 2mom
arr = [False, False, False, False, False]
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

# basedir is the base directory where the nc files are found.
# tag is how to label the output npy.
# f is the fraction of high cloud coveraged required.
# startindx and endinx are the files over which to iterate.
def meanProfile(basedir, tag, f, startindx, endindx, fileprefix):
    # How many vertical levels depends on which set we look at
    if suffix1 == '_PL2':
       c = 120
    elif suffix1 == '_PL':
       c = 18
    H = np.zeros((24,2,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        flx = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
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
    np.save('../output/H_' + tag + suffix2 + '.npy',H)
    return H


if arr[0] == True:
   H_1mom = meanProfile('/scratch/b/b380873/tropic_run2/','1mom',0,1,25,'RAD_3D')
else:
   H_1mom = np.load('../output/H_1mom' + suffix2 + '.npy')


if arr[1] == True:
   H_no2mom = meanProfile('/scratch/b/b380873/tropic_run5_no2mom/','no2mom',0,1,25,'FLX')
else:
   H_no2mom = np.load('../output/H_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   H_novgrid = meanProfile('/scratch/b/b380873/tropic_run5_novgrid/','novgrid',0,1,25)
   H_radnovgrid = meanProfile('/scratch/b/b380873/tropic_run7_radnovgrid/','radnovgrid',0,1,24,'FLX')
else:
   H_novgrid = np.load('../output/H_novgrid' + suffix2 + '.npy')
   H_radnovgrid = np.load('../output/H_radnovgrid' + suffix2 + '.npy')

if arr[3] == True:
   #H_2mom = meanProfile('/scratch/b/b380873/tropic_run5/','2mom',0,49,73)
   H_rad2mom = meanProfile('/scratch/b/b380873/tropic_run7_rad2mom/','rad2mom',0,1,24,'FLX')
else:
   H_2mom = np.load('../output/H_2mom' + suffix2 + '.npy')
   H_rad2mom = np.load('../output/H_rad2mom' + suffix2 + '.npy')


if arr[4] == True:
   #H_Atest = meanProfile('/scratch/b/b380873/tropic_run6/','Atest',0,1,24,'FLX')
   H_PDA = meanProfile('/scratch/b/b380873/tropic_run8_pda/','PDA',0,1,24,'FLX')
else:
   H_Atest = np.load('../output/H_Atest' + suffix2 + '.npy')
   H_PDA = np.load('../output/H_PDA' + suffix2 + '.npy')

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
basedir = '/scratch/b/b380873/tropic_run7_rad2mom/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0001' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,5.5))
ax[0].plot(np.nanmean(H_ERA5[:,0],axis=0),pl/100,color='black',linestyle='--',label='ERA5')
ax[0].plot(np.nanmedian(H_1mom[:,0],axis=0),pl/100,color='red',label='ICON-1mom')
ax[0].plot(np.nanmedian(H_2mom[:,0],axis=0),pl/100,color='green',label='ICON-2mom')
ax[0].plot(np.nanmedian(H_no2mom[:,0],axis=0),pl/100,color='blue',label='ICON-no2mom')
ax[0].plot(np.nanmedian(H_novgrid[:,0],axis=0),pl/100,color='gold',label='ICON-novgrid')
#ax[0].plot(np.nanmedian(H_Atest[:,0],axis=0),pl/100,color='black',label='ICON-Atest')
#ax[0].plot(np.nanmedian(H_rad2mom[:,0],axis=0),pl/100,color='purple',label='ICON-rad2mom')
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'LW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
#ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].set_xlim([-1,1])
ax[0].set_yscale('log')
ax[0].set_ylim([50,800])
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)
ax[0].legend(fontsize=fs-3,loc='upper right')

ax[1].plot(np.nanmean(H_ERA5[:,1],axis=0),pl/100,color='black',linestyle='--',label='ERA5')
ax[1].plot(np.nanmedian(H_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
ax[1].plot(np.nanmedian(H_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
ax[1].plot(np.nanmedian(H_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
ax[1].plot(np.nanmedian(H_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
#ax[1].plot(np.nanmedian(H_Atest[:,1],axis=0),pl/100,color='black',label='ICON-Atest')
#ax[1].plot(np.nanmedian(H_rad2mom[:,1],axis=0),pl/100,color='purple',label='ICON-rad2mom')
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[1].set_xlabel(r'SW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
#ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)
ax[1].set_ylim([50,800])
ax[1].set_xlim([-1,1])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)
#fig.suptitle('Daily-median domain-mean heating profiles',fontsize=fs)

fig.savefig('../output/lw-sw-heating-profiles_1mom+2mom+ERA5.pdf')
plt.show()


fig = plt.figure(figsize=(5.5,6.5))
plt.plot(np.nanmedian(H_ERA5[:,0]+H_ERA5[:,1],axis=0),pl/100,color='black',linestyle='--',label='ERA5')
plt.plot(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
plt.plot(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
plt.plot(np.nanmedian(H_no2mom[:,0]+H_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
plt.plot(np.nanmedian(H_novgrid[:,0]+H_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
#plt.plot(np.nanmedian(H_rad2mom[:,0]+H_rad2mom[:,1],axis=0),pl/100,color='purple',label='ICON-rad2mom')
#plt.plot(np.nanmedian(H_radnovgrid[:,0]+H_radnovgrid[:,1],axis=0),pl/100,color='purple',linestyle='--',label='ICON-radnovgrid')
#plt.plot(np.nanmedian(H_PDA[:,0]+H_PDA[:,1],axis=0),pl/100,color='pink',label='PDA')
jj = np.argwhere((np.abs(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0))<0.05))
kk = np.argwhere((np.abs(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0))<0.05))
print('Altitudes of lowest radiative heating (< 0.05 K per day): ')
print(pl.values[jj[:,0]])
print(pl.values[kk[:,0]])
jj = np.argwhere((np.abs(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0))>0.425))
kk = np.argwhere((np.abs(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0))>0.75))
print('Altitudes of largest radiative heating (> 0.425 / 0.75 K per day): ')
print(pl.values[jj[:,0]])
print(pl.values[kk[:,0]])
print('Largest radiative heatings at those altitudes: ')
print(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)[jj[:,0]])
print(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)[kk[:,0]])
jj = np.argwhere((np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)<-0.065))
kk = np.argwhere((np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)<-0.5))
print('Altitudes of largest radiative cooling (< -0.065 / -0.5 K per day): ')
print(pl.values[jj[:,0]])
print(pl.values[kk[:,0]])
print('Largest radiative coolings at those altitudes: ')
print(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)[jj[:,0]])
print(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)[kk[:,0]])
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative heating [K day$^{-1}$]',fontsize=fs)
#plt.title('Daily-mean domain-mean heating profiles',fontsize=fs)
plt.xlim([-0.75,1.5])
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.gca().tick_params(labelsize=fs)
plt.gca().legend(fontsize=fs-1)
fig.savefig('../output/heating-profiles_1mom+2mom+ERA5.pdf')
plt.show()


