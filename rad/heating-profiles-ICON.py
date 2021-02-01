import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

infile = 'ICON-NWP_prp_AMIP_free_3d_fluxes_heatingrates_mm.nc'
# Heat capacity [J kg-1 K-1]
#cp = 1.08*10**(3)
# Volumetric heat capacity as ICON evaluates on model levels, not pl
cv = 0.718*10**3
# Gravity [m s-2]
g = 9.8
# Array of heating rates
# 60 monthly means, 37 levels, and 2 radiative bands.
H = np.zeros((60,37,2))

# tag is how to label the output npy.
# latmin - latmax is the latitudinal band to extract.
def meanProfile(tag, latmin, latmax):
    flx = xr.open_dataset(infile)
    lats = flx.lat.values
    ii = np.argwhere((lats > latmin) & (lats < latmax))[:,0]
    for i in np.arange(60):
        print(i)
        # The factor of -1 is because pressure decreases upward
        temp = -1.*(flx.lwflxall.isel(time=i,lat=slice(ii[0],ii[-1])) -
                    flx.lwflxclr.isel(time=i,lat=slice(ii[0],ii[-1])))
        temp = temp.assign_coords(lev=temp.lev.values)
        temp = temp.differentiate('lev')
        # Factor of 86400 to convert K s-1 to K day-1.
        #vals = g/cp*86400*temp.mean(dim={'lat','lon'}).values
        #print(vals)
        #time.sleep(5)
        H[i,:,0] = g/cv*86400*temp.mean(dim={'lat','lon'}).values

        temp = -1.*(flx.swflxall.isel(time=i,lat=slice(ii[0],ii[-1]))-
                    flx.swflxclr.isel(time=i,lat=slice(ii[0],ii[-1])))
        #temp = temp.assign_coords(plev=temp.lev.values)
        temp = temp.differentiate('lev')
        H[i,:,1] = g/cv*86400*temp.mean(dim={'lat','lon'}).values
    print('Saving mean heating rates from ' + infile + '...')
    np.save('H_' + tag + '.npy',H)
    return H

#HICON = meanProfile('ICON-AMIP_mm', -90, 90)
HICON_trop = meanProfile('ICON-AMIP_tropic_mm', -30, 30)
HICON = meanProfile('ICON-AMIP_mm',-90,90)
#np.save('lev_ICON-AMIP_mm.npy',xr.open_dataset(infile).lev.values)
sys.exit()

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

#fig.savefig('../output/lw-sw-heating-profiles_1mom+2mom+ERA5.pdf')
plt.show()


fig = plt.figure(figsize=(5.5,6.5))
plt.plot(np.nanmedian(H_ERA5[:,0]+H_ERA5[:,1],axis=0),pl/100,color='black',linestyle='--',label='ERA5')
plt.plot(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0),pl/100,color='red',label='ICON-1mom')
plt.plot(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0),pl/100,color='green',label='ICON-2mom')
plt.plot(np.nanmedian(H_no2mom[:,0]+H_no2mom[:,1],axis=0),pl/100,color='blue',label='ICON-no2mom')
plt.plot(np.nanmedian(H_novgrid[:,0]+H_novgrid[:,1],axis=0),pl/100,color='gold',label='ICON-novgrid')
plt.plot(np.nanmedian(H_rad2mom[:,0]+H_rad2mom[:,1],axis=0),pl/100,color='purple',label='ICON-rad2mom')
#plt.plot(np.nanmedian(H_radnovgrid[:,0]+H_radnovgrid[:,1],axis=0),pl/100,color='purple',linestyle='--',label='ICON-radnovgrid')
plt.plot(np.nanmedian(H_PDA[:,0]+H_PDA[:,1],axis=0),pl/100,color='pink',label='ICON-ACI')
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
print(np.nanmax(np.nanmedian(H_1mom[:,0]+H_1mom[:,1],axis=0)))
print(np.nanmax(np.nanmedian(H_rad2mom[:,0]+H_rad2mom[:,1],axis=0)))
print(np.nanmax(np.nanmedian(H_2mom[:,0]+H_2mom[:,1],axis=0)))
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
fig.savefig('../output/heating-profiles_all.pdf')
plt.show()


