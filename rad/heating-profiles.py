import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8
# Hold all the heating rates here, 18 pressure levels x 4 simulations x (LW+SW) x 2 times
H = np.zeros((16,18))

for i in np.arange(2):
    # .where(clch_1mom > 0.5)
    basedir = '/work/bb1131/b380873/tropic_run2_output/reproduce/'
    flx_1mom = xr.open_dataset(basedir + 'RAD3D_13_19_0.025deg_PL.nc')
    clch_1mom = xr.open_dataset(basedir + 'CLCONV_2D_26_38_0.025deg.nc')
    lwcld_1mom = (flx_1mom.lwflxclr.isel(time=i)-flx_1mom.lwflxall.isel(time=i)).mean(dim={'lat','lon'})
    swcld_1mom = (flx_1mom.swflxclr.isel(time=i)-flx_1mom.swflxall.isel(time=i)).mean(dim={'lat','lon'})

    # .where(clch_2mom > 0.5)
    basedir = '/work/bb1131/b380873/tropic_run5_output/'
    flx_2mom = xr.open_dataset(basedir + 'RAD3D_61_67_0.025deg_PL.nc')
    clch_2mom = xr.open_dataset(basedir + 'CLCONV_2D_122_138_0.025deg.nc').clch.isel(time=i,lev=0)
    lwcld_2mom = (flx_2mom.lwflxclr.isel(time=i)-flx_2mom.lwflxall.isel(time=i)).mean(dim={'lat','lon'})
    swcld_2mom = (flx_2mom.swflxclr.isel(time=i)-flx_2mom.swflxall.isel(time=i)).mean(dim={'lat','lon'})

    # .where(clch_novgrid > 0.5)
    basedir = '/work/bb1131/b380873/tropic_run5_output/novgrid/'
    flx_novgrid = xr.open_dataset(basedir + 'RAD3D_01_07_0.025deg_PL.nc')
    clch_novgrid = xr.open_dataset(basedir + 'CLCONV_2D_02_14_0.025deg.nc').clch.isel(time=i,lev=0)
    lwcld_novgrid = (flx_novgrid.lwflxclr.isel(time=i)-flx_novgrid.lwflxall.isel(time=i)).mean(dim={'lat','lon'})
    swcld_novgrid = (flx_novgrid.swflxclr.isel(time=i)-flx_novgrid.swflxall.isel(time=i)).mean(dim={'lat','lon'})

    basedir = '/work/bb1131/b380873/tropic_run5_output/no2mom/'
    flx_no2mom = xr.open_dataset(basedir + 'RAD3D_01_07_0.025deg_PL.nc')
    clch_no2mom = xr.open_dataset(basedir + 'CLCONV_2D_02_14_0.025deg.nc').clch.isel(time=i,lev=0)
    lwcld_no2mom = (flx_no2mom.lwflxclr.isel(time=i)-flx_no2mom.lwflxall.isel(time=i)).mean(dim={'lat','lon'})
    swcld_no2mom = (flx_no2mom.swflxclr.isel(time=i)-flx_no2mom.swflxall.isel(time=i)).mean(dim={'lat','lon'})

    # Multiply by 86400 to convert [K s-1] to [K day-1]
    H[i] = g/cp*np.gradient(lwcld_1mom,flx_1mom.plev_2)*86400.
    H[i+2] = g/cp*np.gradient(lwcld_2mom,flx_2mom.plev_2)*86400.
    H[i+4] = g/cp*np.gradient(lwcld_novgrid,flx_novgrid.plev_2)*86400.
    H[i+6] = g/cp*np.gradient(lwcld_no2mom,flx_no2mom.plev_2)*86400.
    H[i+8] = g/cp*np.gradient(swcld_1mom,flx_1mom.plev_2)*86400.
    H[i+10] = g/cp*np.gradient(swcld_2mom,flx_2mom.plev_2)*86400.
    H[i+12] = g/cp*np.gradient(swcld_novgrid,flx_novgrid.plev_2)*86400.
    H[i+14] = g/cp*np.gradient(swcld_no2mom,flx_no2mom.plev_2)*86400.


fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(8,5))
ax[0].plot(H[0],flx_1mom.plev_2/100,color='red',label='ICON-1mom')
ax[0].plot(H[2],flx_2mom.plev_2/100,color='green',label='ICON-2mom')
ax[0].plot(H[4],flx_no2mom.plev_2/100,color='blue',label='ICON-no2mom')
ax[0].plot(H[6],flx_novgrid.plev_2/100,color='gold',label='ICON-novgrid')
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'LW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].set_ylim([50,800])
ax[0].set_yscale('log')
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)
ax[0].legend(fontsize=fs-3,loc='center left')

ax[1].plot(H[8],flx_1mom.plev_2/100,color='red',label='ICON-1mom')
ax[1].plot(H[10],flx_2mom.plev_2/100,color='green',label='ICON-2mom')
ax[1].plot(H[12],flx_no2mom.plev_2/100,color='blue',label='ICON-no2mom')
ax[1].plot(H[14],flx_novgrid.plev_2/100,color='gold',label='ICON-novgrid')
#ax[1].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[1].set_xlabel(r'SW cloud radiative heating [K day$^{-1}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)
ax[1].set_ylim([50,800])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)
fig.suptitle('Domain-mean 8 Aug 2017 01:00 UTC',fontsize=fs)

fig.savefig('./output/lw-sw-heating-profiles.pdf')
plt.show()

fs = 14
fig = plt.figure(figsize=(8,4.5))
#plt.plot(H[1]+H[9],flx_1mom.plev_2/100,color='red',label='ICON-1mom')
#plt.plot(H[3]+H[11],flx_2mom.plev_2/100,color='green',label='ICON-2mom')
#plt.plot(H[5]+H[13],flx_no2mom.plev_2/100,color='blue',label='ICON-no2mom')
#plt.plot(H[7]+H[15],flx_novgrid.plev_2/100,color='gold',label='ICON-novgrid')
plt.plot(H[0]+H[8],flx_1mom.plev_2/100,color='red',label='ICON-1mom')
plt.plot(H[2]+H[10],flx_2mom.plev_2/100,color='green',label='ICON-2mom')
plt.plot(H[4]+H[12],flx_no2mom.plev_2/100,color='blue',label='ICON-no2mom')
plt.plot(H[6]+H[14],flx_novgrid.plev_2/100,color='gold',label='ICON-novgrid')
plt.plot([0,0],[50,1000],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Cloud radiative heating [K day$^{-1}$]',fontsize=fs)
#plt.text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
plt.title('Domain-mean 8 Aug 2017 01:00 UTC',fontsize=fs)
plt.ylim([50,1000])
plt.gca().invert_yaxis()
plt.gca().tick_params(labelsize=fs)
plt.gca().legend(fontsize=fs-1)
fig.savefig('./output/heating-profiles2.pdf')
plt.show()


