import numpy as np
import xarray as xr
import time, sys
import matplotlib.pyplot as plt
from z_from_ml import z_from_ml

# Read in the Manus ARM site data for August between 1999 and 2010.
# TWP C1, 2.05S and 147.42E
w_8km_Manus = []
w_10km_Manus = []
w_12km_Manus = []
basedir = '/pf/b/b380873/Kalesse-updrafts-ARM/Manus/kalesse/daily_netcdf_manus/1hr/'
for i in np.arange(1999,2011):
    if i == 2002:
       continue
    else:
       twp = xr.open_dataset(basedir + str(i) + '/' + str(i) + '08/twpcidyn1hC1.c1.' + str(i) + '08.nc')
       w = twp['V_air_m']
       h = twp['height']  # 133 heights
       h1 = np.argmin(np.abs(h.values - 8))
       h2 = np.argmin(np.abs(h.values - 10))
       h3 = np.argmin(np.abs(h.values - 12))

       w_8km_Manus.extend(w[:,h1])
       w_10km_Manus.extend(w[:,h2])
       w_12km_Manus.extend(w[:,h3])

# Read in the updraft velocities from the simulation.
basedir = '/work/bb1131/b380873/tropic_run2_output/'
w_8km_ICON = []
w_10km_ICON = []
w_12km_ICON = []
for i in np.arange(51,62):
    icon = xr.open_dataset(basedir + 'W_3D_icon_tropic_00' + str(i) + '_remapdis_global0.025.nc')
    lat = icon['lat'].values
    lon = icon['lon'].values
    j = np.argwhere((lat > 2.) & (lat < 2.5))      # 20 corresponding lats
    k = np.argwhere((lon > 147.0) & (lon < 147.5)) # 20 corresponding lons
    w = icon['w']

    wsub = w[0,:,j[:,0],k[:,0]]
    zh = z_from_ml(lat[j[0,0]],lon[k[0,0]])/1000.
    h1 = np.argmin(np.abs(zh - 8))
    h2 = np.argmin(np.abs(zh - 10))
    h3 = np.argmin(np.abs(zh - 12))

    w_8km_ICON.extend(wsub[h1])
    w_10km_ICON.extend(wsub[h2])
    w_12km_ICON.extend(wsub[h3])

w_8km_Manus = np.asarray(w_8km_Manus)
w_10km_Manus = np.asarray(w_10km_Manus)
w_12km_Manus = np.asarray(w_12km_Manus)

w_8km_Manus = w_8km_Manus[w_8km_Manus > -10.]
w_10km_Manus = w_10km_Manus[w_10km_Manus > -10.]
w_12km_Manus = w_12km_Manus[w_12km_Manus > -10.]

print(len(w_8km_Manus),np.nanmin(w_8km_Manus),np.nanmean(w_8km_Manus),np.nanmax(w_8km_Manus))  
# len = 5040, <-10 = 4776, min/mean/max = (-0.12472, 0.01025, 0.19188)
print(len(w_10km_Manus),np.nanmin(w_10km_Manus),np.nanmean(w_10km_Manus),np.nanmax(w_10km_Manus))
# len = 5040, <-10 = 4497, min/mean/max = (-0.09257, 0.00429, 0.36036)
print(len(w_12km_Manus),np.nanmin(w_12km_Manus),np.nanmean(w_12km_Manus),np.nanmax(w_12km_Manus))
# len = 5040, <-10 = 4466, min/mean/max = (-0.36263, 3.566e-7, 0.21358)

w_8km_ICON = np.asarray(w_8km_ICON)
w_10km_ICON = np.asarray(w_10km_ICON)
w_12km_ICON = np.asarray(w_12km_ICON)

w_8km_ICON = w_8km_ICON[w_8km_ICON > -10.]
w_10km_ICON = w_10km_ICON[w_10km_ICON > -10.]
w_12km_ICON = w_12km_ICON[w_12km_ICON > -10.]

fields = [w_8km_Manus,w_10km_Manus,w_12km_Manus,w_8km_ICON,w_10km_ICON,w_12km_ICON]
let = ['(a)  8 km','(b)  10 km','(c)  12 km']

fig, ax = plt.subplots(nrows=3,ncols=1,figsize=(8,9))
u = 0.3; d = -0.25; n = 50
fs = 13
for i in np.arange(3):
    wgts = np.ones_like(fields[i])/len(fields[i])*100.
    ax[i].hist(fields[i],bins=np.linspace(d,u,n),weights=wgts,color='red',edgecolor='black',alpha=0.5,label='MMCR')
    wgts = np.ones_like(fields[i+3])/len(fields[i+3])*100.
    ax[i].hist(fields[i+3],bins=np.linspace(d,u,n),weights=wgts,color='blue',edgecolor='black',alpha=0.5,label='ICON')
    ax[i].set_yscale('log')
    #ax[i].set_xscale('symlog')
    ax[i].set_ylim([10**(-2),100])
    ax[i].set_ylabel('Probability [%]',fontsize=fs)
    ax[i].text(0.05,0.9,let[i],weight='bold',color='black',transform=ax[i].transAxes,fontsize=fs+2)
    if i == 0:
       ax[i].legend(loc='upper right',fontsize=fs-3)
    if i == 2:
       ax[i].set_xlabel(r'Vertical velocity [m s^{-1}]',fontsize=fs)


fig.savefig('MMCR_ICON_vert_velocity.pdf',bbox_inches='tight')
plt.show()
