# conda activate pyhdf
from CloudSat_read import read_cloudsatcalipso_hdf_file
from z_from_ml import z_from_ml
from matplotlib import cm
from scipy.interpolate import CubicSpline
import xarray as xr
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime,timedelta
import sys,time

basedir = 'obs/2C-ICE/'
writenpy = True
#modeltimestep = '0043'
#modeltimestep = '0066'
modeltimestep = '0067'  #'0068'
#modeltimestep = '0079'
# All file names at bottom
#fi = '2017219064524_59987_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
#fi = '2017220054945_60001_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
fi = '2017220072838_60002_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
#fi = '2017220190049_60009_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
var = 'IWC'
var2 = 're'
var3 = 'EXT_coef'

iwc, h, lon, lat, zeit = read_cloudsatcalipso_hdf_file(basedir + fi,var)
# lon runs from 180 to -180
re, _, _, _, _ = read_cloudsatcalipso_hdf_file(basedir + fi,var2)
ext,_, _, _, _ = read_cloudsatcalipso_hdf_file(basedir + fi,var3)
# lon, lat are arrays
# dat = (37082,125), lon/lat = (37082,1)

# Find the data in the Asian monsoon region.
ii = np.argwhere((lon > 70) & (lon <= 100) & (lat > 5) & (lat <= 35))
zeit    = np.asarray(zeit)
t_AMA   = zeit[ii[:,0]]

# Create datetimes from the t_AMA values above, assuming they are seconds after start.
startyr = datetime(2017,1,1,0,0,0)
startda = startyr + timedelta(days=int(fi[4:7]),hours=int(fi[7:9]),minutes=int(fi[9:11]),seconds=int(fi[11:13]))
da = np.array([startda + timedelta(seconds=i[0]) for i in t_AMA])

h_AMA   = h[ii[:,0]]   # (3111,125)
iwc_AMA = iwc[ii[:,0]]
re_AMA  = re[ii[:,0]]
ext_AMA = ext[ii[:,0]]
ll_AMA  = np.stack((lat[ii[:,0]],lon[ii[:,0]]))[:,:,0]   # (2,1918)

# Print out typical bounds of the data.
print('Time min, max: ' + str(np.nanmin(t_AMA)) + ' ' + str(np.nanmax(t_AMA)))
print('IWC min, mean, max: ' + str(np.nanmin(iwc_AMA)) + ' ' + str(np.nanmean(iwc_AMA)) + ' ' + str(np.nanmax(iwc_AMA)))
print('reff min, mean, max: ' + str(np.nanmin(re_AMA)) + ' ' + str(np.nanmean(re_AMA)) + ' ' + str(np.nanmax(re_AMA)))

# Create the coordinate mesh and convert from m to km for the height coordinate.
xx, yy  = np.meshgrid(np.arange(ii.shape[0]),h_AMA[0])
yy = yy/1000.

fs = 12
fig, ax = plt.subplots(nrows=4,ncols=1,figsize=(14,6.5))
# Plot only non-zero values.
iwc_AMA_pos = iwc_AMA
iwc_AMA_pos[iwc_AMA_pos <= 0] = np.nan
levels = np.logspace(-3,1,10)
c = ax[0].contourf(xx,yy,iwc_AMA_pos.T,norm=colors.LogNorm(),levels=levels,cmap=cm.Blues)
ax[0].set_ylabel('Height [km]',fontsize=fs)
ax[0].set_ylim([2,20])
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[0],ticks=[0.001,0.01,0.1,1,10])

re_AMA_pos = re_AMA
re_AMA_pos[re_AMA_pos <= 0] = np.nan
c = ax[2].contourf(xx,yy,re_AMA_pos.T,cmap=cm.Reds)
ax[2].set_ylabel('Height [km]',fontsize=fs)
ax[2].set_ylim([2,20])
plt.colorbar(c,label=r'$r_e$ [$\mu$m]',ax=ax[2])

ext_AMA_pos = ext_AMA
ext_AMA_pos[ext_AMA_pos <= 0] = np.nan
c = ax[3].contourf(xx,yy,ext_AMA_pos.T,cmap=cm.Greens,norm=colors.LogNorm())
ax[3].set_ylabel('Height [km]',fontsize=fs)
ax[3].set_xlabel('Scan number',fontsize=fs)
ax[3].set_ylim([2,20])
plt.colorbar(c,label=r'$\epsilon_i$',ax=ax[3])

# Pull the simulated profile 2017-08-07 7:00:00 UTC (43 hours after start + 5 hour time diff)
if writenpy == True:
   fi_icon = xr.open_dataset('../tropic_run2_output/CLCONV_3D_icon_tropic_' + modeltimestep + '_remapdis_global0.025.nc')
   lat_icon = fi_icon.lat
   lon_icon = fi_icon.lon
   # Make sure the modeled and observed lats and lons are the same range.
   #print(fi_icon.lat.values,np.nanmax(fi_icon.lat.values))
   qi_icon = fi_icon.qi.values
   qi_icon_grid = np.zeros((iwc_AMA.shape[0],105))  # iwc_AMA.shape
   ll_ICON = np.zeros((ii.shape[0],2))
   for l in np.arange(ll_AMA.shape[1]):   # 3111
       print(l)
       lat = ll_AMA[0,l]; lon = ll_AMA[1,l]
       alt = z_from_ml(lat,lon)
       # where is the diff between satellite and modeled lat least?
       oo  = np.abs(lat - lat_icon).argmin()
       # where is the diff between satellite and modeled lon least?
       jj = np.abs(lon - lon_icon).argmin()
       # Store these simulated lat / lon pairs.
       ll_ICON[l,0] = lat_icon.values[oo]
       ll_ICON[l,1] = lon_icon.values[jj]
       #mm = np.argwhere(h_AMA[l] < 0) # mm[0,0] is generally 105
       spl = CubicSpline(np.flip(alt),np.flip(qi_icon[0,:105,oo,jj])) #,k=3)
       qi_icon_grid[l] = spl(np.flip(h_AMA[l,:105]))
       #qi_icon_grid[l] = np.interp(h_AMA[l,:105],alt,qi_icon[0,:105,oo,jj])
       #qi_icon_grid[l] = qi_icon[0,:,oo,jj]

   # 1000/1.225 converts kg ice / kg air --> g ice / m3 air
   qi_icon_grid = qi_icon_grid*1000/1.225
   qi_icon_pos = qi_icon_grid
   # Filter out noise.
   qi_icon_pos[qi_icon_pos <= 10**(-7)] = np.nan
   np.save('qi_icon' + modeltimestep + '_CS' + fi[:19] + '_interp.npy',qi_icon_pos)
else:
   qi_icon_pos = np.load('qi_icon' + modeltimestep + '_CS' + fi[:19] + '_interp.npy')

xx1, yy1 = np.meshgrid(np.arange(ii.shape[0]),np.arange(105))
yy1 = yy1/1000.
c = ax[1].contourf(xx1,yy1,qi_icon_pos.T,cmap=cm.Blues,norm=colors.LogNorm(),levels=levels)
ax[1].set_ylabel('Model level',fontsize=fs)
ax[1].set_xlabel('CloudSat scan number',fontsize=fs)
ax[1].set_ylim([2,20])
#ax[1].invert_yaxis()
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[1],ticks=[0.001,0.01,0.1,1,10])

fig.savefig('CS_' + fi[:20] + 'ICON' + modeltimestep + '_comparison_ML.pdf',bbox_inches='tight')
plt.show()
sys.exit()

fig2, ax = plt.subplots(nrows=1,ncols=2,figsize=(7,4))
ax[0].scatter(ll_ICON[:,0],ll_AMA[0,:],s=5)
ax[0].plot(np.linspace(17,35,100),np.linspace(17,35,100),color='k',linewidth=0.75)
ax[0].set_ylabel('CloudSat latitude')
ax[0].set_xlabel('ICON latitude')
ax[0].set_xlim([17,35]); ax[0].set_ylim([17,35])

ax[1].scatter(ll_ICON[:,1],ll_AMA[1,:],s=5)
ax[1].plot(np.linspace(95,100,100),np.linspace(95,100,100),color='k',linewidth=0.75)
ax[1].set_ylabel('CloudSat longitude')
ax[1].set_xlabel('ICON longitude')
ax[1].set_xlim([95,100]); ax[1].set_ylim([95,100])

#fig2.savefig('ll_CloudSat_ICON_validation.pdf',bbox_inches='tight')
plt.show()

# All file names.
#2017218060209_59972_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf  2017219181735_59994_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf
#2017218074102_59973_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf  2017220054945_60001_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf
#2017218173420_59979_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf  2017220072838_60002_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf
#2017218191313_59980_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf  2017220172156_60008_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf
#2017219064524_59987_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf  2017220190049_60009_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf
