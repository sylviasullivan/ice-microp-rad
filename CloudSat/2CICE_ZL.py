# conda activate pyhdf
from CloudSat_read import read_cloudsatcalipso_hdf_file
from matplotlib import cm
from scipy.interpolate import CubicSpline
import xarray as xr
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime,timedelta
import sys,time

basedir = '/work/bb1018/b380873/tropic_vis/obs/2C-ICE/'
writenpy = True
# Which model timesteps are being used for comparison?
dt1 = '0043'
dt2 = '0066'

# All CloudSat filenames are at the bottom.
fi1 = '2017219064524_59987_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
fi2 = '2017220054945_60001_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
var = 'IWC'

# To read in effective radius, replace iwc with re and 'var' with 're'.
# Extinction is also available in the CloudSat files.
# Longitudes run from 180 to -180. lon and lat are arrays.
# dat = (37082,125), lon/lat = (37082,1)
iwc1, h1, lon1, lat1, zeit1 = read_cloudsatcalipso_hdf_file(basedir + fi1,var)
iwc2, h2, lon2, lat2, zeit2 = read_cloudsatcalipso_hdf_file(basedir + fi2,var)
#print(type(zeit1))
#print(type(lon1))

# Find the data in the Asian monsoon region.
ii1 = np.argwhere((lon1 > 70) & (lon1 <= 100) & (lat1 > 5) & (lat1 <= 35))
zeit1 = np.asarray(zeit1)
t_AMA1 = zeit1[ii1[:,0]]
ii2 = np.argwhere((lon2 > 70) & (lon2 <= 100) & (lat2 > 5) & (lat2 <= 35))
zeit2 = np.asarray(zeit2)
t_AMA2 = zeit2[ii2[:,0]]

# Create datetimes from the t_AMA values above. Values are seconds after start.
startyr = datetime(2017,1,1,0,0,0)
startda = startyr + timedelta(days=int(fi1[4:7]),hours=int(fi1[7:9]),minutes=int(fi1[9:11]),seconds=int(fi1[11:13]))
da1 = np.array([startda + timedelta(seconds=i[0]) for i in t_AMA1])
startda = startyr + timedelta(days=int(fi2[4:7]),hours=int(fi2[7:9]),minutes=int(fi2[9:11]),seconds=int(fi2[11:13]))
da2 = np.array([startda + timedelta(seconds=i[0]) for i in t_AMA2])

h_AMA1   = h1[ii1[:,0]]   # (3111,125)
iwc_AMA1 = iwc1[ii1[:,0]]
ll_AMA1  = np.stack((lat1[ii1[:,0]],lon1[ii1[:,0]]))[:,:,0]   # (2,1918)
h_AMA2   = h2[ii2[:,0]]
iwc_AMA2 = iwc2[ii2[:,0]]
ll_AMA2  = np.stack((lat2[ii2[:,0]],lon2[ii2[:,0]]))[:,:,0]

# Print out typical bounds of the data.
print('Time 1 min, max: ' + str(np.nanmin(t_AMA1)) + ' ' + str(np.nanmax(t_AMA1)))
print('IWC 1 min, mean, max: ' + str(np.nanmin(iwc_AMA1)) + ' ' + str(np.nanmean(iwc_AMA1)) + ' ' + str(np.nanmax(iwc_AMA1)))
print('Time 2 min, max: ' + str(np.nanmin(t_AMA2)) + ' ' + str(np.nanmax(t_AMA2)))
print('IWC 2 min, mean, max: ' + str(np.nanmin(iwc_AMA2)) + ' ' + str(np.nanmean(iwc_AMA2)) + ' ' + str(np.nanmax(iwc_AMA2)))

# Create the coordinate mesh and convert from m to km for the height coordinate.
xx1, yy1 = np.meshgrid(np.arange(ii1.shape[0]),h_AMA1[0])
yy1 = yy1/1000.
xx2, yy2 = np.meshgrid(np.arange(ii2.shape[0]),h_AMA2[0])
yy2 = yy2/1000.

fs = 12
fig, ax = plt.subplots(nrows=4,ncols=1,figsize=(14,6.5))
# Plot only non-zero values.
iwc_AMA1[iwc_AMA1 <= 0] = np.nan
levels = np.logspace(-3,1,10)
c = ax[0].contourf(xx1,yy1,iwc_AMA1.T,norm=colors.LogNorm(),levels=levels,cmap=cm.Blues)
ax[0].set_ylabel('Height [km]',fontsize=fs)
ax[0].set_ylim([2,20])
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[0],ticks=[0.001,0.01,0.1,1,10])

iwc_AMA2[iwc_AMA2 <= 0] = np.nan
c = ax[2].contourf(xx2,yy2,iwc_AMA2.T,norm=colors.LogNorm(),levels=levels,cmap=cm.Greens)
ax[2].set_ylabel('Height [km]',fontsize=fs)
ax[2].set_ylim([2,20])
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[2],ticks=[0.001,0.01,0.1,1,10])

# Pull the simulated profile 2017-08-07 7:00:00 UTC (43 hours after start + 5 hour time diff)
if writenpy == True:
   basedir = '/work/bb1131/b380873/tropic_run5_output/'
   fi_icon1 = xr.open_dataset(basedir + 'CLCONV_3D_icon_tropic_' + dt1 + '_remapdis_0.025_HL.nc')
   fi_icon2 = xr.open_dataset(basedir + 'CLCONV_3D_icon_tropic_' + dt2 + '_remapdis_0.025_HL.nc')
   lat_icon1 = fi_icon1.lat.values
   lon_icon1 = fi_icon1.lon.values
   h_AMA1 = fi_icon1.height
   np.save('height_0043_0066.npy',h_AMA1)

   # Make sure the modeled and observed lats and lons are the same range.
   qi_icon1 = fi_icon1.qi
   qs_icon1 = fi_icon1.qs
   qg_icon1 = fi_icon1.qg
   qi_icon2 = fi_icon2.qi
   qs_icon2 = fi_icon2.qs
   qg_icon2 = fi_icon2.qg

   qice_icon_grid = np.zeros((2,iwc_AMA1.shape[0],120))
   for l in np.arange(ll_AMA1.shape[1]):
       print(l)
       lat1 = ll_AMA1[0,l]; lon1 = ll_AMA1[1,l]
       # Extract the value for which the diff between satellite and modeled lat / lon is least.
       v1 = qi_icon1.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       v2 = qs_icon1.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       v3 = qg_icon1.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       qice_icon_grid[0,l] = v1 + v2 + v3

   for l in np.arange(ll_AMA2.shape[1]):
       print(l)
       lat1 = ll_AMA2[0,l]; lon1 = ll_AMA2[1,l]
       # Extract the value for which the diff between satellite and modeled lat / lon is least.
       v1 = qi_icon2.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       v2 = qs_icon2.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       v3 = qg_icon2.sel(time=0,lat=lat1,lon=lon1,method='nearest')
       qice_icon_grid[1,l] = v1 + v2 + v3

   # 1000/1.225 converts kg ice / kg air --> g ice / m3 air
   qice_icon_grid = qice_icon_grid*1000/1.225
   # Filter out noise.
   qice_icon_grid = np.where(qice_icon_grid > 10**(-7),qice_icon_grid,np.nan)
   print(qice_icon_grid)
   np.save('qi+qs+qg_icon0043_icon0066_59987_60001.npy',qice_icon_grid)
else:
   h_AMA1 = np.load('height_0043_0066.npy')
   qice_icon_grid = np.load('qi+qs+qg_icon0043_icon0066_59987_60001.npy')


xx1, yy1 = np.meshgrid(np.arange(qice_icon_grid.shape[1]),h_AMA1)
yy1 = yy1/1000.
c = ax[1].contourf(xx1,yy1,qice_icon_grid[0].T,cmap=cm.Blues,norm=colors.LogNorm(),levels=levels)
ax[1].set_ylabel('Height [km]',fontsize=fs)
ax[1].set_xlabel('CloudSat scan number',fontsize=fs)
ax[1].set_ylim([2,20])
#ax[1].invert_yaxis()
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[1],ticks=[0.001,0.01,0.1,1,10])

xx2, yy2 = np.meshgrid(np.arange(qice_icon_grid.shape[1]),h_AMA1)
yy2 = yy2/1000.
c = ax[3].contourf(xx2,yy2,qice_icon_grid[1].T,cmap=cm.Greens,norm=colors.LogNorm(),levels=levels)
ax[3].set_ylabel('Height [km]',fontsize=fs)
ax[3].set_xlabel('CloudSat scan number',fontsize=fs)
ax[3].set_xlim([0,500])
ax[3].set_ylim([2,20])
#ax[3].invert_yaxis()
plt.colorbar(c,label=r'IWC [g m$^{-3}$]',ax=ax[3],ticks=[0.001,0.01,0.1,1,10])

#fig.savefig('CS_' + fi[:20] + 'ICON' + modeltimestep + '_comparison_ZL.pdf',bbox_inches='tight')
fig.savefig('CloudSat_ICON_comparison.pdf',bbox_inches='tight')
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
