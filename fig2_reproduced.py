import matplotlib.pyplot as plt
from netCDF4 import Dataset,num2date
from datetime import datetime
from functools import reduce
import numpy.ma as ma
import scipy.interpolate
import numpy as np
import xarray as xr
import sys,time,os
from z_from_ml import z_from_ml

#full_level_heights = np.asarray([73681,71102,68600,66152,63775,61470,59233,57054,\
#      54927,52851,50825,48858,46952,45106,43318,41587,39912,38293,36727,35215,33755,\
#      32346,30988,29680,28421,27210,26053,24952,23912,22934,22015,21150,20336,19572,\
#      18853,18176,17540,16942,16380,15851,15354,14886,14447,14027,13621,13221,12821,\
#      12421,12021,11621,11221,10821,10421,10021,9621,9221,8821,8421,8021,7621,7221,\
#      6821,6421,6021,5621,5227,4846,4480,4127,3788,3462,3151,2853,2570,2300,2043,\
#      1801,1573,1358,1158,972,801,644,503,377,267,174,98,42,10])

# top_height parameter in run_TROPIC.run was 30 km.
full_level_heights = np.asarray([30988,29680,28421,27210,26053,24952,23912,22934,\
      22015,21150,20336,19572,18853,18176,17540,16942,16380,15851,15354,14886,14447,\
      14027,13621,13221,12821,12421,12021,11621,11221,10821,10421,10021,9621,9221,\
      8821,8421,8021,7621,7221,6821,6421,6021,5621,5227,4846,4480,4127,3788,3462,3151,\
      2853,2570,2300,2043,1801,1573,1358,1158,972,801,644,503,377,267,174,98,42,10])

# Small function to reorder the elements of arr by the indices of indx.
def reorder(arr,indx):
    print('Original array shape: ' + str(arr.shape))
    temp = np.zeros((arr.shape[0],))
    for i in np.arange(arr.shape[0]):
        temp[i] = arr[indx[i]]
    print('Reordered array shape: ' + str(temp.shape))
    return temp

daten = Dataset('obs/stratoclim2017.geophysika.0808_1.master.ci_eval.nc','r')
zeit = num2date(daten.variables['time'][:],'seconds since 2000-01-01 00:00:0.0')

# In Figure 2 of Lee et al., only values between 6:20 and 6:48 UTC are used.
# These indices correspond exactly to 2017-08-08 06:20:00 and 2017-08-08 06:48:01.
ii = np.argwhere(zeit >= datetime(2017,8,8,6,20))[0,0]
jj = np.argwhere(zeit > datetime(2017,8,8,6,48))[0,0]

# Vertical profile of water vapour [ppmv] for first panel.
# Filter for lack of altitude data, recorded as -99999.
alt = np.asarray(daten.variables['BEST:ALT'][ii:jj])

# Some weirdness here using 0 versus O in the variable name.
# Filter for lack of altitude data, recorded as -99999.
h2o_flash = np.asarray(daten.variables['BEST:H2O_gas'][ii:jj])
h2o_fish = np.asarray(daten.variables['BEST:H2O_enh'][ii:jj])
iwc = np.asarray(daten.variables['BEST:IWC'][ii:jj])
temp = np.asarray(daten.variables['BEST:TEMP'][ii:jj])
theta = np.asarray(daten.variables['BEST:THETA'][ii:jj])
rhice_flash = np.asarray(daten.variables['BEST:RH_ice_gas'][ii:jj])
rhice_fish = np.asarray(daten.variables['BEST:RH_ice_enh'][ii:jj])

# Filter for positive values of vapor mixing ratio.
i1 = np.argwhere(alt > 0)[:,0]
i2 = np.argwhere(h2o_flash > 0)[:,0]
i3 = np.argwhere(h2o_fish > 0)[:,0]
indx = reduce(np.intersect1d,(i1,i2,i3))
alt_filtered1 = alt[indx]
h2o_flash_filtered = h2o_flash[indx]
h2o_fish_filtered = h2o_fish[indx]

# Filter for positive values of IWC.
i2 = np.argwhere(iwc > 0)[:,0]
indx = reduce(np.intersect1d,(i1,i2))
alt_filtered2 = alt[indx]
iwc_filtered = iwc[indx]

# Filter for positive temperatures.
i2 = np.argwhere((temp > 0) & (theta > 0))[:,0]
indx = reduce(np.intersect1d,(i1,i2))
alt_filtered3 = alt[indx]
temp_filtered = temp[indx]
theta_filtered = theta[indx]

# Filter for positive RHice values.
i2 = np.argwhere((rhice_flash > 0) & (rhice_fish > 0))
indx = reduce(np.intersect1d,(i1,i2))
alt_filtered4 = alt[indx]
rhice_flash_filtered = rhice_flash[indx]
rhice_fish_filtered = rhice_fish[indx]

# Order the altitudes from lowest to highest and align the h2o_flash and h2o_fish.
#alt_sort = np.sort(alt_filtered)
#alt_sort_i = np.argsort(alt_filtered)
#h2o_flash_sort = reorder(h2o_flash,alt_sort_i)

# If binning between up and down with <binnum> bins, which elements go in which bin?
# Make a multidimensional list of alt and h2o values in each. 
binnum = 100
up = 14000
down = 22000
indx1 = np.digitize(alt_filtered1,bins=np.linspace(up,down,binnum))
indx2 = np.digitize(alt_filtered2,bins=np.linspace(up,down,binnum))

# All elements
alt_list1 = [[] for i in np.arange(binnum)]
h2o_flash_list = [[] for i in np.arange(binnum)]
h2o_fish_list = [[] for i in np.arange(binnum)]
alt_list2 = [[] for i in np.arange(binnum)]
iwc_list = [[] for i in np.arange(binnum)]
# Their mean
h2o_flash_m = np.zeros((binnum,))
h2o_fish_m = np.zeros((binnum,))
iwc_m = np.zeros((binnum,))
for i in np.arange(binnum):
    indx_subset1 = np.argwhere((indx1 == i))[:,0]
    indx_subset2 = np.argwhere((indx2 == i))[:,0]
    # All elements
    alt_list1[i].extend(alt_filtered1[indx_subset1])
    h2o_flash_list[i].extend(h2o_flash_filtered[indx_subset1])
    h2o_fish_list[i].extend(h2o_fish_filtered[indx_subset1])
    alt_list2[i].extend(alt_filtered2[indx_subset2])
    iwc_list[i].extend(iwc_filtered[indx_subset2])
    # Their mean, only if you have at least 5 elements
    if (len(h2o_flash_list[i]) > 5):
       h2o_flash_m[i] = sum(h2o_flash_list[i]) / len(h2o_flash_list[i])
    if (len(h2o_fish_list[i]) > 5):
       h2o_fish_m[i] = sum(h2o_fish_list[i]) / len(h2o_fish_list[i])
    if (len(iwc_list[i]) > 5):
       iwc_m[i] = sum(iwc_list[i]) / len(iwc_list[i])    

# Output the number of measurements in each bin
#for i in np.arange(binnum):
#    print(len(h2o_flash_list[i]))

# Read in the simulation values. 
basedir = '/scratch/b/b380873/tropic_run2_restart/output_2017080800-2017080806/'
sims = xr.open_dataset(basedir + 'CLCONV_3D_icon_tropic_0066_remap_global0.025.nc')
qi_s = sims.qi
qv_s = sims.qv    # dims are ('time','height','lat','lon')
lat_s = sims.lat  # shape is (1800,) and type is DataArray
lon_s = sims.lon  # shape is (4600,) and type is DataArray

# Filter the simulation specific humidities for 25 to 26.5 N and 85 to 85.5 E.
i1 = np.abs(lat_s - 25).argmin()
lat_i1 = np.unravel_index(i1,lat_s.shape)[0]
i2 = np.abs(lat_s - 26.5).argmin()
lat_i2 = np.unravel_index(i2,lat_s.shape)[0]
i1 = np.abs(lon_s - 85).argmin()
lon_i1 = np.unravel_index(i1,lon_s.shape)[0]
i2 = np.abs(lon_s - 85.5).argmin()
lon_i2 = np.unravel_index(i2,lon_s.shape)[0]

lat_s_filtered = lat_s[lat_i1:lat_i2]
lon_s_filtered = lon_s[lon_i1:lon_i2]
qi_s_filtered = qi_s[0,:,lat_i1:lat_i2,lon_i1:lon_i2]
qv_s_filtered = qv_s[0,:,lat_i1:lat_i2,lon_i1:lon_i2]   # shape (90, 60, 20) dims ('height', 'lat', 'lon')

# Create an array of associated heights if it does not already exist.
if os.path.isfile('z_lat25-26.5_lon85-85.5.npy') == False:
   z = np.zeros((68,60,20))
   for lat in np.arange(60):
       for lon in np.arange(20):
           lala = lat_s_filtered[lat].item()
           lolo = lon_s_filtered[lon].item()
           zprofil = z_from_ml(lala,lolo)
           z[:,lat,lon] = zprofil
   np.save('z_lat25-26.5_lon85-85.5.npy',z)
else:
   z = np.load('z_lat25-26.5_lon85-85.5.npy')

# Interpolate the values to the fixed levels.
qi_s_filtered_regridded = np.zeros((68,60,20))
qv_s_filtered_regridded = np.zeros((68,60,20))
for lat in np.arange(60):
    for lon in np.arange(20):
        qi_s_filtered_regridded[:,lat,lon] = scipy.interpolate.griddata(z[:,lat,lon],\
              qi_s_filtered[:,lat,lon],full_level_heights,method='cubic')
        qv_s_filtered_regridded[:,lat,lon] = scipy.interpolate.griddata(z[:,lat,lon],\
              qv_s_filtered[:,lat,lon],full_level_heights,method='cubic')

# Multiply the spatially-averaged qv_s_filtered_regridded by 10^6 to translate kg kg-1 to ppmv.
mw_dryair = 28.97*1000    # kg air (mol air)-1
mw_watvap = 18.02*1000    # kg wv (mol wv)-1
conv = mw_dryair / mw_watvap
# Where the qi values < 0, set to nan
qi_s_filtered_regridded_masked = ma.masked_where(qi_s_filtered_regridded < 0,qi_s_filtered_regridded)
qv_s_final = np.nanmean(np.nanmean(qv_s_filtered_regridded*conv,2),1)*10**6
qi_s_final = np.nanmean(np.nanmean(qi_s_filtered_regridded_masked*conv,2),1)*10**6

#fig = plt.figure()
#plt.subplot2grid((1,1),(0,0))
#plt.scatter(np.nanmean(np.nanmean(qv_s_filtered_regridded,2),1),full_level_heights)
#plt.show()

fs = 12
fig = plt.figure()
plt.subplot2grid((1,4),(0,0))
plt.grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
plt.scatter(h2o_flash_filtered,alt_filtered1/1000,color='k',s=10,marker='*')
plt.scatter(h2o_fish_filtered,alt_filtered1/1000,color='gray',s=10,marker='d')
plt.scatter(qv_s_final,full_level_heights/1000,marker='o',color='r',s=15)
#plt.plot(h2o_flash_m,np.linspace(up,down,binnum)/1000,linewidth=1.25,color='k')
#plt.plot(h2o_fish_m,np.linspace(up,down,binnum)/1000,linewidth=1.25,color='b')
plt.plot([4,4],[14,21],linewidth=0.75,linestyle='--',color='k')
plt.xlim([2,10])
#plt.ylim([14,21])
plt.xlabel('Water vapour [ppmv]',fontsize=fs)
plt.ylabel('Altitude [km]',fontsize=fs)
plt.tick_params(labelsize=fs)

plt.subplot2grid((1,4),(0,1))  #5
plt.grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
plt.scatter(iwc_filtered,alt_filtered2/1000,color='k',s=10,marker='*')
plt.scatter(qi_s_final,full_level_heights/1000,marker='o',color='r',s=15)
#plt.plot(iwc_m,np.linspace(up,down,binnum)/1000,linewidth=1.25,color='k')
#plt.ylim([14,21])
plt.xlim([0,2])
plt.xlabel('Ice [eq.ppmv]',fontsize=fs)
plt.gca().set_yticklabels([])

ax1 = plt.subplot2grid((1,4),(0,2))  
ax2 = ax1.twiny()
ax1.grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
ax1.scatter(temp_filtered-273,alt_filtered3/1000,color='k',s=10,marker='*')
ax2.scatter(theta_filtered,alt_filtered3/1000,color='gray',s=10,marker='*')
ax1.set_ylim([14,21])
ax1.set_xlim([-85,-65]); ax2.set_xlim([350,430])
ax1.set_xlabel('Temperature [deg C]',fontsize=fs)
ax2.set_xlabel('Potential temperature [K]',fontsize=fs)
ax1.set_yticklabels([])
ax2.set_yticklabels([])

plt.subplot2grid((1,4),(0,3))
plt.grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
plt.scatter(rhice_fish_filtered,alt_filtered4/1000,color='k',s=10,marker='*')
plt.scatter(rhice_flash_filtered,alt_filtered4/1000,color='gray',s=10,marker='*')
plt.ylim([14,21])
plt.xlim([40,130])
plt.xlabel('RHice [%]',fontsize=fs)
plt.gca().set_yticklabels([])

plt.tight_layout()
plt.savefig('fig2_reproduced.pdf',dpi=200,bbox_inches='tight')
plt.show()
