import matplotlib.pyplot as plt
from netCDF4 import Dataset,num2date
from datetime import datetime
from functools import reduce
import numpy.ma as ma
import scipy.interpolate
import numpy as np
import xarray as xr
import sys,time,os

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

daten = Dataset('../obs/stratoclim2017.geophysika.0808_1.master.ci_eval.nc','r')
zeit = num2date(daten.variables['time'][:],'seconds since 2000-01-01 00:00:0.0')
lat = np.asarray(daten.variables['BEST:LAT'])
lon = np.asarray(daten.variables['BEST:LON'])
ii = np.argwhere(lat > 0)

# Vertical profile of water vapour [ppmv] for first panel.
# Filter for lack of altitude data, recorded as -99999.
alt = np.asarray(daten.variables['BEST:ALT'])

fs = 12
fig = plt.figure()
#plt.subplot2grid((1,4),(0,0))
#plt.grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
plt.plot(zeit[ii[:,0]],lat[ii[:,0]],color='red')
plt.ylabel('Latitude',fontsize=fs)
plt.xlabel('Time',fontsize=fs)
plt.tick_params(labelsize=fs)
#plt.savefig('flight_track.pdf',dpi=200,bbox_inches='tight')
plt.show()
