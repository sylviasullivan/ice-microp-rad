import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import xarray as xr
import sys

#gridfi = '/work/bb1131/b380459/TROPIC/grids/icon-grid_tropic_55e170e5s40n_R2500m.nc'
#grid = xr.open_dataset(gridfi)
#landsea = grid.cell_sea_land_mask.values
xpfi = '/work/bb1131/b380459/TROPIC/extpar/extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap.nc' 
xp = xr.open_dataset(xpfi)
landsea = xp.FR_LAND.values   # we are still on an unstructured grid here...
print(landsea.shape)
print(np.nanmax(landsea),np.nanmean(landsea),np.nanmedian(landsea),np.nanmin(landsea))

basedir = '/work/bb1131/b380873/tropic_run2_output/unstructured_land_sea/'
tqi_fi = xr.open_dataset(basedir + 'tqi_unstructured.nc')
tqi_vals = tqi_fi.tqi.values
landpts = np.argwhere(landsea > 0.9)[:,0]
tqi_vals_land = tqi_vals[:,landpts]*1000
print(np.nanmax(tqi_vals_land),np.nanmean(tqi_vals_land),np.nanmin(tqi_vals_land))
seapts = np.argwhere(landsea < 0.1)[:,0]
tqi_vals_sea = tqi_vals[:,seapts]*1000
print(np.nanmax(tqi_vals_sea),np.nanmean(tqi_vals_sea),np.nanmin(tqi_vals_sea))
del tqi_vals


fs = 13
fig = plt.figure()
#wgts = np.ones_like(tqi_vals_land)/len(tqi_vals_land)*100
d = -3.5; u = 3.5; n = 30
h1, bin_edges = np.histogram(tqi_vals_land,bins=np.logspace(d,u,n))
h2, _ = np.histogram(tqi_vals_sea,bins=np.logspace(d,u,n))
print(h1,h2)

bin_center = np.zeros((bin_edges.shape[0]-1,))
for i,b in enumerate(bin_edges):
    if i != len(bin_edges)-1:
       bin_center[i] = (b + bin_edges[i+1])/2
    else:
       break

l1 = tqi_vals_land.shape[0]*tqi_vals_land.shape[1]
plt.step(bin_center,h1/l1*100,color='r',label='Land',linewidth=0.75)
plt.plot([np.nanmean(tqi_vals_land),np.nanmean(tqi_vals_land)],[0,6.5],color='r',linewidth=1.25)
l2 = tqi_vals_sea.shape[0]*tqi_vals_sea.shape[1]
plt.step(bin_center,h2/l2*100,color='b',label='Sea',linewidth=0.75)
plt.plot([np.nanmean(tqi_vals_sea),np.nanmean(tqi_vals_sea)],[0,6.5],color='b',linewidth=1.25)
plt.ylim([0,7]); plt.ylabel('Probability [%]',fontsize=fs)
#plt.xlim([0.001,3000]); 
#plt.xlim([0.1,1000]); 
plt.legend(loc='upper right')
plt.xlabel(r'IWP [kg m$^{-2}$]',fontsize=fs)
plt.gca().set_xscale('log')
#fig.savefig('IWPdist-log-land_ocean.pdf',bbox_inches='tight')
plt.show()
