import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import cm,colors
import xarray as xr
import sys
from datetime import datetime
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

execfile('../MidpointNormalize.py')

# Read in the ERA5 values.
basedir = '/work/bb1018/b380873/tropic_vis/obs/ERA5/'
w_file = basedir + 'ERA5_w_0.25deg-20170807-20170808.nc'
w_data = xr.open_dataset(w_file)
w250_era5 = w_data.w.values[:,0]   # The first level is 250 hPa. The second level is 500 hPa.
w500_era5 = w_data.w.values[:,1]
lons_era5 = w_data.longitude
lats_era5 = w_data.latitude
zeit_era5 = w_data.time.values

# Extract 8 August 2017 at 1 am
exttime = datetime(2017,8,8,1,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
w250_era5_sub = w250_era5[ii]
w500_era5_sub = w500_era5[ii]
zeit_era5_sub = zeit_era5[ii]
print(zeit_era5_sub)

# Read in the ICON-1mom simulation at 8 August 2017 1 am.
basedir = '/work/bb1018/b380873/tropic_run2_output/'
w_file = basedir + 'W_3D_icon_tropic_0061_remapdis_global0.025.nc'
#w_file = basedir + 'w250_ALL-0051-0061.nc'
w_data = xr.open_dataset(w_file)
w250_icon_sub = w_data.w.values[0,53]
# Index 64 is roughly 500 hPa. Index 53 is roughly 250 hPa.
#w_file = basedir + 'w500_ALl-0051-0061.nc'
w500_icon_sub = w_data.w.values[0,64]
lons_icon = w_data.lon
lats_icon = w_data.lat
zeit_icon_sub = w_data.time.values[0]
print(zeit_icon_sub)

# Read in the ICON-2mom simulation at 8 August 2017 1 am.
basedir = '/work/bb1018/b380873/tropic_run5_output/'
#basedir = '/scratch/b/b380873/tropic_run5/'
w_file = basedir + 'OMEGA_60-63_0.025deg_PL.nc'
#w_file = basedir + 'W_icon_tropic_0061remapdis_0.025_PL.nc'
w_data = xr.open_dataset(w_file) # Pressure levels correspond to 900,850,750,500,250 hPa
w250_icon_sub2 = w_data.omega.isel(time=1,plev=4) # plev_2
w500_icon_sub2 = w_data.omega.isel(time=1,plev=3)
print(w500_icon_sub2.time)
lons_icon2 = w_data.lon
lats_icon2 = w_data.lat
zeit_icon_sub2 = w_data.time.values[1]
print(zeit_icon_sub)

print('ERA5     ICON-1mom     ICON-2mom    (w250)')
print('Shape: ' + str(w250_era5_sub.shape) + ' ' + str(w250_icon_sub.shape) + ' ' + str(w250_icon_sub2.shape))
print('Max: ' + str(np.nanmax(w250_era5_sub)) + ' ' + str(np.nanmax(w250_icon_sub)) + ' ' + str(np.nanmax(w250_icon_sub2)))
print('Mean: ' + str(np.nanmean(w250_era5_sub)) + ' ' + str(np.nanmean(w250_icon_sub)) + ' ' + str(np.nanmean(w250_icon_sub2)))
print('Std: ' + str(np.nanstd(w250_era5_sub)) + ' ' + str(np.nanstd(w250_icon_sub)) + ' ' + str(np.nanstd(w250_icon_sub2)))
print('Median: ' + str(np.nanmedian(w250_era5_sub)) + ' ' + str(np.nanmedian(w250_icon_sub)) + ' ' + str(np.nanmedian(w250_icon_sub2)))
print('Min: ' + str(np.nanmin(w250_era5_sub)) + ' ' + str(np.nanmin(w250_icon_sub)) + ' ' + str(np.nanmin(w250_icon_sub2)))
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colorbar(event):
    plt.draw()
    posn = ax[0,1].get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

levs = np.linspace(-2,3,30)
titre = [r'ERA5 $\omega_{250}$: '+ str(zeit_era5_sub)[:16], r'ICON-2mom $\omega_{250}$: '+ str(zeit_icon_sub2)[:16],\
         r'ERA5 $\omega_{500}$: ' + str(zeit_era5_sub)[:16], r'ICON-2mom $\omega_{500}$: '+ str(zeit_icon_sub2)[:16]]
        #r'ICON-1mom $\omega_{250}$: ' + str(zeit_icon_sub)[:16], r'ICON-1mom $\omega_{250}$: ' + str(zeit_icon_sub)[:16],
lons = [lons_era5, lons_icon2, lons_era5, lons_icon2] # lons_icon, lons_icon
lats = [lats_era5, lats_icon2, lats_era5, lats_icon2] # lats_icon, lats_icon
w = [w250_era5_sub, w250_icon_sub2, w500_era5_sub, w500_icon_sub2] # w250_icon_sub, w500_icon_sub

fig, ax = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(11,11))
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0.1, wspace=0.05, top=0.925, left=0.1)

c = 0
for j in np.arange(2):
    for i in np.arange(2):
        ax[i,j].set_title(titre[c])
        im = ax[i,j].contourf(lons[c],lats[c],w[c],cmap=cm.seismic,transform=ccrs.PlateCarree(),levels=levs,
              norm=MidpointNormalize(midpoint=0.,vmin=-3,vmax=3.0))
        gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                   linewidth=1,color='gray',alpha=0.5,linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        if i < 2:
           gl.xlabels_bottom = False
        if j == 1:
           gl.ylabels_left = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        ax[i,j].set_extent([55,115,-5,38],crs=ccrs.PlateCarree())
        #im.set_clim([-1,2])
        ax[i,j].coastlines()
        #plt.colorbar(im,ax=ax[i,j])
        c = c + 1

fig.tight_layout()
fig.canvas.mpl_connect('resize_event', resize_colorbar)
plt.colorbar(im,label=r'm s$^{-1}$',cax=cbar_ax)

resize_colorbar(None)
#fig.savefig('w-comparison_115e_0.25-0.025deg.png',bbox_inches='tight')
plt.show()
