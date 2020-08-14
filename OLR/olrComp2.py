import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys
from datetime import datetime
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Read in the CERES data.
basedir = '/work/bb1131/b380873/tropic_vis/obs/CERES/'
olr_file = basedir + 'CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20170801-20170831_full-domain.nc'
olr_data = xr.open_dataset(olr_file)
olr_ceres = np.abs(olr_data.toa_lw_all_1h.values)
lons_ceres = olr_data.lon
lats_ceres = olr_data.lat
zeit_ceres = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub = olr_ceres[ii]   # shape = 35, 115
zeit_ceres_sub = zeit_ceres[ii]

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub1 = olr_ceres[ii]   # shape = 35, 115
zeit_ceres_sub1 = zeit_ceres[ii]

# Read in the ERA5 values.
basedir = '/work/bb1131/b380873/tropic_vis/obs/ERA5/'
olr_file = basedir + 'ERA5_OLR_1deg[55-170]-20170805-20170809.nc'
olr_data = xr.open_dataset(olr_file)
olr_era5 = np.abs(olr_data.mtnlwrf.values)
lons_era5 = olr_data.longitude
lats_era5 = olr_data.latitude
zeit_era5 = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub = olr_era5[ii]
# shape = 120, 721, 1440
zeit_era5_sub = zeit_era5[ii]

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub1 = olr_era5[ii]
# shape = 120, 721, 1440
zeit_era5_sub1 = zeit_era5[ii]

# Read in the ICON 1-mom simulation.
basedir = '/work/bb1131/b380873/tropic_run2_output/'
olr_file = basedir + 'OLR_TOA_all_1deg.nc'
olr_data = xr.open_dataset(olr_file)
olr_icon = np.abs(olr_data.lwflxall.values)
lons_icon = olr_data.lon
lats_icon = olr_data.lat
zeit_icon = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub = olr_icon[ii,0]
zeit_icon_sub = zeit_icon[ii]

exttime = datetime(2017,8,8,0,6)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub1 = olr_icon[ii,0]
zeit_icon_sub1 = zeit_icon[ii]

# Read in the ICON 2-mom simulation.
basedir = '/work/bb1131/b380873/tropic_run5_output/'
olr_file = basedir + 'OLR_120-141_1.0deg.nc'
olr_data = xr.open_dataset(olr_file)
olr_icon2 = np.abs(olr_data.thb_t.values)
lons_icon2 = olr_data.lon
lats_icon2 = olr_data.lat
zeit_icon2 = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon2 >= np.datetime64(exttime))[0,0]
olr_icon_sub2 = olr_icon2[ii]
zeit_icon_sub2 = zeit_icon2[ii]

exttime = datetime(2017,8,8,0,6)
ii = np.argwhere(zeit_icon2 >= np.datetime64(exttime))[0,0]
olr_icon_sub12 = olr_icon2[ii]
zeit_icon_sub12 = zeit_icon2[ii]

print('       CERES    ERA5     ICON-1mom    ICON-2mom')
print('Shape: ' + str(olr_ceres_sub.shape) + ' ' + str(olr_era5_sub.shape) + ' ' + \
      str(olr_icon_sub.shape) + ' ' + str(olr_icon_sub2.shape))
print('Mean: ' + str(np.nanmean(olr_ceres_sub)) + '  ' + str(np.nanmean(olr_era5_sub)) + ' ' + \
      str(np.nanmean(olr_icon_sub)) + ' ' + str(np.nanmean(olr_icon_sub2)))
print('Median: ' + str(np.nanmedian(olr_ceres_sub)) + '  ' + str(np.nanmedian(olr_era5_sub)) + \
      ' ' + str(np.nanmedian(olr_icon_sub)) + ' ' + str(np.nanmedian(olr_icon_sub2)))
print('Std: ' + str(np.nanstd(olr_ceres_sub)) + ' ' + str(np.nanstd(olr_era5_sub)) + ' ' + \
      str(np.nanstd(olr_icon_sub)) + ' ' + str(np.nanstd(olr_icon_sub2)))

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()
    posn = ax[0,1].get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

levs = np.linspace(90,375,15)
titre = ['CERES OLR: ' + str(zeit_ceres_sub),'ERA5 OLR: '+ str(zeit_era5_sub), \
         'ICON-1mom OLR: ' + str(zeit_icon_sub),'ICON-2mom OLR: ' + str(zeit_icon_sub2), \
         'CERES OLR: ' + str(zeit_ceres_sub1),'ERA5 OLR: ' + str(zeit_era5_sub1), \
         'ICON-1mom OLR: ' + str(zeit_icon_sub1),'ICON-2mom OLR: ' + str(zeit_icon_sub12)]
lons = [lons_ceres, lons_era5, lons_icon, lons_icon2, lons_ceres, lons_era5, lons_icon, lons_icon2]
lats = [lats_ceres, lats_era5, lats_icon, lats_icon2, lats_ceres, lats_era5, lats_icon, lats_icon2]
olr = [olr_ceres_sub, olr_era5_sub, olr_icon_sub, olr_icon_sub2, olr_ceres_sub1, olr_era5_sub1, \
       olr_icon_sub1, olr_icon_sub12]

fig, ax = plt.subplots(nrows=4,ncols=2,figsize=(11,11.1),subplot_kw={'projection':ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0.1, wspace=0, top=0.925, left=0.1)

c = 0
for j in np.arange(2):
    for i in np.arange(4):
        print(i,j,c)
        ax[i,j].set_title(titre[c])
        im = ax[i,j].contourf(lons[c],lats[c],olr[c],levels=levs,cmap=cm.viridis,\
                  transform=ccrs.PlateCarree())
        gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                linewidth=1,color='gray',alpha=0.5,linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabels_bottom = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        ax[i,j].set_extent([55,115,-5,38],crs=ccrs.PlateCarree())
        im.set_clim([90,375])
        ax[i,j].coastlines()
        c += 1

fig.canvas.mpl_connect('resize_event', resize_colobar)
plt.colorbar(im,label=r'W m$^{-2}$',cax=cbar_ax)

resize_colobar(None)
#fig.savefig('olr-comparison_115e_1deg.pdf',bbox_inches='tight')
plt.show()
