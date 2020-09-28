import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys
from datetime import datetime
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# The CERES and ERA5 block can be read in and arranged from olrComp2.py
# Read in the ICON no2mom simulation.
basedir = '/work/bb1131/b380873/tropic_run5_output/no2mom/'
olr_file = basedir + 'OLR_120-140_0.025deg.nc'
olr_data = xr.open_dataset(olr_file)
olr_no2mom = np.abs(olr_data.lwflxall.values)
lons_no2mom = olr_data.lon
lats_no2mom = olr_data.lat
zeit_no2mom = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_no2mom >= np.datetime64(exttime))[0,0]
olr_no2mom_sub = olr_no2mom[ii,0]
zeit_no2mom_sub = zeit_no2mom[ii]

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_no2mom >= np.datetime64(exttime))[0,0]
olr_no2mom_sub1 = olr_no2mom[ii,0]
zeit_no2mom_sub1 = zeit_no2mom[ii]

# Read in the ICON novgrid simulation.
basedir = '/work/bb1131/b380873/tropic_run5_output/novgrid/'
olr_file = basedir + 'OLR_120-140_0.025deg.nc'
olr_data = xr.open_dataset(olr_file)
olr_novgrid = np.abs(olr_data.lwflxall.values)
lons_novgrid = olr_data.lon
lats_novgrid = olr_data.lat
zeit_novgrid = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_novgrid >= np.datetime64(exttime))[0,0]
olr_novgrid_sub = olr_novgrid[ii,0]
zeit_novgrid_sub = zeit_novgrid[ii]

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_novgrid >= np.datetime64(exttime))[0,0]
olr_novgrid_sub1 = olr_novgrid[ii,0]
zeit_novgrid_sub1 = zeit_novgrid[ii]

# Read in the ICON 1-mom simulation.
basedir = '/work/bb1131/b380873/tropic_run2_output/'
olr_file = basedir + 'OLR_TOA_all.nc' #'OLR_TOA_all_1deg.nc'
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

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub1 = olr_icon[ii,0]
zeit_icon_sub1 = zeit_icon[ii]

# Read in the ICON 2-mom simulation.
basedir = '/work/bb1131/b380873/tropic_run5_output/'
olr_file = basedir + 'OLR_120-141_0.025deg.nc' # 'OLR_120-141_1.0deg.nc'
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

exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon2 >= np.datetime64(exttime))[0,0]
olr_icon_sub12 = olr_icon2[ii]
zeit_icon_sub12 = zeit_icon2[ii]

print('       no2mom    novgrid     ICON-1mom    ICON-2mom')
print('Shape: ' + str(olr_no2mom_sub.shape) + ' ' + str(olr_novgrid_sub.shape) + ' ' + \
      str(olr_icon_sub.shape) + ' ' + str(olr_icon_sub2.shape))
print('Mean: ' + str(np.nanmean(olr_no2mom_sub)) + '  ' + str(np.nanmean(olr_novgrid_sub)) + ' ' + \
      str(np.nanmean(olr_icon_sub)) + ' ' + str(np.nanmean(olr_icon_sub2)))
print('Median: ' + str(np.nanmedian(olr_no2mom_sub)) + '  ' + str(np.nanmedian(olr_novgrid_sub)) + \
      ' ' + str(np.nanmedian(olr_icon_sub)) + ' ' + str(np.nanmedian(olr_icon_sub2)))
print('Std: ' + str(np.nanstd(olr_no2mom_sub)) + ' ' + str(np.nanstd(olr_novgrid_sub)) + ' ' + \
      str(np.nanstd(olr_icon_sub)) + ' ' + str(np.nanstd(olr_icon_sub2)))

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()
    posn = ax[0,1].get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

levs = np.linspace(80,375,15)
titre = ['no2mom OLR: ' + str(zeit_no2mom_sub)[:16],'novgrid OLR: '+ str(zeit_novgrid_sub)[:16], \
         'ICON-1mom OLR: ' + str(zeit_icon_sub)[:16],'ICON-2mom OLR: ' + str(zeit_icon_sub2)[:16], \
         'no2mom OLR: ' + str(zeit_no2mom_sub1)[:16],'novgrid OLR: ' + str(zeit_novgrid_sub1)[:16], \
         'ICON-1mom OLR: ' + str(zeit_icon_sub1)[:16],'ICON-2mom OLR: ' + str(zeit_icon_sub12)[:16]]
lons = [lons_no2mom, lons_novgrid, lons_icon, lons_icon2, lons_no2mom, lons_novgrid, lons_icon, lons_icon2]
lats = [lats_no2mom, lats_novgrid, lats_icon, lats_icon2, lats_no2mom, lats_novgrid, lats_icon, lats_icon2]
olr = [olr_no2mom_sub, olr_novgrid_sub, olr_icon_sub, olr_icon_sub2, olr_no2mom_sub1, olr_novgrid_sub1, \
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
        im.set_clim([80,375])
        ax[i,j].coastlines()
        c += 1

fig.canvas.mpl_connect('resize_event', resize_colobar)
plt.colorbar(im,label=r'W m$^{-2}$',cax=cbar_ax)

resize_colobar(None)
#fig.savefig('../output/olr-comparison_115e_novgrid_no2mom.png',bbox_inches='tight')
plt.show()
