import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from cartopy import config
from matplotlib import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

basedir = '/scratch/b/b380873/tropic_run2_restart/output_2017080800-2017080806/'
step66 = basedir + 'CLCONV_3D_icon_tropic_0066_remapdis_global0.025.nc'
step66_nc = xr.open_dataset(step66)
lons = step66_nc.lon
lats = step66_nc.lat
qi = step66_nc.qi
qi_ml = qi[0,36,:,:]  # take the first time step and the nth model level

fig, ax = plt.subplots(1,1,figsize=(6,5),
                       subplot_kw={'projection': ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])


fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)
levs = np.logspace(-2,-0.5,10)
im = ax.contourf(lons,lats,qi_ml*1000,levels=levs,cmap=cm.viridis,\
                  transform=ccrs.PlateCarree())

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()

    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

fig.canvas.mpl_connect('resize_event', resize_colobar)
plt.colorbar(im,label=r'g kg$^{-1}$',cax=cbar_ax)

gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                  linewidth=1,color='gray',alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

ax.set_extent([55,170,-5,38],crs=ccrs.PlateCarree())
im.set_clim([0,0.2])

ax.coastlines()
resize_colobar(None)
#plt.savefig('qi_field.pdf',dpi=200,bbox_inches='tight')
plt.show()
