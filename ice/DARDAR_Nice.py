import sys
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import cm

basedir = '/work/bb1131/b380873/tropic_vis/obs/DARDAR/'
globe = xr.open_dataset(basedir + 'globalmap.nc')

def resize_colorbar(event):
    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

ta = globe.ta
# Factor of 1000 to convert from L-1 to m-3.
# Then Sara shows 10**(8) m-2 in Fig S1.
ICNC = globe.icnc.sum('ta')/10**(5)
lat = globe.lat
lon = globe.lon

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(9,7),
              subplot_kw={'projection':ccrs.PlateCarree()})
levs = np.linspace(0,3000,200)/10**(5)
im = ax.contourf(lon,lat,ICNC,levs,cmap=cm.viridis,transform=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,color='gray',alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax.coastlines()

cbar_ax = fig.add_axes([0,0,0.1,0.1])
fig.canvas.mpl_connect('resize_event',resize_colorbar)
plt.colorbar(im,label=r'10$^{8}$ m$^{-3}$',cax=cbar_ax)
resize_colorbar(None)

fig.savefig('DARDAR_Nice_sum_over_ta.pdf')
plt.show()
