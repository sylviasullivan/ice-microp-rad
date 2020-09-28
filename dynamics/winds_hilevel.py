import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sys
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import cm

basedir = '/work/bb1131/b380873/tropic_run2_output/'
fi = basedir + 'UV_3D_icon_tropic_0051_remapdis_global2.5.nc'
uv51 = xr.open_dataset(fi)
lat = uv51.lat
lon = uv51.lon
u = uv51.u
v= uv51.v

l = [75, 60, 30]
windMAG = np.sqrt(u**2 + v**2)  # gives a field with dims (1, 90, 1800, 4600)
for i in np.arange(3):
    print(np.nanmin(windMAG[0,l[i]]),np.nanmean(windMAG[0,l[i]]),np.nanmax(windMAG[0,l[i]]))
windDIR = np.arctan2(v,u)

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
#def resize_colobar(event):
#    plt.draw()
#    posn = ax[0].get_position()
#    cbar_ax[0].set_position([posn.x0 + posn.width + 0.01, posn.y0,
#                          0.04, posn.height])
#    posn = ax[1].get_position()
#    cbar_ax[1].set_position([posn.x0 + posn.width + 0.01, posn.y0,
#                          0.04, posn.height])
#    posn = ax[2].get_position()
#    cbar_ax[2].set_position([posn.x0 + posn.width + 0.01, posn.y0,
#                          0.04, posn.height])


fig, ax = plt.subplots(1,1,figsize=(11,4),subplot_kw={'projection':ccrs.PlateCarree()})
#cbar_ax = [fig.add_axes([0, 0, 0.7, 0.7]), fig.add_axes([0, 0, 0.4, 0.4]), \
#           fig.add_axes([0, 0, 0.1, 0.1])]

fig.subplots_adjust(hspace=0.1, wspace=0, top=0.925, left=0.1)
lons, lats = np.meshgrid(lon, lat)
clevs = [np.linspace(0,30,40),np.linspace(0,35,40),np.linspace(0,40,40)]

i = 2
im = ax.contourf(lon, lat, windMAG[0,l[i]], clevs[i], transform=ccrs.PlateCarree(),\
               cmap=cm.jet)

# Get the coordinate reference system used by the data
ax.quiver(lons,lats,u[0,l[i]], v[0,l[i]], transform=ccrs.PlateCarree())

gl = ax.gridlines()
gl.xlabels_bottom = True
gl.ylabels_left = True
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#    fig.canvas.mpl_connect('resize_event', resize_colobar)
plt.colorbar(im,label=r'm s$^{-1}$') #,cax=cbar_ax[i])

ax.set_extent([55,170,-5,40],crs=ccrs.PlateCarree())
ax.coastlines()
fig.savefig('wind_hilevel.pdf',bbox_inches='tight')
plt.show()

