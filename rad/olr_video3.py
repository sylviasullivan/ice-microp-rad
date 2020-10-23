import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import glob
import xarray as xr
import sys,os,subprocess
import numpy as np
from cartopy import config
from matplotlib import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

basedir = '/work/bb1018/b380873/tropic_run2_output/'
olr_file = basedir + 'OLR_TOA_all.nc'
olr_data = xr.open_dataset(olr_file)
olr_vals = np.abs(olr_data.lwflxall.values)
lons = olr_data.lon
lats = olr_data.lat
zeit = olr_data.time
print(zeit[0])
print(zeit[-1])
sys.exit()
# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()

    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])


fig, ax = plt.subplots(1,1,figsize=(11,3.7),
                       subplot_kw={'projection': ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

#levs = np.logspace(-1,2.5,15)
#levs = np.linspace(0,800,15)
levs = np.linspace(100,390,15)
def animate(i):
    print(i)
    ax.clear()
    print(olr_vals[i,0].shape)
    print(np.nanmin(olr_vals[i,0]),np.nanmean(olr_vals[i,0]),np.nanmax(olr_vals[i,0]))
    im = ax.contourf(lons,lats,olr_vals[i,0],levels=levs,cmap=cm.viridis,\
                  transform=ccrs.PlateCarree()) 

    fig.canvas.mpl_connect('resize_event', resize_colobar)
    plt.colorbar(im,label=r'W m$^{-2}$',cax=cbar_ax)

    ax.set_title('TOA OLR: '+ str(zeit[i].values))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                linewidth=1,color='gray',alpha=0.5,linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.set_extent([55,170,-5,38],crs=ccrs.PlateCarree())
    im.set_clim([100,400])
    ax.coastlines()
    resize_colobar(None)


ananas = ani.FuncAnimation(fig,animate,35,interval=700,blit=False)
ananas.save('olr_full.mp4')
plt.show()
