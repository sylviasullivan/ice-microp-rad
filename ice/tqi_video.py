import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import glob
import xarray as xr
import sys,os,subprocess
import numpy as np
from matplotlib.colors import LogNorm
from cartopy import config
from matplotlib import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#basedir = '/work/bb1131/b380873/tropic_run2_output/'
#tqi_file = basedir + 'tqi_all.nc'
sim = '0V2M0A0R'
flighttrack = False

if flighttrack == True:
   basedir = '/scratch/b/b380873/' + sim + '/flight-track/'
   tqi_file = basedir + 'CLCONV_2D_F10MIN_icon_tropic_flight-track.nc'
   tag = 'flight-track'
   figsz = (4,4.5)
else:
   basedir = '/scratch/b/b380873/' + sim + '/'
   tqi_file = basedir + 'CLCONV_2D_F10MIN_icon_tropic.nc' #'TQI_
   tag = 'full'
   figsz = (7.5,4.5)

tqi_data = xr.open_dataset(tqi_file)
tqi_vals = tqi_data.tqi.values
lons = tqi_data.lon
lats = tqi_data.lat
zeit = tqi_data.time

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()

    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])


fig, ax = plt.subplots(1,1,figsize=figsz,
                       subplot_kw={'projection': ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

#levs = np.logspace(-1,2.5,15)
#levs = np.linspace(1,500,15)
levs = np.logspace(0,3,15)
def animate(i):
    print(i)
    ax.clear()
    print(np.nanmin(tqi_vals[i]*10**3),np.nanmean(tqi_vals[i]*10**3),np.nanmax(tqi_vals[i]*10**3))
    im = ax.contourf(lons,lats,tqi_vals[i]*10**3,levels=levs,cmap=cm.viridis,\
                  transform=ccrs.PlateCarree(),norm=LogNorm())

    fig.canvas.mpl_connect('resize_event', resize_colobar)
    plt.colorbar(im,label=r'g m$^{-2}$',cax=cbar_ax,ticks=[1,10,100,1000])

    ax.set_title('IWP: '+ str(zeit[i].values))
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                linewidth=1,color='gray',alpha=0.5,linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'rotation': 30}

    if flighttrack == True:
       ax.set_extent([83,89,18,30],crs=ccrs.PlateCarree())
    else:
       ax.set_extent([55,115,-5,38],crs=ccrs.PlateCarree())
    im.set_clim([1,1000])
    ax.coastlines()
    resize_colobar(None)


ananas = ani.FuncAnimation(fig,animate,18,interval=700,blit=False)
ananas.save('../output/TQI_F10MIN_' + tag + '_' + sim + '.mp4')
plt.show()
