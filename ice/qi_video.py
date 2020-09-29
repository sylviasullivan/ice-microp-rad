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

basedir = '/work/bb1131/b380873/tropic_run2_output/'
qi_file = basedir + 'qi35_all_try2.nc'
qi_data = xr.open_dataset(qi_file)
qi_vals = qi_data.qi.values
lons = qi_data.lon
lats = qi_data.lat
zeit = qi_data.time

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()

    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])


fig, ax = plt.subplots(1,1,figsize=(11.5,3.7),
                       subplot_kw={'projection': ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

levs = np.logspace(-2,-0.5,10)
def animate(i):
    print(i)
    ax.clear()
    print(np.nanmin(qi_vals[i,0]*1000),np.nanmean(qi_vals[i,0]*1000),np.nanmax(qi_vals[i,0]*1000))
    im = ax.contourf(lons,lats,qi_vals[i,0]*1000,levels=levs,cmap=cm.viridis,\
                  transform=ccrs.PlateCarree()) 

    fig.canvas.mpl_connect('resize_event', resize_colobar)
    plt.colorbar(im,label=r'g kg$^{-1}$',cax=cbar_ax)

    ax.set_title(r'$q_i$: '+ str(zeit[i].values))
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


ananas = ani.FuncAnimation(fig,animate,35,interval=700,blit=False)
#ananas.save('qi_full2.mp4')
plt.show()
