import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import sys
from cartopy import config
from matplotlib import cm,colorbar,colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

fs = 12
basedir = '/scratch/b/b380873/tropic_run2_restart/output_2017080800-2017080806/'
fi1 = basedir + 'CLCONV_3D_icon_tropic_0066_remapdis_global0.025.nc'
fi2 = basedir + 'CLCONV_3D_icon_tropic_0066_remapcon_global0.025.nc'
fi1_nc = xr.open_dataset(fi1); fi2_nc = xr.open_dataset(fi2)
lons1 = fi1_nc.lon; lats1 = fi1_nc.lat
lons2 = fi2_nc.lon; lats2 = fi2_nc.lat
qi1 = fi1_nc.qi; qi2 = fi2_nc.qi
qi1_ml = qi1[0,35,:,:]  # take the first time step and the nth model level
qi2_ml = qi2[0,35,:,:]

fig, ax = plt.subplots(2,2,figsize=(12,9),
                       subplot_kw={'projection': ccrs.PlateCarree()})
print(ax.shape)
cbar_ax = fig.add_axes([0.5, 0.5, 0.1, 0.1])  # upper half

fig.subplots_adjust(hspace=0, wspace=0.1, top=0.925, left=0.1)
im1 = ax[0,0].contourf(lons1,lats1,qi1_ml*1000,30,cmap=cm.Purples,\
                  transform=ccrs.PlateCarree())

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()

    posn = ax[0,1].get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

fig.canvas.mpl_connect('resize_event', resize_colobar)
plt.colorbar(im1,label=r'g kg$^{-1}$',cax=cbar_ax)

gl = ax[0,0].gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                  linewidth=1,color='gray',alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

ax[0,0].set_extent([55,170,-5,38],crs=ccrs.PlateCarree())
im1.set_clim([0,0.2])
ax[0,0].coastlines()
ax[0,0].set_title('(a) regriddis',fontsize=fs)

im2 = ax[0,1].contourf(lons2,lats2,qi2_ml*1000,30,cmap=cm.Purples,\
                  transform=ccrs.PlateCarree())
gl = ax[0,1].gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                  linewidth=1,color='gray',alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
im2.set_clim([0,0.2])
ax[0,1].coastlines()
ax[0,1].set_title('(b) regridcon',fontsize=fs)

im3 = ax[1,0].contourf(lons2,lats2,(qi1_ml-qi2_ml)*1000,30,cmap=cm.coolwarm,\
                  transform=ccrs.PlateCarree(),vmin=-0.025,vmax=0.025)
gl = ax[1,0].gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                  linewidth=1,color='gray',alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
im3.set_clim([-0.025,0.025])
ax[1,0].set_title(r'(c) $\Delta$regriddes-regridcon',fontsize=fs)

cbar_ax2 = fig.add_axes([0.5, 0.1, 0.04, 0.25])
cbar = colorbar.ColorbarBase(cbar_ax2,cmap=plt.get_cmap('coolwarm'),
       norm=colors.Normalize(vmin=-0.025,vmax=0.025))
cbar.set_clim(-0.025,0.025)
#def resize_colobar(event):
#    plt.draw()
#
#    posn = ax[1,0].get_position()
#    cbar_ax2.set_position([posn.x0 + posn.width + 0.01, posn.y0,
#                          0.04, posn.height])

fig.canvas.mpl_connect('resize_event', resize_colobar)
#plt.colorbar(im3,label=r'g kg$^{-1}$',cax=cbar_ax2)

plt.savefig('grid_comp.pdf',bbox_inches='tight')
plt.show()
