import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys
from datetime import datetime
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Read in the CERES data.
basedir = '/work/bb1018/b380873/tropic_vis/obs/CERES/'
olr_data = xr.open_dataset(basedir + 'CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20170801-20170831_full-domain.nc')
olr_ceres = olr_data.toa_lw_all_1h.sel(time=slice(datetime(2017,8,7,12,0,0),datetime(2017,8,8,12,0,0)))
olr_ceres = olr_ceres.mean(dim={'time'})
lons_ceres = olr_data.lon
lats_ceres = olr_data.lat

# Read in the ERA5 values.
basedir = '/work/bb1018/b380873/tropic_vis/obs/ERA5/'
olr_data = xr.open_dataset(basedir + 'ERA5_OLR-20170805-20170809.nc')
olr_era5 = olr_data.ttr.sel(time=slice(datetime(2017,8,7,12,0,0),datetime(2017,8,8,12,0,0)))
olr_era5 = -1.*olr_era5.mean(dim={'time'})/3600
lons_era5 = olr_data.longitude
lats_era5 = olr_data.latitude

# Read in the ICON 1-mom simulation.
basedir = '/work/bb1018/b380873/model_output/ICON/'
olr_data = xr.open_dataset(basedir + 'OLR_24h_tropic_run2_ll.nc')
olr_icon_1mom = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)
lons_icon = olr_data.lon
lats_icon = olr_data.lat

# Read in the ICON 2-mom simulation.
olr_data = xr.open_dataset(basedir + 'OLR_24h_tropic_run5_ll.nc')
olr_icon_2mom = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)

# Read in the ICON no2mom simulation.
olr_data = xr.open_dataset(basedir + 'OLR_12h_tropic_run5_no2mom_ll.nc')
olr_icon_no2mom = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)

# Read in the ICON novgrid simulation.
olr_data = xr.open_dataset(basedir + 'OLR_12h_tropic_run5_novgrid_ll.nc')
olr_icon_novgrid = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)

# Read in the ICON rad2mom simulation.
olr_data = xr.open_dataset(basedir + 'OLR_24h_tropic_run7_rad2mom_ll.nc')
olr_icon_rad2mom = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)

# Read in the ICON PDA simulation.
olr_data = xr.open_dataset(basedir + 'OLR_24h_tropic_run8_pda_ll.nc')
olr_icon_pda = -1.*olr_data.lwflxall.mean(dim={'time'}).sel(height=1)

print('CERES   ERA5   1mom   2mom   no2mom   novgrid   rad2mom   pda')
#print('Shape: ' + str(olr_ceres_sub.shape) + ' ' + str(olr_era5_sub.shape) + ' ' + \
#      str(olr_icon_sub.shape) + ' ' + str(olr_icon_sub2.shape))
print('Mean: ' + str(np.nanmean(olr_ceres)) + '  ' + str(np.nanmean(olr_era5)) + ' ' + \
      str(np.nanmean(olr_icon_1mom)) + ' ' + str(np.nanmean(olr_icon_2mom)) + ' ' + \
      str(np.nanmean(olr_icon_no2mom)) + ' ' + str(np.nanmean(olr_icon_novgrid)) + ' ' + \
      str(np.nanmean(olr_icon_rad2mom)) + ' ' + str(np.nanmean(olr_icon_pda)))
print('Median: ' + str(np.nanmedian(olr_ceres)) + '  ' + str(np.nanmedian(olr_era5)) + ' ' + \
      str(np.nanmedian(olr_icon_1mom)) + ' ' + str(np.nanmedian(olr_icon_2mom)) + ' ' + \
      str(np.nanmedian(olr_icon_no2mom)) + ' ' + str(np.nanmedian(olr_icon_novgrid)) + ' ' + \
      str(np.nanmedian(olr_icon_rad2mom)) + ' ' + str(np.nanmedian(olr_icon_pda)))
print('Std: ' + str(np.nanstd(olr_ceres)) + ' ' + str(np.nanstd(olr_era5)) + ' ' + \
      str(np.nanstd(olr_icon_1mom)) + ' ' + str(np.nanstd(olr_icon_2mom)) + ' ' + \
      str(np.nanstd(olr_icon_no2mom)) + ' ' + str(np.nanstd(olr_icon_novgrid)) + ' ' + \
      str(np.nanstd(olr_icon_rad2mom)) + ' ' + str(np.nanstd(olr_icon_pda)))

# pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colobar(event):
    plt.draw()
    posn = ax[0,1].get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])

levs = np.linspace(80,375,15)
titre = ['24h-mean CERES OLR', '24h-mean ERA5 OLR', '24h-mean ICON 1mom', '24h-mean ICON 2mom',\
         '24h-mean ICON no2mom', '24h-mean ICON novgrid', '24h-mean ICON rad2mom',\
         '24h-mean PDA']
farbe = ['gray', 'black', 'red', 'green', 'blue', 'gold', 'purple', 'pink']
lons = [lons_ceres, lons_era5, lons_icon, lons_icon, lons_icon, lons_icon, lons_icon, lons_icon]
lats = [lats_ceres, lats_era5, lats_icon, lats_icon, lats_icon, lats_icon, lats_icon, lats_icon]
olr = [olr_ceres, olr_era5, olr_icon_1mom, olr_icon_2mom, olr_icon_no2mom, olr_icon_novgrid,\
       olr_icon_rad2mom, olr_icon_pda]

fig, ax = plt.subplots(nrows=4,ncols=2,figsize=(9,11.1),subplot_kw={'projection':ccrs.PlateCarree()})
cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0.1, wspace=0, top=0.925, left=0.1)

c = 0
for j in np.arange(2):
    for i in np.arange(4):
        print(i,j,c)
        ax[i,j].set_title(titre[c],color=farbe[c])
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
#fig.savefig('../output/olr-comparison_24h.png',bbox_inches='tight')
plt.show()

fig2, ax2 = plt.subplots(nrows=1,ncols=3,figsize=(6.5,6.5))
fs = 16
lw = 2.5
ax2[0].plot([1,2],[np.nanmean(olr_ceres),np.nanmean(olr_ceres)],linewidth=lw,color='gray')
ax2[0].plot([1,2],[np.nanmean(olr_era5),np.nanmean(olr_era5)],linewidth=lw,color='black')
ax2[0].plot([1,2],[np.nanmean(olr_icon_1mom),np.nanmean(olr_icon_1mom)],linewidth=lw,color='red')
ax2[0].plot([1,2],[np.nanmean(olr_icon_no2mom),np.nanmean(olr_icon_no2mom)],linewidth=lw,linestyle='--',color='red')
ax2[0].plot([1,2],[np.nanmean(olr_icon_2mom),np.nanmean(olr_icon_2mom)],linewidth=lw,color='blue')
#ax2[0].plot([1,2],[np.nanmean(olr_icon_no2mom),np.nanmean(olr_icon_no2mom)],linewidth=lw,color='blue')
#ax2[0].plot([1,2],[np.nanmean(olr_icon_novgrid),np.nanmean(olr_icon_novgrid)],linewidth=lw,color='gold')
ax2[0].plot([1,2],[np.nanmean(olr_icon_rad2mom),np.nanmean(olr_icon_rad2mom)],linewidth=lw,color='goldenrod')
ax2[0].plot([1,2],[np.nanmean(olr_icon_pda),np.nanmean(olr_icon_pda)],linewidth=lw,color='green')
ax2[0].set_xlim([0,4])
ax2[0].set_ylim([200,240])
ax2[0].tick_params('both',labelsize=fs,rotation=45)
ax2[0].spines["top"].set_visible(False)
ax2[0].spines["bottom"].set_visible(False)
ax2[0].spines["right"].set_visible(False)
ax2[0].set_ylabel(r'Mean OLR [W m$^{-2}$]',fontsize=fs) #Spatial mean OIR
ax2[0].get_xaxis().set_visible(False)

ax2[1].plot([1,2],[np.nanmedian(olr_ceres),np.nanmedian(olr_ceres)],linewidth=lw,color='gray')
ax2[1].plot([1,2],[np.nanmedian(olr_era5),np.nanmedian(olr_era5)],linewidth=lw,color='black')
print(np.nanmedian(olr_ceres),np.nanmedian(olr_era5))
ax2[1].plot([1,2],[np.nanmedian(olr_icon_1mom),np.nanmedian(olr_icon_1mom)],linewidth=lw,color='red')
ax2[1].plot([1,2],[np.nanmedian(olr_icon_2mom),np.nanmedian(olr_icon_2mom)],linewidth=lw,color='green')
ax2[1].plot([1,2],[np.nanmedian(olr_icon_no2mom),np.nanmedian(olr_icon_no2mom)],linewidth=lw,color='blue')
ax2[1].plot([1,2],[np.nanmedian(olr_icon_novgrid),np.nanmedian(olr_icon_novgrid)],linewidth=lw,color='gold')
ax2[1].plot([1,2],[np.nanmedian(olr_icon_rad2mom),np.nanmedian(olr_icon_rad2mom)],linewidth=lw,color='purple')
ax2[1].plot([1,2],[np.nanmedian(olr_icon_pda),np.nanmedian(olr_icon_pda)],linewidth=lw,color='pink')
ax2[1].set_xlim([0,4])
ax2[1].set_ylim([200,240])
ax2[1].tick_params('both',labelsize=fs,rotation=45)
ax2[1].spines["top"].set_visible(False)
ax2[1].spines["bottom"].set_visible(False)
ax2[1].spines["right"].set_visible(False)
ax2[1].set_ylabel(r'Median OLR [W m$^{-2}$]',fontsize=fs)  #Median OIR
ax2[1].get_xaxis().set_visible(False)

ax2[2].plot([1,2],[np.nanstd(olr_ceres),np.nanstd(olr_ceres)],linewidth=lw,color='gray')
ax2[2].plot([1,2],[np.nanstd(olr_era5),np.nanstd(olr_era5)],linewidth=lw,color='black')
ax2[2].plot([1,2],[np.nanstd(olr_icon_1mom),np.nanstd(olr_icon_1mom)],linewidth=lw,color='red')
ax2[2].plot([1,2],[np.nanstd(olr_icon_2mom),np.nanstd(olr_icon_2mom)],linewidth=lw,color='green')
ax2[2].plot([1,2],[np.nanstd(olr_icon_no2mom),np.nanstd(olr_icon_no2mom)],linewidth=lw,color='blue')
ax2[2].plot([1,2],[np.nanstd(olr_icon_novgrid),np.nanstd(olr_icon_novgrid)],linewidth=lw,color='gold')
ax2[2].plot([1,2],[np.nanstd(olr_icon_rad2mom),np.nanstd(olr_icon_rad2mom)],linewidth=lw,color='purple')
ax2[2].plot([1,2],[np.nanstd(olr_icon_pda),np.nanstd(olr_icon_pda)],linewidth=lw,color='pink')
ax2[2].set_xlim([0,4])
ax2[2].set_ylim([38,60])
ax2[2].tick_params('both',labelsize=fs,rotation=45)
ax2[2].spines["top"].set_visible(False)
ax2[2].spines["bottom"].set_visible(False)
ax2[2].spines["right"].set_visible(False)
ax2[2].set_ylabel(r'Stdev OLR [W m$^{-2}$]',fontsize=fs)  #Stdev OIR
ax2[2].get_xaxis().set_visible(False)

#fig2.savefig('../output/OLR_stats_all2.pdf',bbox_inches='tight')  #OIR
plt.show()
