# run in conda environment pyhdf
# /work/bb1018/b380459/TROPIC/extpar/extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap.nc
# as input
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys
sys.path.insert(1,'/work/bb1018/b380873/tropic_vis/CloudSat/')
import xarray as xr
import matplotlib.colors as colors
from CloudSat_read import read_cloudsatcalipso_hdf_file

file = sys.argv[1:]
dset = xr.open_dataset(str(file[0]))
lon_min = np.min(dset['lon']).item()
lon_max = np.max(dset['lon']).item()
lat_min = np.min(dset['lat']).item()
lat_max = np.max(dset['lat']).item()
#print('Longitude min. %7.4f, max. %7.4f' % (lon_min, lon_max) )
#print('Latitude min. %7.4f, max. %7.4f' % (lat_min, lat_max) )

def get_projection(dset,axes):
    x, y = dset['lon'].values, dset['lat'].values
    m = Basemap(projection='mill', llcrnrlon=dset['lon'].min()-1.,
             llcrnrlat=dset['lat'].min()-1.,urcrnrlon=dset['lon'].max()+1.,
             urcrnrlat=dset['lat'].max()+1.,resolution='i',ax=axes)
    m.drawcountries(linewidth=0.5, color='black')
    m.shadedrelief()
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.0, 90.0, 5.), linewidth=0.5, color='white',
              labels=[True, False, False, True], fontsize=14)
    m.drawmeridians(np.arange(0.0, 360.0, 5.), linewidth=0.5, color='white',
              labels=[True, False, False, True], fontsize=14)
    x, y = m(x, y)
    return(m, x, y)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    """Truncate a coloramp by specifying the start and endpoint."""
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n}.{a:.2f},{b:.2f})'.format(n=cmap.name,a=minval,b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return(new_cmap)

# Read in the flight track from Martina Kraemer's data
basedir = '/work/bb1018/b380873/tropic_vis/'
scfi = basedir + 'obs/stratoclim2017.geophysika.0808_1.master.ci_eval.nc'
sc_data = xr.open_dataset(scfi)
lat_sc = sc_data['BEST:LAT'].values
lon_sc = sc_data['BEST:LON'].values
t_sc = sc_data['time'].values
i_sc = np.argwhere((~np.isnan(lat_sc)) & (~np.isnan(lon_sc)) & (lat_sc > 0) & (lon_sc > 0))

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,10))
m, x, y = get_projection(dset,axes=ax)
cmap = truncate_colormap(plt.get_cmap('terrain'), 0.2, 1.)
m.contourf(x, y, dset['topography_c'], tri=True, cmap=cmap,
     levels=np.arange(-100.,3000.,10.))

#lon_center = np.rad2deg(dset.clon.values).mean()
#lat_center = np.rad2deg(dset.clat.values).mean()
#x, y = np.rad2deg(dset.clon.values), np.rad2deg(dset.clat.values)
#ax = plt.axes(projection=ccrs.PlateCarree())   # false_easting, false_northing
#ax.add_feature(cfeature.COASTLINE,lw=0.5)
#ax.add_feature(cfeature.BORDERS,lw=0.5)
#ax.add_patch(mpatches.Rectangle(xy=[55,-5],width=115,height=45,color='red',alpha=0.4,
#             transform=ccrs.PlateCarree()))


# Pulling from https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
# Create a set of line segments so that we can color them individually
xx,yy = m(lon_sc[i_sc[:,0]],lat_sc[i_sc[:,0]])
points = np.array([xx,yy]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]],axis=1)

# Convert the times from np.datetime64 to float
t_sc = t_sc[i_sc[:,0]]
t_sc_f = t_sc.astype("float")/1000000000.0
t_sc_f = t_sc_f - np.nanmin(t_sc_f)
norm = plt.Normalize(t_sc_f.min(),t_sc_f.max())
lc = LineCollection(segments,cmap=cm.autumn,norm=norm)
lc.set_array(t_sc_f)
lc.set_linewidth(2)
ax.add_collection(lc)

#m.plot(xx,yy,color='red',linewidth=2)
#m.plot(lat_sc[i_sc[:,0]],lon_sc[i_sc[:,0]],color='red',linewidth=2,latlon=True)

# Read in the hdf5 CloudSat data corresponding roughly to StratoClim Flight 7.
CSfi = basedir + 'obs/2C-ICE/2017219064524_59987_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
_, _, lon_2cice, lat_2cice, zeit_2cice = read_cloudsatcalipso_hdf_file(CSfi,'IWC')
xx2,yy2 = m(lon_2cice, lat_2cice)
m.plot(xx2,yy2,color='purple')

#CSfi = 'obs/2C-ICE/2017219181735_59994_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
#_, _, lon_2cice, lat_2cice, zeit_2cice = read_cloudsatcalipso_hdf_file(CSfi,'IWC')
#xx2,yy2 = m(lon_2cice, lat_2cice)
#m.plot(xx2,yy2,color='magenta')

CSfi = basedir + 'obs/2C-ICE/2017220054945_60001_CS_2C-ICE_GRANULE_P1_R05_E06_F01.hdf'
_, _, lon_2cice, lat_2cice, zeit_2cice = read_cloudsatcalipso_hdf_file(CSfi,'IWC')
xx2,yy2 = m(lon_2cice, lat_2cice)
m.plot(xx2,yy2,color='red')

fig.savefig('../output/topography_flight_track_CloudSat-115E.png',bbox_inches='tight',dpi=100)
plt.show(block=True)
