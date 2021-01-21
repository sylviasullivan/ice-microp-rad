import numpy as np
import xarray as xr
import sys, time, os, glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import cartopy.feature
import cartopy.crs as ccrs


timestep = sys.argv[1]
directory = sys.argv[2]
basedir = '/work/bb1018/b380873/traj_output/' + directory + '/'
pi = 3.141592653589793238
os.environ["CARTOPY_USER_BACKGROUNDS"] = "/pf/b/b380873/conda-envs/ncplot/lib/python3.7/site-packages/cartopy/data/raster/natural_earth"

# Define a function to convert radians to degrees.
def rad2deg(x):
    return x*180/pi

# Dimensions here are [patches], [timesteps], [trajs]
fi_list = glob.glob(basedir + 'traj_tst00000' + timestep + '_p*.nc')
# Tells you how many of the files in the directory to visualize.
p1 = 6; p2 = 7
fi_list = fi_list[p1:p2]
patches = len(fi_list)
if directory == 'test2h':
   timesteps = xr.open_dataset(fi_list[0]).dims['time_step']
   numtraj = xr.open_dataset(fi_list[0]).dims['traj_id']
elif directory == 'test24h':
   timesteps = xr.open_dataset(fi_list[0]).dims['time']
   numtraj = xr.open_dataset(fi_list[0]).dims['id']
else:
   timesteps = 3306
   numtraj = xr.open_dataset(fi_list[0]).dims['id']
traj_alt = np.zeros((patches,timesteps,numtraj))
traj_lat = np.zeros((patches,timesteps,numtraj))
traj_lon = np.zeros((patches,timesteps,numtraj))

for j,f in enumerate(fi_list):
    print('patch ' + str(j))
    fi = xr.open_dataset(f)
    alt = fi.alt.values
    lon = fi.lon.values
    lat = fi.lat.values
    t = fi.t.values
    rtime = fi.rtime.values
    #print(alt.shape)

    # Find indices where the matrix != 0.
    xs, ys = np.where(alt != 0)
    # Extract the square with extreme limits.
    # In limited testing, this seems always to generate [=] (88,5308)
    alt = alt[:max(xs)+1,:max(ys)+1]
    lon = lon[:max(xs)+1,:max(ys)+1]
    lat = lat[:max(xs)+1,:max(ys)+1]
    rtime = rtime[:max(xs)+1]
    #print(alt.shape)

    # Store the trimmed matrices.
    print(traj_alt[j].shape)
    print(alt.shape)
    traj_alt[j] = alt/1000.
    traj_lat[j] = rad2deg(lat)
    traj_lon[j] = rad2deg(lon)

plotornot = True
if(plotornot):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(11,11),subplot_kw={'projection':\
         ccrs.PlateCarree()})
    fs = 15
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,color='gray')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':fs}
    gl.ylabel_style = {'size':fs}

    ax.set_title('File p00' + str(p1+1) + ' - 24h',y=1.01)  # str(j+1)
    ax.set_xlabel(r'Latitude [$^{\circ}$N]',fontsize=fs)
    ax.set_ylabel(r'Longitude [$^{\circ}$E]',fontsize=fs)
    #ax.set_extent([76,86,25.5,35],crs=ccrs.PlateCarree())
    #ax.set_extent([80,90,20,30],crs=ccrs.PlateCarree()) # small domain
    ax.set_extent([80,115,18,36],crs=ccrs.PlateCarree()) # large domain
    ax.coastlines()
    ax.background_img(name='BM',resolution='high')
    norm = plt.Normalize(5,22)
    for j in np.arange(patches):
        # How many trajectories to plot?
        n = 200
        # Create a colormap based on altitude.
        #cmap = lambda x : cm.rainbow((x-5.)/17.)
               #(x-np.nanmin(x))/(np.nanmax(x)-np.nanmin(x)))

        for i in np.arange(n):
            # Create a set of line segments to color individually. Points in N x 1 x 2 array.
            points = np.array([traj_lon[j-1,:,i],traj_lat[j-1,:,i]]).T.reshape(-1,1,2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            lc = LineCollection(segments,cmap='rainbow',norm=norm)
            lc.set_array(traj_alt[j-1,:,i])
            lc.set_linewidth(0.5)
            line = ax.add_collection(lc)

    # Read in the flight track from Martina Kraemer's data
    basedir = '/work/bb1018/b380873/tropic_vis/'
    scfi = basedir + 'obs/stratoclim2017.geophysika.0808_1.master.ci_eval.nc'
    sc_data = xr.open_dataset(scfi)
    lat_sc = sc_data['BEST:LAT'].values
    lon_sc = sc_data['BEST:LON'].values
    t_sc = sc_data['time'].values
    i_sc = np.argwhere((~np.isnan(lat_sc)) & (~np.isnan(lon_sc)) & (lat_sc > 0) & (lon_sc > 0))

    # Pulling from https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
    # Create a set of line segments so that we can color them individually
    #xx,yy = m(lon_sc[i_sc[:,0]],lat_sc[i_sc[:,0]])
    #points = np.array([xx,yy]).T.reshape(-1,1,2)
    points = np.array([lon_sc[i_sc[:,0]],lat_sc[i_sc[:,0]]]).T.reshape(-1,1,2)
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


sm = plt.cm.ScalarMappable(cmap='rainbow',norm=norm)
sm.set_array([])
c = plt.colorbar(sm)
c.set_label('Traj. altitude [km]',fontsize=fs)
c.ax.tick_params(labelsize=fs)

fig.savefig('../output/traj_full60h_fast_24hvis' + str(p1+1) + '.png',bbox_inches='tight')
plt.show()
