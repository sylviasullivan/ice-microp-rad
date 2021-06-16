import numpy as np
import xarray as xr
from datetime import datetime
from syn_traj_stats import syn_traj_stats

# Binning values - These come from the in-situ data but we have fewer model levels.
#u = 14000
#d = 22000
#n = 70

# Binning values from the vertical grid file
vgrid = xr.open_dataset('/work/bb1018/b380873/vgrid_icon-grid_tropic_55e115e5s40n.nc')
alt = vgrid.vct_a.values[:,0]
j = np.argwhere((alt >= 14000) & (alt <= 22000))
bins_sims = alt[j[:,0]]

# Time range from Lee et al. 2019 (6:20-6:48 UTC)
time0 = datetime(2017, 8, 8, 6, 20)
timef = datetime(2017, 8, 8, 6, 48)

# Do the same for the synthetic simulation trajectories
basedir = '/work/bb1018/b380873/model_output/ICON/'
syntraj1 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V1M0A0R.nc')
syntraj2 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V2M0A0R.nc')
syntraj3 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V2M1A1R.nc')

alt_ICON1 = syntraj1['alt'].sel( time=slice(time0, timef) )
alt_ICON2 = syntraj2['alt'].sel( time=slice(time0, timef) )
alt_ICON3 = syntraj3['alt'].sel( time=slice(time0, timef) )

# Binning in altitude between <u> and <d> with <n> bins, which elements go in which bin?
# np.digitize returns the indices of the bins to which each element in alt* belongs.
# 3 sets of simulations, 1681 times, 625 trajectories
idx = np.zeros((3, alt_ICON1.shape[0], alt_ICON1.shape[1]))
for i in np.arange(alt_ICON1.shape[1]):
    idx[0,:,i] = np.digitize( alt_ICON1[:,i], bins=bins_sims )
    idx[1,:,i] = np.digitize( alt_ICON2[:,i], bins=bins_sims )
    idx[2,:,i] = np.digitize( alt_ICON3[:,i], bins=bins_sims )
    #idx[0,:,i] = np.digitize( alt_ICON1[:,i], bins=np.linspace(u,d,n) )
    #idx[1,:,i] = np.digitize( alt_ICON2[:,i], bins=np.linspace(u,d,n) )
    #idx[2,:,i] = np.digitize( alt_ICON3[:,i], bins=np.linspace(u,d,n) )

stats1 = syn_traj_stats( alt_ICON1, syntraj1.sel( time=slice(time0, timef) ), idx[0], bins_sims, 'qi' )
np.save( '../../output/ICON_syntrajs_0V1M0A0R_stats_qi.npy', stats1 )
stats2 = syn_traj_stats( alt_ICON2, syntraj2.sel( time=slice(time0, timef) ), idx[1], bins_sims, 'qi' )
np.save( '../../output/ICON_syntrajs_0V2M0A0R_stats_qi.npy', stats2 )
stats3 = syn_traj_stats( alt_ICON3, syntraj3.sel( time=slice(time0, timef) ), idx[2], bins_sims, 'qi' )
np.save( '../../output/ICON_syntrajs_0V2M1A1R_stats_qi.npy', stats3 )
