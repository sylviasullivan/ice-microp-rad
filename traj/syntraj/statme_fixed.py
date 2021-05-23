import numpy as np
import xarray as xr
from datetime import datetime
from syn_traj_stats_fixed import syn_traj_stats_fixed

# Binning values
u = 14000
d = 22000
n = 70

# Time range from Lee et al. 2019 (6:20-6:48 UTC)
time0 = datetime(2017, 8, 8, 6, 20)
timef = datetime(2017, 8, 8, 6, 48)

# Do the same for the synthetic simulation trajectories
basedir = '/work/bb1018/b380873/model_output/ICON/'
syntraj1 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V1M0A0R.nc')
syntraj2 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V2M0A0R.nc')
syntraj3 = xr.open_dataset(basedir + 'ICON_synthetic_trajs_0V2M1A1R.nc')

alt_ICON1 = syntraj1['alt'].sel( time=slice(time0, timef) )
temp_ICON1 = syntraj1['temp'].sel( time=slice(time0, timef) )
press_ICON1 = syntraj1['air_pressure'].sel( time=slice(time0, timef) )

alt_ICON2 = syntraj2['alt'].sel( time=slice(time0, timef) )
temp_ICON2 = syntraj2['temp'].sel( time=slice(time0, timef) )
press_ICON2 = syntraj2['air_pressure'].sel( time=slice(time0, timef) )

alt_ICON3 = syntraj3['alt'].sel( time=slice(time0, timef) )
temp_ICON3 = syntraj3['temp'].sel( time=slice(time0, timef) )
press_ICON3 = syntraj3['air_pressure'].sel( time=slice(time0, timef) )

# Convert kg kg-1 to ppmv for the mixing ratios
mw_dryair = 28.97*1000    # kg air (mol air)-1
mw_watvap = 18.02*1000    # kg wv (mol wv)-1
conv = mw_dryair / mw_watvap

qv_ICON1 = syntraj1['qv'].sel( time=slice(time0, timef) ) * conv * 10**6
qv_ICON2 = syntraj2['qv'].sel( time=slice(time0, timef) ) * conv * 10**6
qv_ICON3 = syntraj3['qv'].sel( time=slice(time0, timef) ) * conv * 10**6

# Binning in altitude between <u> and <d> with <n> bins, which elements go in which bin?
# np.digitize returns the indices of the bins to which each element in alt* belongs.
# 3 sets of simulations, 1681 times, 625 trajectories
idx = np.zeros((3, alt_ICON1.shape[0], alt_ICON1.shape[1]))
for i in np.arange(alt_ICON1.shape[1]):
    idx[0,:,i] = np.digitize( alt_ICON1[:,i], bins=np.linspace(u,d,n) )
    idx[1,:,i] = np.digitize( alt_ICON2[:,i], bins=np.linspace(u,d,n) )
    idx[2,:,i] = np.digitize( alt_ICON3[:,i], bins=np.linspace(u,d,n) )

stats1_fixed = syn_traj_stats_fixed( alt_ICON1, temp_ICON1, press_ICON1, qv_ICON1, idx[0] )
np.save('../../output/ICON_syntrajs_0V1M0A0R_stats_fixed.npy', stats1_fixed)
stats2_fixed = syn_traj_stats_fixed( alt_ICON2, temp_ICON2, press_ICON2, qv_ICON2, idx[1] )
np.save('../../output/ICON_syntrajs_0V2M0A0R_stats_fixed.npy', stats2_fixed)
stats3_fixed = syn_traj_stats_fixed( alt_ICON3, temp_ICON3, press_ICON3, qv_ICON3, idx[2] )
np.save('../../output/ICON_syntrajs_0V2M1A1R_stats_fixed.npy', stats3_fixed)
