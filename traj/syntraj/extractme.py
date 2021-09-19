import numpy as np
import xarray as xr
import sys, time
import pandas as pd
import datetime
from extractSim import extractSim

# Global variables in the script form
sim_acronym = '0V1M0A0R'
n = 5**4
ll_interval = 0.75
alt_interval = 1

# Load the observational data
basedir = '/work/bb1018/b380873/tropic_vis/obs/'
fi = basedir + 'stratoclim2017.geophysika.0808_1.filtered_per_sec.nc'
Stratoclim = xr.open_dataset(fi)
flight_times = Stratoclim['time']

# <j> is the first iteration for which there are ICON high-resolution values available.
j = 1942 # 1342
#tt = flight_times.shape[0] - j
tt = 0

# Load the ICON values
ICON = xr.open_dataset('/work/bb1018/b380873/model_output/ICON/ICON_3D_F10MIN_icon_tropic_' + sim_acronym + '_PL2.nc')

# Initiate the synthetic trajectory Dataset
syn_traj = xr.Dataset( data_vars=dict(
                            omega=( ["time", "ntraj"], np.empty([tt, n]) ),
                            temp=( ["time", "ntraj"], np.empty([tt, n]) ),
                            air_pressure=( ["time", "ntraj"], np.empty([tt, n]) ),
                            qv=( ["time", "ntraj"], np.empty([tt, n]) ),
                            qi=( ["time", "ntraj"], np.empty([tt, n]) ),
                            qc=( ["time", "ntraj"], np.empty([tt, n]) ),
                            qs=( ["time", "ntraj"], np.empty([tt, n]) ),
                            qg=( ["time", "ntraj"], np.empty([tt, n]) ),
                            clc=( ["time", "ntraj"], np.empty([tt, n]) ),
                            lat=( ["time", "ntraj"], np.empty([tt, n]) ),
                            lon=( ["time", "ntraj"], np.empty([tt, n]) ),
                            alt=( ["time", "ntraj"], np.empty([tt, n]) ),
                        ),
                       coords=dict(
                           time=[], ntraj=np.arange(1, n+1)) #time=flight_times[j:]
                     )

# Set the variable attributes as in the standard ICON output file.
for v in syn_traj.data_vars:
    if v != 'alt':
       syn_traj[v].attrs["long_name"] = ICON[v].long_name
       syn_traj[v].attrs["units"] = ICON[v].units
       syn_traj[v].attrs["standard_name"] = ICON[v].standard_name
    else:
       syn_traj[v].attrs["long_name"] = "Altitude"
       syn_traj[v].attrs["units"] = "m"

syn_traj['ntraj'].attrs["long_name"] = 'Trajectory ID'

for flight_iter, flight_time in enumerate(flight_times[j:]):
    if flight_iter%500 == 0:
       print(flight_iter)

    # timestamp-datetime-timedoodadoo conversion from hell....
    flight_time_not_np = pd.to_datetime(flight_time.values).to_pydatetime()
    flight_pressure = Stratoclim['BEST:PRESS'].sel(time=flight_time_not_np).values*100 # [Pa]
    flight_lat = Stratoclim['BEST:LAT'].sel(time=flight_time_not_np).values
    flight_lon = Stratoclim['BEST:LON'].sel(time=flight_time_not_np).values
    flight_alt = Stratoclim['BEST:ALT'].sel(time=flight_time_not_np).values

    # Based on the flight values, load the relevant chunk of simulations
    syn_traj = extractSim(syn_traj, ICON, flight_time_not_np, flight_pressure, flight_lat, flight_lon, flight_alt)

syn_traj.to_netcdf(path='/work/bb1018/b380873/model_output/ICON/ICON_synthetic_trajs_' + sim_acronym + '_3.nc')
