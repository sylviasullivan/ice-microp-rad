#def lagrangianComp(obsFile, simFile):
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys, time
import datetime
import pandas as pd

# Which simulation and variable do you want to look at?
global sim_acronym
global sim_var

sim_acronym = '0V2M0A0R'
sim_var = 'qv'

# Which file to open based upon which variable you want to look at
def vardictionary(var_key):
    vd = {'temp'         : 'TPOMEGA_3D_F10MIN_icon_tropic_' + sim_acronym + '_PL2.nc',
          'qv'           : 'CLCONV_3D_F10MIN_icon_tropic_' + sim_acronym + '_PL2.nc',
          'qi'           : 'CLCONV_3D_F10MIN_icon_tropic_' + sim_acronym + '_PL2.nc',
          'air_pressure' : 'TPOMEGA_3D_F10MIN_icon_tropic_' + sim_acronym + '_PL2.nc',
         }
    return vd[var_key]

# Input a np.datetime64 and round it to the nearest 10 minutes.
def timeround10(dt):
    dt_not_np = pd.to_datetime(dt)
    b = round(dt_not_np.minute,-1)
    if b == 60:
        return_time = datetime.datetime(2017, 8, 8, dt_not_np.hour + 1, 0)
    else:
        return_time = datetime.datetime(2017, 8, 8, dt_not_np.hour, int(b))
    return return_time

# Input the time, pressure, lat, and lon from the flight track
def extractSim(flight_time, flight_pressure, flight_lat, flight_lon):
    # Translate the variable key to its corresponding file
    fi = vardictionary(sim_var)
    sim_dir = '/work/bb1018/b380873/model_output/ICON/'
    var_ICON = xr.open_dataset(sim_dir + fi)[sim_var]

    # Find the nearest whole 10-min time.
    flight_time_approx = timeround10(flight_time.values)

    # Construct the time window to extract.
    early_time = flight_time_approx - datetime.timedelta(minutes=30)
    late_time = flight_time_approx + datetime.timedelta(minutes=30)

    # Find the indices for levels above and below the closest match.
    basedir = '/work/bb1018/b380873/tropic_vis/remapping/'
    sim_pressures = np.loadtxt(basedir + 'PMEAN_48-72.txt')
    i = np.argmin(np.abs(flight_pressure - sim_pressures))
    if i < 1 or i > 117:
        raise Exception('Flight pressure outside of simulation range.')
    var_ICON = var_ICON.isel(plev=slice(i-1,i+2))

    # Define the lat-lon interval to extract
    ll_interval = 0.75
    var_ICON = var_ICON.sel(time=slice(early_time, late_time),
                            lat=slice(flight_lat-ll_interval, flight_lat+ll_interval),
                             lon=slice(flight_lon-ll_interval, flight_lon+ll_interval))

    return var_ICON


# Load the observational data
basedir = '/work/bb1018/b380873/tropic_vis/obs/'
fi = basedir + 'stratoclim2017.geophysika.0808_1.filtered_per_sec.nc'
Stratoclim = xr.open_dataset(fi)
flight_times = Stratoclim['time']

# Initialize <n> number of synthetic trajectories. These must hold at least altitude and 
n = 1000
first_val = False
syn_traj = np.zeros((2,flight_times.shape[0], n))
for flight_iter, flight_time in enumerate(flight_times):
    print(flight_iter)
    flight_pressure = Stratoclim['BEST:PRESS'].sel(time=flight_time).values*100 # [Pa]
    flight_lat = Stratoclim['BEST:LAT'].sel(time=flight_time).values
    flight_lon = Stratoclim['BEST:LON'].sel(time=flight_time).values

    # Based on the flight values, load the relevant chunk of simulations
    var_ICON = extractSim(flight_time, flight_pressure, flight_lat, flight_lon)

    # Unwrap the simulation values and extract a random sample of n values
    if var_ICON.shape[0] != 0:
       # This statement stores the iteration from which we should save syn_trajs
       if first_val == False:
            sim_iter = flight_iter
            print('sim_iter: ' + str(sim_iter))
            first_val = True

       syn_traj[flight_iter] = np.random.choice(np.ravel(var_ICON), size=n, replace=False)

# Trim the leading zeros
syn_traj = syn_traj[sim_iter:]
np.save('../output/' + sim_var + '-syntraj-' + sim_acronym + '.npy', syn_traj)
# Output these synthetic trajectories to a nc file
#ds = xr.Dataset(
#     { sim_var: (("n", "t")), syn_traj) },
#     coords = {
#          "n": np.arange(n)

