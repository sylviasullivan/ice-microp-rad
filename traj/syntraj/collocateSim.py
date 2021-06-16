# Input the existing syn_traj Dataset along with time, pressure, lat, and lon from the flight track
#@profile
def collocateSim(syn_traj, var_ICON, flight_time, flight_pressure, flight_lat,
                  flight_lon, flight_alt, sim_acronym):
    import datetime, sys, time
    import numpy as np
    import xarray as xr
    from timeround10 import timeround10

    # Find the closest whole 10-min time from the simulation.
    flight_time_approx = timeround10(flight_time)
    if (flight_time_approx < datetime.datetime(2017, 8, 8, 4, 40)):
        flight_time_approx = datetime.datetime(2017, 8, 8, 4, 40)

    # Find the index corresponding to the ICON pressure level closest to that of the flight track
    basedir = '/work/bb1018/b380873/tropic_vis/remapping/'
    press_ICON = np.loadtxt(basedir + 'PMEAN_48-72.txt')
    i = np.argmin(np.abs(flight_pressure - press_ICON))
    if i < 1 or i > 117:
        raise Exception('Flight pressure outside of simulation range.')

    # Find the index corresponding to the ICON lat/lon closest to that of the flight track
    lat_ICON = var_ICON.lat
    lon_ICON = var_ICON.lon
    j = np.argmin(np.abs(flight_lat - lat_ICON))
    k = np.argmin(np.abs(flight_lon - lon_ICON))

    # Extract the time and lat-lon intervals
    var_ICON = var_ICON.isel( plev=i, lat=j, lon=k )
    var_ICON = var_ICON.sel( time=flight_time_approx )

    # Move the lat, lon, and plev values from coordinates / dimensions to variables
    var_ICON = var_ICON.reset_coords(names=['lat', 'lon', 'plev'], drop=False)

    # Create a new 'time' dimension that only contains the single flight time.
    var_ICON = var_ICON.expand_dims("time").assign_coords(time=("time", [flight_time_approx]))
    syn_traj = xr.concat([syn_traj, var_ICON], dim='time')
    return syn_traj

#profilewrapper()
