# Input the existing syn_traj Dataset along with time, pressure, lat, and lon from the flight track
#@profile
def extractSim(syn_traj, var_ICON, flight_time, flight_pressure, flight_lat, flight_lon, flight_alt):
    import datetime, sys
    import numpy as np
    import xarray as xr
    from timeround10 import timeround10
    ri = np.random.randint

    # Global variables in the script form
    sim_acronym = '0V1M0A0R'
    n = 5
    ll_interval = 0.75
    alt_interval = 1

    # Find the nearest whole 10-min time.
    flight_time_approx = timeround10(flight_time)

    # Construct the time window to extract.
    early_time = flight_time_approx - datetime.timedelta(minutes=30)
    late_time = flight_time_approx + datetime.timedelta(minutes=30)

    # Find indices corresponding to ICON pressure levels above + below the closest flight track match.
    basedir = '/work/bb1018/b380873/tropic_vis/remapping/'
    sim_pressures = np.loadtxt(basedir + 'PMEAN_48-72.txt')
    i = np.argmin(np.abs(flight_pressure - sim_pressures))
    if i < 1 or i > 117:
        raise Exception('Flight pressure outside of simulation range.')
    var_ICON = var_ICON.isel( plev=slice(i-alt_interval, i+alt_interval+1) )

    # Extract the time and lat-lon intervals
    var_ICON = var_ICON.sel( time=slice(early_time, late_time),
                             lat=slice(flight_lat-ll_interval, flight_lat+ll_interval),
                             lon=slice(flight_lon-ll_interval, flight_lon+ll_interval) )

    # Randomly generate <n> indices along each axis and save these from the var_ICON structure into syn_traj.
    tv = var_ICON['qv'].shape
    guy = var_ICON.isel(time=ri(0,tv[0],size=n), plev=ri(0,tv[1],size=n), lat=ri(0,tv[2],size=n), lon=ri(0,tv[3],size=n))

    # Merge the lat, lon, plev, and time dimensions into a new 'ntraj' dimension.
    guy = guy.stack(ntraj=("lat", "lon", "plev", "time"))
    guy['ntraj'] = np.arange(1, n**4+1)

    # Create a new 'time' dimension that only contains the single flight time.
    guy = guy.expand_dims("time").assign_coords(time=("time", [flight_time]))

    # Add lat, lon, and alt fields to guy that contain the single flight lat, lon, alt.
    guy['lat'] = (["ntraj", "time"], flight_lat*np.ones((n**4, 1)))
    guy['lon'] = (["ntraj", "time"], flight_lon*np.ones((n**4, 1)))
    guy['alt'] = (["ntraj", "time"], flight_alt*np.ones((n**4, 1)))
    syn_traj = xr.concat([syn_traj, guy], dim='time')
    return syn_traj

#profilewrapper()
