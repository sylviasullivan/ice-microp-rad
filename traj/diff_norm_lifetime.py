import xarray as xr
import numpy as np
import time

# Define bins in normalized time
m = 100
tbins_norm = np.linspace( 0., 1., m )

basedir = '/work/bb1018/b380873/traj_output/'
clams_fi = xr.open_dataset( basedir + 'CLAMS-Tf_0V2M0A0R_tst00000450_trim_extract_dt_iwc.nc' )
icon_fi = xr.open_dataset( basedir + 'ICON_0V2M0A0R_tst00000450_trim_extract_dt.nc' )

# How many trajectories are we handling?
n = clams_fi['id'].shape[0]
diff_vals = np.zeros( (n, m) )
diff_vals[:] = np.nan

# Which variable?
#icon_traj = icon_fi['qi']
#clams_traj = clams_fi['qi']
#factor = 10**6

icon_traj0 = icon_fi['Ni'] * icon_fi['rho'] / 1000.
clams_traj0 = clams_fi['Ni'] * clams_fi['rho'] / 1000.
icon_traj = icon_traj0.where( (icon_traj0 > 0) & (clams_traj0 > 0))
clams_traj = clams_traj0.where( (icon_traj0 > 0) & (clams_traj0 > 0))
factor = 1

c = 0
for j in np.arange(n):
    if j%100 == 0:
       print( j )
       np.save( '/scratch/b/b380873/Ni_norm_diff_vals' + str(c) + '_2M.npy', diff_vals )
       c = c + 1

    # Reinitialize this list of differences at every iterationwith the bin number <m>
    diff_list = [ [] for i in np.arange(m) ]
    icon_Ni = icon_traj[:,j]
    clams_Ni = clams_traj[:,j]

    # Find the first and last non-zero element (continuously) along the CLaMS trajectory
    i = np.argwhere( (np.array(clams_Ni) > 0) )

    # Extract the Ni values in this continuous 'time patch' for both CLaMS and ICON
    # only if such a patch exists
    if len(i) > 3:
       start = i[1][0]
       end = i[-1][0]

       # Calculate the Ni differences along this 'time patch'
       clams_Ni = clams_Ni[start:end]
       icon_Ni = icon_Ni[start:end]
       diff = icon_Ni - clams_Ni

       # Calculate a normalized time coordinate
       t_norm = np.arange(len(diff))/float(len(diff))

       # Digitize these t_norm values into bins
       #k = np.digitize( x=t_norm, bins=tbins_norm )
       k = xr.apply_ufunc( np.digitize, t_norm, tbins_norm )

       for elem_idx, group_idx in enumerate(k):
           diff_list[group_idx].append( diff[elem_idx].values*factor )

       for l in np.arange(m):
           if len(diff_list[l]):
              diff_vals[j,l] = np.nanmean( diff_list[l] )
    else:
       continue


np.save( '/scratch/b/b380873/Ni_norm_diff_vals' + str(c) + '_2M.npy', diff_vals )
