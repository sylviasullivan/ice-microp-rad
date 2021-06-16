# Calculate the mean, median, and standard deviation of <var> in a given set of simulation trajs
# This function does not fix the number of elements per bin.
# var_names = ['temp', 'air_pressure', 'qv', 'qc', 'qi', 'qs', 'qg', 'clc', 'omega', 'alt', 'lon', 'lat']

def syn_traj_stats( alt_ICON, fi_ICON, indx, bins_sim, var ):

    import numpy as np

    # Store 3 statistics over <n> bins for 625 trajectories
    n = bins_sim.shape[0]
    stats = np.empty((3, n, alt_ICON.shape[1]))
    stats[:] = np.nan

    for i in np.arange(alt_ICON.shape[1]):
        print(i)
        var_list = [ [] for k in np.arange(n) ]

        # Group <var> values along this trajectory into bins
        for elem_idx, group_idx in enumerate( indx[:, i] ):
            var_list[int(group_idx-1)].append( fi_ICON[var][elem_idx, i].item() )

        # Calculate statistics in each bin
        for k in np.arange(n):
            stats[0, k, i] = np.nanmean( var_list[k] )
            stats[1, k, i] = np.nanmedian( var_list[k] )
            stats[2, k, i] = np.nanstd( var_list[k] )

    return stats
