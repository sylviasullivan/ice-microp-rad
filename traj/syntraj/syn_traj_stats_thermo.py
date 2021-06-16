# Calculate the mean, median, and standard deviation of T, P, and qv for a given set of simulation trajs
# This function does not fix the number of elements per bin.
def syn_traj_stats( alt_ICON, t_ICON, p_ICON, qv_ICON, indx, bins_sim ):
    import numpy as np
    n = bins_sim.shape[0]

    # Store 3 variables and 3 statistics over 70 levels for 625 trajectories
    stats = np.empty((3, 3, n, alt_ICON.shape[1]))
    for i in np.arange(alt_ICON.shape[1]):
        if i%100 == 0:
           print(i)
        t_list = [ [] for k in np.arange(n) ]
        p_list = [ [] for k in np.arange(n) ]
        qv_list = [ [] for k in np.arange(n) ]

        # Group values for this trajectory into bins.
        for elem_idx, group_idx in enumerate( indx[:,i] ):
            t_list[int(group_idx)-1].append( t_ICON[elem_idx, i].values )
            p_list[int(group_idx)-1].append( p_ICON[elem_idx, i].values )
            qv_list[int(group_idx)-1].append( qv_ICON[elem_idx, i].values )

        for j in np.arange(n):
            stats[0, 0, j, i] = np.nanmean( t_list[j] )
            stats[0, 1, j, i] = np.nanmedian( t_list[j] )
            stats[0, 2, j, i] = np.nanstd( t_list[j] )

            stats[1, 0, j, i] = np.nanmean( p_list[j] )
            stats[1, 1, j, i] = np.nanmedian( p_list[j] )
            stats[1, 2, j, i] = np.nanstd( p_list[j] )

            stats[2, 0, j, i] = np.nanmean( qv_list[j] )
            stats[2, 1, j, i] = np.nanmedian( qv_list[j] )
            stats[2, 2, j, i] = np.nanstd( qv_list[j] )

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    return stats
