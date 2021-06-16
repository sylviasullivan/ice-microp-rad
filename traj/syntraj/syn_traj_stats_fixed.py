# Calculate the mean, median, and standard deviation of <var> in a given set of simulation trajs
# This function fixes the number of elements per bin.
# var_names = ['temp', 'air_pressure', 'qv', 'qc', 'qi', 'qs', 'qg', 'clc', 'omega', 'alt', 'lon', 'lat']

def syn_traj_stats_fixed( alt_ICON, fi_ICON, indx, bins_sim, var ):

    import numpy as np
    import random

    # Store 3 statistics over <n> bins for 625 trajectories
    n = bins_sim.shape[0]
    stats = np.empty((3, n, alt_ICON.shape[1]))
    stats[:] = np.nan

    # Read in the number of in-situ measurements in each bin
    basedir = '/work/bb1018/b380873/tropic_vis/'
    Stratoclim_temp_len = np.load(basedir + 'output/Stratoclim_temp_len.npy')
    Stratoclim_qv_len = np.load(basedir + 'output/Stratoclim_qv_len.npy')

    for i in np.arange(alt_ICON.shape[1]):
        #if i%50 == 0:
        print(i)
        var_list = [ [] for i in np.arange(n) ]
        #t_list = [ [] for i in np.arange(n) ]
        #p_list = [ [] for i in np.arange(n) ]
        #qv_list = [ [] for i in np.arange(n) ]

        # Group <var> values along this trajectory into bins.
        for elem_idx, group_idx in enumerate( indx[:,i] ):
            var_list[int(group_idx-1)].append( fi_ICON[var][elem_idx, i].item() )
            #t_list[int(group_idx)-1].append( t_ICON[elem_idx, i].item() )
            #p_list[int(group_idx)-1].append( p_ICON[elem_idx, i].item() )
            #qv_list[int(group_idx)-1].append( qv_ICON[elem_idx, i].item() )

        # Only retain the corresponding number of elements from the measurements
        # We do this calculation only for n-1 bins as the last bin contains the whole "upper-level" trajectory set (z ~ 22 km)
        var_list_fixed = [ [] for i in np.arange(n) ]
        #t_list_fixed = [ [] for i in np.arange(n) ]
        #p_list_fixed = [ [] for i in np.arange(n) ]
        #qv_list_fixed = [ [] for i in np.arange(n) ]
        for j in np.arange(n-1):
            #print(len(t_list[j]))
            #print(Stratoclim_temp_len[j])
            if var_list[j]:
                var_list_fixed[j] = random.sample( var_list[j], Stratoclim_temp_len[j] )
                # To be revisited above, whether Stratoclim_temp_len[j] is relevant for all variables
            #if t_list[j]:
            #    t_list_fixed[j] = random.sample( t_list[j], Stratoclim_temp_len[j] )
            #if p_list[j]:
            #    p_list_fixed[j] = random.sample( p_list[j], Stratoclim_temp_len[j] )
            #if qv_list[j]:
            #    qv_list_fixed[j] = random.sample( qv_list[j], Stratoclim_qv_len[j] )

            if var_list_fixed[j]:
                stats[0, j, i] = np.nanmean( var_list_fixed[j] )
                stats[1, j, i] = np.nanmedian( var_list_fixed[j] )
                stats[2, j, i] = np.nanstd( var_list_fixed[j] )
            #if t_list_fixed[j]:
            #    stats[0, 0, j, i] = np.nanmean( t_list_fixed[j] )
            #    stats[0, 1, j, i] = np.nanmedian( t_list_fixed[j] )
            #    stats[0, 2, j, i] = np.nanstd( t_list_fixed[j] )
            #if p_list_fixed[j]:
            #    stats[1, 0, j, i] = np.nanmean( p_list_fixed[j] )
            #    stats[1, 1, j, i] = np.nanmedian( p_list_fixed[j] )
            #    stats[1, 2, j, i] = np.nanstd( p_list_fixed[j] )
            #if qv_list_fixed[j]:
            #    stats[2, 0, j, i] = np.nanmean( qv_list_fixed[j] )
            #    stats[2, 1, j, i] = np.nanmedian( qv_list_fixed[j] )
            #    stats[2, 2, j, i] = np.nanstd( qv_list_fixed[j] )

    return stats
