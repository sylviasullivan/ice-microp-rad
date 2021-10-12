# Generate a 2-d density plot
# histvals is 'qiTh', 'qiRHih', 'NiTh_flight', 'NiRHih_instiu', etc
# xvar is the variable along the x-axis ('T' or 'RHi') and xmin / xmax are its bounds
# yvar is the variable along the y-axis ('qi' or 'Ni')
# tiwc_line is an optional boolean to include the in-situ climatological T-IWC line of Kraemer et al
# centroid is an optional boolean to include the centroid of the density plot
# figsave is an optional boolean to save the output plot to a pdf

def densityPlot( datasets, histvals, xvar, xmin, xmax, yvar, tiwc_line=False,
                    centroid=False, figsave=False ):
    import matplotlib.pyplot as plt
    import sys, os
    import numpy as np
    from matplotlib import cm, colors
    from plotting_utilities import sexy_axes
    sys.path.append(os.path.abspath("/work/bb1018/b380873/tropic_vis/traj/"))
    from icetraj import martina_T_IWC_line, martina_T_qi_perc_tropics

    fs = 15
    ymin = -3
    ymax = 3.5
    alf = 0.5
    fig, ax = plt.subplots( nrows=3, ncols=2, figsize=(11,8) )

    # Load the 2D histogram values for 5 trajectory simulations.
    lbl = [ r"(a) CLaMS 0V1M0A0R-$T'$", '(b) ICON 0V1M0A0R', r"(c) CLaMS 0V2M0A0R-$T'$",
             '(d) ICON 0V2M0A0R', '(e) CLaMS 0V2M0A0R' ]
    farbe = [ cm.Oranges, cm.Oranges, cm.Blues, cm.Blues, cm.Blues ]

    # If this input boolean is True, we add the in-situ IWC-T climatology line.
    if tiwc_line == True:
       T_MK, IWC_10_MK, _, IWC_50_MK, _, IWC_90_MK = martina_T_qi_perc_tropics()
       #T_MK, IWC_min_MK, IWC_max_MK, IWC_med_MK = martina_T_IWC_line()

    # Iterate through the five simulations
    for d, a, l, f, j in zip( datasets, ax.reshape(-1), lbl, farbe, np.arange(5) ):
        if 'T' in histvals:
           h = a.imshow( d[histvals], origin='lower', cmap=f, extent=[190, 240, -3, 3.5],
               vmin=0, vmax=1, aspect=1.65 )
        elif 'RHi' in histvals:
           h = a.imshow( d[histvals], origin='lower', cmap=f, extent=[60, 120, -3, 3.5],
                norm=colors.LogNorm(vmin=10**(-2),vmax=1), aspect=3 )

        if tiwc_line == True:
           a.plot( T_MK, np.log10(IWC_10_MK), linewidth=1.25, linestyle='--', color='gray', alpha=alf )
           a.plot( T_MK, np.log10(IWC_50_MK), linewidth=2, linestyle='-', color='gray', alpha=alf )
           a.plot( T_MK, np.log10(IWC_90_MK), linewidth=1.25, linestyle='--', color='gray', alpha=alf )
           #a.plot( T_MK, np.log10(IWC_min_MK), linewidth=1.25, linestyle='--', color='gray' )
           #a.plot( T_MK, np.log10(IWC_med_MK), linewidth=2, linestyle='-', color='gray' )
           #a.plot( T_MK, np.log10(IWC_max_MK), linewidth=1.25, linestyle='--', color='gray' )

        # If this input boolean is True, we add centroids to the 2D pdf.
        if centroid == True:
           ym, _, xm, _ = centroids( histvals, j )
           # First and third outputs are means. Second and fourth are medians.
           a.scatter( xm, np.log10(ym), marker='x', color='k', s=70, zorder=10 )

        a.set_xlim( [xmin, xmax] )
        a.set_ylim( [-3, 3.5] )
        a.set_yticks( [ymin, ymin+2, ymin+4, ymin+6] )
        a.set_yticklabels( [r'10$^{-3}$', r'10$^{-1}$', r'10$^{:1d}$'.format(ymin+4),
                r'10$^{:1d}$'.format(ymin+6)] )
        a.text( 0.05, 0.93, l, weight='bold', fontsize=fs, transform=a.transAxes )
        sexy_axes( a, fs+2 )

    if xvar == 'T':
       xl = "Temperature [K]"
    if xvar == 'RHi':
       xl = r"RH$_i$ [%]"
    if yvar == 'qi':
       yl = r"$q_i$ [ppmv]"
    if yvar == 'Ni':
       yl = r"$N_i$ [L$^{-1}$]"

    ax[0,0].set_ylabel( yl, fontsize=fs+3 )
    ax[1,0].set_ylabel( yl, fontsize=fs+3 )
    ax[1,1].set_xlabel( xl, fontsize=fs+3 )
    ax[2,0].set_ylabel( yl, fontsize=fs+3 )
    ax[2,0].set_xlabel( xl, fontsize=fs+3 )
    ax[2,1].set_visible('')

    cb = fig.colorbar(h, ax=ax[2,1])
    cb.ax.tick_params(labelsize=fs+3)
    cb.ax.set_ylabel('Probability [%]', fontsize=fs+3)
    plt.tight_layout()

    if figsave == True:
       if centroid == True:
          cstr = '_centroid'
       else:
          cstr = ''
       fig.savefig( 'CLaMS-ICON-' + yvar + '-' + xvar + '-density-norm' + cstr + suffix + '.pdf' )


# Return the mean or median variables for a variable (or two) denoted by histvals
# and a simulation denoted by sim_tag. sim_tag = [0, 5) with the usual ordering:
# 0 = CLaMS-Tf-0V1M0A0R, 1 = ICON-0V1M0A0R, 2 = CLaMS-Tf-0V2M0A0R, 3 = ICON-0V2M0A0R,
# 4 = CLaMS-0V2M0A0R
def centroids( histvals, sim_tag ):
    import numpy as np

    # One-dimensional
    qi_mean = [ 30.1573, 12.815, 50.063, 17.56, 49.5122, 53.43395 ]
    qi_med = [ 3.999, 2.513, 10.706, 0.0731, 9.4459, 16.6216 ]

    T_mean = [ 229.6377, 228.16, 225.469, 220.01, 225.433 ]
    T_med = [ 229.808, 230.37, 225.975, 227.54, 225.906 ]

    qi_outflow_mean = [ 25.339, 13.007, 47.548, 24.020, 46.835, 52.0328 ]
    qi_outflow_med = [ 3.402, 3.349, 10.345, 1.400, 9.0512, 16.7797 ]

    T_outflow_mean = [ 228.793, 229.767, 225.244, 229.05, 225.204 ]
    T_outflow_med = [ 229.2277, 230.4299, 225.77, 230.159, 225.699 ]

    qi_insitu_mean = [ 163.78, 0.0384, 46.214, 0.4429, 44.947, 36.62998 ]
    qi_insitu_med = [ 56.115, 0.00135, 7.389, 0.0004155, 5.233, 9.1968 ]

    T_insitu_mean = [ 206.966, 194.991, 207.71, 198.039, 207.8228 ]
    T_insitu_med = [ 207.95, 195.19, 208.349, 198.1322, 208.428 ]

    qi_flight_mean = [ 21.8905, 12.62, 10.213, 4.248, 9.4398, 14.7689 ]
    qi_flight_med = [ 2.6432, 3.799, 1.476, 0.0544, 1.1227, 5.8066 ]

    T_flight_mean = [ 225.061, 226.3257, 220.506, 217.581, 220.4409 ]
    T_flight_med = [ 225.061, 226.4384, 220.640, 223.613, 220.603 ]

    # The second element of the Ni lists should always be an np.nan as ICNC are not output from ICON-0V1M0A0R
    Ni_mean = [ 555.956, np.nan, 692.665, 1132.75, 685.620, 536.199 ]
    Ni_med = [ 2.869, np.nan, 17.6179, 0.7372, 13.2077, 24.2849 ]

    # Fourth element here does not seem right..
    RHi_mean = [ 97.639, np.nan, 99.225, 90, 99.510 ]
    RHi_med = [ 98.073, np.nan, 99.3655, 98, 99.538 ]

    Ni_outflow_mean = [ 585.261, np.nan, 664.354, 1632.38, 658.524, 522.13871 ]
    Ni_outflow_med = [ 2.478, np.nan, 17.2353, 4.839, 12.761, 23.9566 ]

    RHi_outflow_mean = [ 97.552, np.nan, 99.273, 89.304, 99.573 ]
    RHi_outflow_med = [ 97.843, np.nan, 99.314, 92.848, 99.474 ]

    Ni_insitu_mean = [ 44797.53, np.nan, 3379.066, 55.224, 3251.917, 1575.3982 ]
    Ni_insitu_med = [ 20825.68, np.nan, 154.44, 0.01459, 72.17, 61.09877 ]

    RHi_insitu_mean = [ 67.069, np.nan, 95.849, 32.7418, 95.179 ]
    RHi_insitu_med = [ 74.091, np.nan, 98.545, 28.254, 98.493 ]

    Ni_flight_mean = [ 1834.39, np.nan, 3379.600, 61.868, 366.348, 242.3075 ]
    Ni_flight_med = [ 5.993, np.nan, 2.510, 0.4592, 1.5065, 11.0267 ]

    RHi_flight_mean = [ 96.699, np.nan, 101.545, 74.808, 102.00 ]
    RHi_flight_med = [ 97.727, np.nan, 99.9767, 85.416, 100.125 ]

    s = sim_tag
    whichfields = {'qi': [qi_mean[s], qi_med[s]], 'qi_outflow': [qi_outflow_mean[s], qi_outflow_med[s]],
       'qi_insitu': [qi_insitu_mean[s], qi_insitu_med[s]], 'qi_flight': [qi_flight_mean[s], qi_flight_med[s]],
       'T': [T_mean[s], T_med[s]], 'T_outflow': [T_outflow_mean[s], T_outflow_med[s]],
       'T_insitu': [T_insitu_mean[s], T_insitu_med[s]], 'T_flight': [T_flight_mean[s], T_flight_med[s]],
       'Ni': [Ni_mean[s], Ni_med[s]], 'Ni_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s]],
       'Ni_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s]], 'Ni_flight': [Ni_flight_mean[s], Ni_flight_med[s]],
       'RHi': [RHi_mean[s], RHi_med[s]], 'RHi_outflow': [RHi_outflow_mean[s], RHi_outflow_med[s]],
       'RHi_insitu': [RHi_insitu_mean[s], Ni_insitu_med[s]], 'RHi_flight': [RHi_flight_mean[s], RHi_flight_med[s]],
       'qiTh': [qi_mean[s], qi_med[s], T_mean[s], T_med[s]],
       'qiTh_outflow': [qi_outflow_mean[s], qi_outflow_med[s], T_outflow_mean[s], T_outflow_med[s]],
       'qiTh_insitu': [qi_insitu_mean[s], qi_insitu_med[s], T_insitu_mean[s], T_insitu_med[s]],
       'qiTh_flight': [qi_flight_mean[s], qi_flight_med[s], T_flight_mean[s], T_flight_med[s]],
       'qiRHih': [qi_mean[s], qi_med[s], RHi_mean[s], RHi_med[s]],
       'qiRHih_outflow': [qi_outflow_mean[s], qi_outflow_med[s], RHi_outflow_mean[s], RHi_outflow_med[s]],
       'qiRHih_insitu': [qi_insitu_mean[s], qi_insitu_med[s], RHi_insitu_mean[s], RHi_insitu_med[s]],
       'qiRHih_flight': [qi_flight_mean[s], qi_flight_med[s], RHi_flight_mean[s], RHi_flight_med[s]],
       'NiTh': [Ni_mean[s], Ni_med[s], T_mean[s], T_med[s]],
       'NiTh_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s], T_outflow_mean[s], T_outflow_med[s]],
       'NiTh_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s], T_insitu_mean[s], T_insitu_med[s]],
       'NiTh_flight': [Ni_flight_mean[s], Ni_flight_med[s], T_flight_mean[s], T_flight_med[s]],
       'NiRHih': [Ni_mean[s], Ni_med[s], RHi_mean[s], RHi_med[s]],
       'NiRHih_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s], RHi_outflow_mean[s], RHi_outflow_med[s]],
       'NiRHih_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s], RHi_insitu_mean[s], RHi_insitu_med[s]],
       'NiRHih_flight': [Ni_flight_mean[s], Ni_flight_med[s], RHi_flight_mean[s], RHi_flight_med[s]] }

    return whichfields[histvals]

