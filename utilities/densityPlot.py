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
    #sys.path.append(os.path.abspath("/work/bb1018/b380873/tropic_vis/traj/"))
    sys.path.append(os.path.abspath("/xdisk/sylvia/tropic_vis/traj/"))
    sys.path.append(os.path.abspath("/xdisk/sylvia/tropic_vis/utilities/"))
    from icetraj import martina_T_IWC_line, martina_T_qi_perc_tropics
    from plotting_utilities import confidence_ellipse

    fs = 15
    ymin = -3
    ymax = 3.5
    alf = 0.5
    #fig, ax = plt.subplots( nrows=3, ncols=2, figsize=(11,8) )
    fig, ax = plt.subplots( nrows=2, ncols=2, figsize=(13,5.5) )  # (12, 5.5) for qiT

    # Load the 2D histogram values for 5 trajectory simulations.
    lbl = [ r"$\bf{(a)}$ $\bf{C}$_1M1T1S", r'$\bf{(b)}$ $\bf{I}$_1M0T0S', r"$\bf{(c)}$ $\bf{C}$_2M1T1S",
            r'$\bf{(d)}$ $\bf{I}$_2M0T0S', r'$\bf{(e)}$ $\bf{C}$_2M0T0S' ]
    farbe = [ cm.Oranges, cm.Oranges, cm.Blues, cm.Blues, cm.Blues ]

    # If this input boolean is True, we add the in-situ IWC-T climatology line.
    if tiwc_line == True:
       T_MK, IWC_10_MK, _, IWC_50_MK, _, IWC_90_MK = martina_T_qi_perc_tropics()
       #T_MK, IWC_min_MK, IWC_max_MK, IWC_med_MK = martina_T_IWC_line()

    # Iterate through the five simulations
    for d, a, l, f, j in zip( datasets, ax.reshape(-1), lbl, farbe, np.arange(4) ):
        if 'T' in histvals:
           h = a.imshow( d[histvals], origin='lower', cmap=f, extent=[190, 240, -3, 3.5],
               vmin=0, vmax=1, aspect=1.65 )
           # four quadrant lines and labels
           a.plot( [210, 237], [np.log10(10), np.log10(10)], linewidth=0.75, linestyle='--', color='red' )
           a.plot( [225, 225], [np.log10(10**(-3)), np.log10(10**3.5)], linewidth=0.75, ls='--', color='red' )
           a.text( 0.9, 0.9, 'I', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.05, 0.73, 'II', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.05, 0.05, 'III', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.9, 0.05, 'IV', c='red', fontsize=fs+2, transform=a.transAxes )
        elif 'RHi' in histvals:
           h = a.imshow( d[histvals], origin='lower', cmap=f, extent=[60, 120, -3, 3.5],
                norm=colors.LogNorm(vmin=10**(-2),vmax=1), aspect=3.8 )
           # four quadrant lines and labels
           a.plot( [60, 120], [np.log10(10), np.log10(10)], linewidth=0.75, linestyle='--', color='red' )
           a.plot( [100, 100], [np.log10(10**(-3)), np.log10(10**3.5)], linewidth=0.75, ls='--', color='red' )
           a.text( 0.9, 0.9, 'I', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.05, 0.73, 'II', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.05, 0.05, 'III', c='red', fontsize=fs+2, transform=a.transAxes )
           a.text( 0.9, 0.05, 'IV', c='red', fontsize=fs+2, transform=a.transAxes )

        if tiwc_line == True:
           a.plot( T_MK, np.log10(IWC_10_MK), linewidth=1.25, linestyle='--', color='gray', alpha=alf )
           a.plot( T_MK, np.log10(IWC_50_MK), linewidth=2, linestyle='-', color='gray', alpha=alf )
           a.plot( T_MK, np.log10(IWC_90_MK), linewidth=1.25, linestyle='--', color='gray', alpha=alf )
           #a.plot( T_MK, np.log10(IWC_min_MK), linewidth=1.25, linestyle='--', color='gray' )
           #a.plot( T_MK, np.log10(IWC_med_MK), linewidth=2, linestyle='-', color='gray' )
           #a.plot( T_MK, np.log10(IWC_max_MK), linewidth=1.25, linestyle='--', color='gray' )

        # If centroid == 0, no metric of central tendency added to the 2D pdf.
        # If centroid == 1, we add means to the 2D pdf.
        if centroid == 1:
           ym, _, xm, _ = centroids( histvals, j )
           # First and third outputs are means. Second and fourth are medians.
           a.scatter( xm, np.log10(ym), marker='x', color='k', s=70, zorder=10 )
        # If centroid == 2, we add medians to the 2D pdf.
        elif centroid == 2:
           _, ym, _, xm = centroids( histvals, j )
           a.scatter( xm, np.log10(ym), marker='x', color='k', s=70, zorder=10 )

        if j == 1 or j == 3:
           cb = fig.colorbar( h, ax=ax[int(j/2),:], shrink=0.75 )
           cb.ax.tick_params(labelsize=fs+3)
           cb.ax.set_ylabel('Probability [%]', fontsize=fs+2)

        a.set_xlim( [xmin, xmax] )
        a.set_ylim( [-3, 3.5] )
        if j == 0 or j == 2:
           a.set_yticks( [ymin, ymin+2, ymin+4, ymin+6] )
           a.set_yticklabels( [r'10$^{-3}$', r'10$^{-1}$', r'10$^{:1d}$'.format(ymin+4),
                   r'10$^{:1d}$'.format(ymin+6)] )
        else:
           a.set_yticks( [ymin+2, ymin+4, ymin+6] )
           a.set_yticklabels( [r'10$^{-1}$', r'10$^{:1d}$'.format(ymin+4),
                   r'10$^{:1d}$'.format(ymin+6)] )
        a.text( 0.02, 1.05, l, fontsize=fs, transform=a.transAxes )
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
#    ax[2,0].set_ylabel( yl, fontsize=fs+3 )
#    ax[2,0].set_xlabel( xl, fontsize=fs+3 )
    ax[1,0].set_xlabel( xl, fontsize=fs+3 )
#    ax[2,1].set_visible('')

#    cb = fig.colorbar(h, ax=ax[2,1])
#    cb.ax.tick_params(labelsize=fs+3)
#    cb.ax.set_ylabel('Probability [%]', fontsize=fs+3)
#    plt.tight_layout()

    if figsave == True:
       print('figsave')
       suffix = ''
       if centroid == True:
          cstr = '_centroid'
       else:
          cstr = ''
       #bd = '/work/bb1018/b380873/tropic_vis/output/'
       bd = '/xdisk/sylvia/tropic_vis/output/'
       fig.savefig( bd + 'CLaMS-ICON-' + yvar + '-' + xvar + '-density-norm' + cstr + suffix + '.pdf' )


# Return the mean or median variables for a variable (or two) denoted by histvals
# and a simulation denoted by sim_tag. sim_tag = [0, 5) with the usual ordering:
# 0 = CLaMS-Tf-0V1M0A0R, 1 = ICON-0V1M0A0R, 2 = CLaMS-Tf-0V2M0A0R, 3 = ICON-0V2M0A0R,
# 4 = CLaMS-0V2M0A0R
def centroids( histvals, sim_tag ):
    import numpy as np

    # One-dimensional
    # These values were calculated from all qi > 10^(-8) g kg-1
    #qi_mean = [ 30.1573, 12.815, 50.063, 17.56, 49.5122, 53.43395 ]
    #qi_med = [ 3.999, 2.513, 10.706, 0.0731, 9.4459, 16.6216 ]

    # These values are calculated only within the bin range in Fig 3
    qi_mean = [ 29.966, 17.825, 56.830, 33.809, 68.958, 67.154 ]
    qi_med = [ 11.308, 11.086, 22.287, 10.035, 26.635, 28.468 ]

    T_mean = [ 229.6377, 228.16, 225.469, 220.01, 225.433, 224.494 ]
    T_med = [ 229.808, 230.37, 225.975, 227.54, 225.906, 225.192 ]

    # These values were calculated from all qi > 10^(-8) g kg-1
    #qi_outflow_mean = [ 25.339, 13.007, 47.548, 24.020, 46.835, 52.0328 ]
    #qi_outflow_med = [ 3.467, 3.349, 10.345, 1.400, 9.0512, 16.7797 ]

    # These values are calculated only within the bin range in Fig 3
    qi_outflow_mean = [ 27.284, 17.3988, 54.297, 32.366, 65.387, 65.477 ]
    qi_outflow_med = [ 10.637, 10.974, 21.759, 9.799, 25.502, 28.628 ]

    T_outflow_mean = [ 228.793, 229.767, 225.244, 229.05, 225.204, 224.513 ]
    T_outflow_med = [ 229.2277, 230.4299, 225.77, 230.159, 225.699, 225.16 ]

    # These values were calculated from all qi > 10^(-8) g kg-1
    #qi_insitu_mean = [ 163.78, 0.0384, 46.214, 0.4429, 44.947, 36.62998 ]
    #qi_insitu_med = [ 56.115, 0.00135, 7.389, 0.0004155, 5.233, 9.1968 ]

    # These values are calculated only within the bin range in Fig 3
    qi_insitu_mean = [ 50.99, 12.967, 68.6266, 48.608, 70.01, 43.801 ]
    qi_insitu_med = [ 37.188, 6.25, 24.265, 12.344, 29.146, 12.386 ]

    T_insitu_mean = [ 206.966, 194.991, 207.71, 198.039, 207.8228, 207.892 ]
    T_insitu_med = [ 207.95, 195.19, 208.349, 198.1322, 208.428, 208.4407 ]

    qi_flight_mean = [ 21.8905, 12.62, 10.213, 4.248, 9.4398, 14.7689 ]
    qi_flight_med = [ 2.6432, 3.799, 1.476, 0.0544, 1.1227, 5.8066 ]

    T_flight_mean = [ 225.061, 226.3257, 220.506, 217.581, 220.4409, 219.554 ]
    T_flight_med = [ 225.061, 226.4384, 220.640, 223.613, 220.603, 219.429 ]

    # The second element of the Ni lists should always be an np.nan as ICNC are not output from ICON-0V1M0A0R
    # These values were calculated from all Ni > 10^(-8) L-1
    #Ni_mean = [ 555.956, np.nan, 692.665, 1132.75, 685.620, 536.199 ]
    #Ni_med = [ 2.869, np.nan, 17.6179, 0.7372, 13.2077, 24.2849 ]

    # These values are calculated only within the bin range in Fig 3
    Ni_mean = [ 871.364, np.nan, 1068.07, 1370.56, 985.43, 683.65 ]
    Ni_med = [ 22.478, np.nan, 64.041, 10.156, 60.737, 57.090 ]

    RHi_mean = [ 96.566, 91.621, 98.556, 74.803, 98.831, 102.335 ]
    RHi_med = [ 97.998, 96.7814, 99.335, 85.245, 99.505, 100.421 ]

    # These values were calculated from all Ni > 10^(-8) L-1
    #Ni_outflow_mean = [ 585.261, np.nan, 664.354, 1632.38, 658.524, 522.13871 ]
    #Ni_outflow_med = [ 2.478, np.nan, 17.2353, 4.839, 12.761, 23.9566 ]

    # These values are calculated only within the bin range in Fig 3
    Ni_outflow_mean = [ 941.58, np.nan, 1019.94, 1476.31, 950.0, 666.59 ]
    Ni_outflow_med = [ 23.15, np.nan, 63.956, 11.162, 59.891, 56.523 ]

    RHi_outflow_mean = [ 97.552, np.nan, 99.273, 89.304, 99.573, 0 ]
    RHi_outflow_med = [ 97.843, np.nan, 99.314, 92.848, 99.474,0 ]

    # These values were calculated from all Ni > 10^(-8) L-1
    #Ni_insitu_mean = [ 44797.53, np.nan, 3379.066, 55.224, 3251.917, 1575.3982 ]
    #Ni_insitu_med = [ 20825.68, np.nan, 154.44, 0.01459, 72.17, 61.09877 ]

    # These values are calculated only within the bin range in Fig 3
    Ni_insitu_mean = [ 46340.43, np.nan, 4281.91, 498.16, 4292.477, 1799.07 ]
    Ni_insitu_med = [ 2.929*10**4, np.nan, 520.436, 2.776, 438.82, 97.4098 ]

    RHi_insitu_mean = [ 67.069, np.nan, 95.849, 32.7418, 95.179, 0 ]
    RHi_insitu_med = [ 74.091, np.nan, 98.545, 28.254, 98.493, 0 ]

    Ni_flight_mean = [ 1834.39, np.nan, 337.60, 61.868, 366.348, 242.3075 ]
    Ni_flight_med = [ 5.993, np.nan, 2.510, 0.4592, 1.5065, 11.0267 ]

    RHi_flight_mean = [ 96.699, np.nan, 101.545, 74.808, 102.00, 0 ]
    RHi_flight_med = [ 97.727, np.nan, 99.9767, 85.416, 100.125, 0 ]

    # Two-dimensional
    qiT_mean = [ [17.64, 10.03, 39.57, 21.36, np.nan, np.nan], [229.15, 227.76, 225.08, 220.56, np.nan, np.nan] ]
    qiT_med = [ [2.18, 1.27, 7.50, 0.064, np.nan, np.nan], [229.13, 230.3, 225.63, 227.63, np.nan, np.nan] ]

    qiT_outflow_mean = [ [15.89, 10.16, 37.72, 22.3, np.nan, np.nan], [228.43, 229.74, 224.95, 228.9, np.nan, np.nan] ]
    qiT_outflow_med = [ [2.033, 2.16, 7.32, 1.18, np.nan, np.nan], [228.6, 230.2, 225.5, 230.04, np.nan, np.nan] ]

    qiT_insitu_mean = [ [46.4, 0.0105, 46.09, 0.436, np.nan, np.nan], [208.98, 194.98, 207.71, 198.04, np.nan, np.nan] ]
    qiT_insitu_med = [ [33.5, 0.00134, 6.43, 0.000416, np.nan, np.nan], [209.3, 195.19, 208.36, 198.13, np.nan, np.nan] ]

    qiT_flight_mean = [ [18.48, 10.0, 9.40, 3.244, np.nan, np.nan], [225.5, 226.98, 220.4677, 217.29, np.nan, np.nan] ]
    qiT_flight_med = [ [1.81, 2.078, 1.409, 0.0324, np.nan, np.nan], [225.4, 226.9, 220.6, 223.12, np.nan, np.nan] ]

    qiRHi_mean = [ [17.64, 10.03, 39.57, 21.36, np.nan, np.nan], [98.66, 88.89, 100.07, 81.17, np.nan, np.nan] ]
    qiRHi_med = [ [2.18, 1.27, 7.50, 0.064, np.nan, np.nan], [98.12, 95.24, 99.43, 81.8, np.nan, np.nan] ]

    qiRHi_outflow_mean = [ [15.89, 10.16, 37.72, 22.3, np.nan, np.nan], [98.75, 92.55, 100.15, 88.05, np.nan, np.nan] ]
    qiRHi_outflow_med = [ [2.033, 2.16, 7.32, 1.18, np.nan, np.nan], [98.08, 96.44, 99.40, 91.94, np.nan, np.nan] ]

    qiRHi_insitu_mean = [ [46.4, 0.0105, 46.09, 0.436, np.nan, np.nan], [95.27, 53.24, 97.55, 33.63, np.nan, np.nan] ]
    qiRHi_insitu_med = [ [33.5, 0.00134, 6.43, 0.000416, np.nan, np.nan], [95.21, 51.4, 98.77, 29.34, np.nan, np.nan] ]

    qiRHi_flight_mean = [ [18.48, 10.0, 9.40, 3.244, np.nan, np.nan], [98.19, 92.39, 101.58, 73.55, np.nan, np.nan] ]
    qiRHi_flight_med = [ [1.81, 2.078, 1.409, 0.0324, np.nan, np.nan], [98.13, 96.22, 99.96, 83.5, np.nan, np.nan] ]

    s = sim_tag
    whichfields = {'qi': [qi_mean[s], qi_med[s]], 'qi_outflow': [qi_outflow_mean[s], qi_outflow_med[s]],
       'qi_insitu': [qi_insitu_mean[s], qi_insitu_med[s]], 'qi_flight': [qi_flight_mean[s], qi_flight_med[s]],
       'T': [T_mean[s], T_med[s]], 'T_outflow': [T_outflow_mean[s], T_outflow_med[s]],
       'T_insitu': [T_insitu_mean[s], T_insitu_med[s]], 'T_flight': [T_flight_mean[s], T_flight_med[s]],
       'Ni': [Ni_mean[s], Ni_med[s]], 'Ni_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s]],
       'Ni_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s]], 'Ni_flight': [Ni_flight_mean[s], Ni_flight_med[s]],
       'RHi': [RHi_mean[s], RHi_med[s]], 'RHi_outflow': [RHi_outflow_mean[s], RHi_outflow_med[s]],
       'RHi_insitu': [RHi_insitu_mean[s], Ni_insitu_med[s]], 'RHi_flight': [RHi_flight_mean[s], RHi_flight_med[s]],
       'qiTh': [qiT_mean[0][s], qiT_med[0][s], qiT_mean[1][s], qiT_med[1][s]],
       'qiTh_outflow': [qiT_outflow_mean[0][s], qiT_outflow_med[0][s], qiT_outflow_mean[1][s], qiT_outflow_med[1][s]],
       'qiTh_insitu': [qiT_insitu_mean[0][s], qiT_insitu_med[0][s], qiT_insitu_mean[1][s], qiT_insitu_med[1][s]],
       'qiTh_flight': [qiT_flight_mean[0][s], qiT_flight_med[0][s], qiT_flight_mean[1][s], qiT_flight_med[1][s]],
       'qiRHih': [qiRHi_mean[0][s], qiRHi_med[0][s], qiRHi_mean[1][s], qiRHi_med[1][s]],
       'qiRHih_outflow': [qiRHi_outflow_mean[0][s], qiRHi_outflow_med[0][s], qiRHi_outflow_mean[1][s], qiRHi_outflow_med[1][s]],
       'qiRHih_insitu': [qiRHi_insitu_mean[0][s], qiRHi_insitu_med[0][s], qiRHi_insitu_mean[1][s], qiRHi_insitu_med[1][s]],
       'qiRHih_flight': [qiRHi_flight_mean[0][s], qiRHi_flight_med[0][s], qiRHi_flight_mean[1][s], qiRHi_flight_med[1][s]],
       'NiTh': [Ni_mean[s], Ni_med[s], T_mean[s], T_med[s]],
       'NiTh_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s], T_outflow_mean[s], T_outflow_med[s]],
       'NiTh_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s], T_insitu_mean[s], T_insitu_med[s]],
       'NiTh_flight': [Ni_flight_mean[s], Ni_flight_med[s], T_flight_mean[s], T_flight_med[s]],
       'NiRHih': [Ni_mean[s], Ni_med[s], RHi_mean[s], RHi_med[s]],
       'NiRHih_outflow': [Ni_outflow_mean[s], Ni_outflow_med[s], RHi_outflow_mean[s], RHi_outflow_med[s]],
       'NiRHih_insitu': [Ni_insitu_mean[s], Ni_insitu_med[s], RHi_insitu_mean[s], RHi_insitu_med[s]],
       'NiRHih_flight': [Ni_flight_mean[s], Ni_flight_med[s], RHi_flight_mean[s], RHi_flight_med[s]] }

    return whichfields[histvals]

