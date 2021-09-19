# Generate a 2-d density plot
# suffix is '_outflow' or '_insitu' or blank ('') for both
# xvar is the variable along the x-axis ('T' or 'RHi') and xmin / xmax are its bounds
# yvar is the variable along the y-axis ('qi' or 'Ni') and ymin / ymax are its bounds
# vmax is the upper bound for the colorbar
# flight_track is an optional boolean to use data only around the Flight 7 track
# tiwc_line is an optional boolean to include the in-situ climatological T-IWC line of Kraemer et al
# centroid is an optional boolean to include the centroid of the density plot
# figsave is an optional boolean to save the output plot to a pdf

def densityPlot( suffix, xvar, xmin, xmax, xbin, yvar, ymin, ymax, ybin, vmax,
    flight_track=False, tiwc_line=False, centroid=False, figsave=False ):
    import matplotlib.pyplot as plt
    import sys, os
    import numpy as np
    from matplotlib import cm, colors
    from plotting_utilities import sexy_axes
    sys.path.append(os.path.abspath("/work/bb1018/b380873/tropic_vis/traj/"))
    from icetraj import martina_T_IWC_line

    fs = 18
    fig, ax = plt.subplots( nrows=3, ncols=2, figsize=(13.5,10.5) )

    # We will load values and make density plots for 5 trajectory simulations.
    directories = [ 'traj_CLAMS-Tf_0V1M0A0R', 'traj_ICON_0V1M0A0R', 'traj_CLAMS-Tf_0V2M0A0R',
                    'traj_ICON_0V2M0A0R', 'traj_CLAMS_0V2M0A0R' ]
    timepoints = [ 1350, 1350, 450, 450, 450 ]
    clams_boolean = [ True, False, True, False, True ]
    lbl = [ r"(a) CLaMS 0V1M0A0R-$T'$", '(b) ICON 0V1M0A0R', r"(c) CLaMS 0V2M0A0R-$T'$",
             '(d) ICON 0V2M0A0R', '(e) CLaMS 0V2M0A0R' ]
    farbe = [ cm.Oranges, cm.Oranges, cm.Blues, cm.Blues, cm.Blues ]

    # If this input boolean is True, we add the in-situ IWC-T climatology line.
    if tiwc_line == True:
       T_MK, IWC_min_MK, IWC_max_MK, IWC_med_MK = martina_T_IWC_line()

    # Iterate through the five simulations
    for d, tp, b, a, l, f in zip( directories, timepoints, clams_boolean, ax.reshape(-1), lbl, farbe ):
        ## To recalculate the values being loaded, uncomment the line below
        #xx, rhrh, i = read_iwctraj( d, tp, b, outflow=True )
        #g = np.stack((xx, rhrh, i))
        #np.save('output/qi-RHi-T_' + d + suffix + '.npy',g)
        #

        if flight_track == True:
           xx = np.load('output/qi-RHi-T_flight-track_' + d + suffix + '.npy')
        else:
           xx = np.load('output/qi-RHi-T_' + d + suffix + '.npy')
        yy = np.log10( xx[2] )
        yym = np.log10( np.nanmean( xx[2] ))
        xx = xx[0]
        wgts = np.ones_like(xx) / float(len(xx))

        h = a.hist2d( x=xx, y=yy, cmap=f, bins=[ np.linspace(xmin, xmax, xbin), np.linspace(ymin, ymax, ybin)],
                      weights=wgts )#, norm=colors.Normalize(vmin=0,vmax=vmax))

        print(xx.min(),xx.max(),yy.min(),yy.mean(),yy.max())
        if tiwc_line == True:
           a.plot( T_MK, np.log10(IWC_min_MK), linewidth=1.25, linestyle='--', color='gray' )
           a.plot( T_MK, np.log10(IWC_med_MK), linewidth=2, linestyle='-', color='gray' )
           a.plot( T_MK, np.log10(IWC_max_MK), linewidth=1.25, linestyle='--', color='gray' )

        a.set_yticks([ymin, ymin+2, ymin+4, ymin+6])
        a.set_yticklabels([r'10$^{:2d}$'.format(ymin), r'10$^{:2d}$'.format(ymin+2),
                           r'10$^{:1d}$'.format(ymin+4), r'10$^{:1d}$'.format(ymin+6)])
        a.text(0.05, 0.93, l, weight='bold', fontsize=fs+2, transform=a.transAxes)
        sexy_axes(a, fs+4)

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

    cb = fig.colorbar(h[3], ax=ax[2,1])
    cb.ax.tick_params(labelsize=fs+3)
    cb.ax.set_ylabel('Occurrence density', fontsize=fs+3)
    plt.tight_layout()

    if figsave == True:
       if centroid == True:
          cstr = '_centroid'
       else:
          cstr = ''
       fig.savefig( 'CLaMS-ICON-' + yvar + '-' + xvar + '-density-norm' + cstr + suffix + '.pdf' )
