# Define the colors and linestyles to use for different simulations
# Color-blind colormap used from https://www.nature.com/articles/nmeth.1618/figures/2
def sim_colors():
    return {'0V1M0A0R': (0.88, 0.3686, 0), '0V2M0A0R': (0, 0.4470, 0.698), #0.8353
            '0V2M0A1R': (0.9411, 0.8941, 0.2588), '0V2M1A0R': (0, 0.6196, 0.4509),
            '0V2M1A1R': (0.9019, 0.6235, 0), '1V1M0A0R': (0.88, 0.3686, 0),
            '1V2M0A0R': (0, 0.4470, 0.698), '1V2M0A1R': (0.9411, 0.8941, 0.2588),
            '1V2M1A0R': (0, 0.6196, 0.4509), '1V2M1A1R': (0.9019, 0.6235, 0),
            'MPI': (0.8, 0.4745, 0.6549), 'ERA': (0, 0, 0),
            'CloudSat': (0.3373, 0.7058, 0.9137), 'CERES': (0.3373, 0.7058, 0.9137),
            'POSIDON': (0.3373, 0.7058, 0.9137), 'ATTREX': (0, 0, 0)}


def sim_ls():
    return {'1V': '-', '0V': '--', 'CE': '-', 'ER': '-'}


# Read in the ICON simulation output files from the file number
# This function prepends the appropriate number of zeros
def file_prefix(j):
    if len(str(j)) == 1:
       return '000'
    elif len(str(j)) == 2:
       return '00'
    elif len(str(j)) == 3:
       return '0'
    else:
       return 'Inappropriate length of input to file_prefix'


# Read in the trajectory module output files from the file number
# This function prepends the appropriate number of zeros
def traj_prefix(j):
    if len(str(j)) == 1:
       return '00'
    elif len(str(j)) == 2:
       return '0'
    else:
       return 'Inappropriate length of input to traj_prefix'


# Helper function to - you guessed it - make sexy axes for generic values
def sexy_axes(ax,fs):
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params('both',labelsize=fs,rotation=45)


# Helper function to make sexy axes for pressure levels
def sexy_axes2(ax, fs, ylab):
    ax.set_ylim([50,800])
    ax.set_yscale('log')
    ax.set_yticks([800,500,300,100])
    ax.set_yticklabels(['800','500','300','100'])
    if ylab == True:
        ax.set_ylabel('Pressure [hPa]',fontsize=fs)
    ax.invert_yaxis()

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params('both',labelsize=fs,rotation=45)


# Define a function to convert radians to degrees.
def rad2deg(x):
    return x*180/pi


# More helper functions
# Pulling this from the following stackoverflow
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
def resize_colorbar(event):
    plt.draw()
    posn = ax.get_position()
    # left, bottom, width, height
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0 - 0.018,
                          0.025, posn.height + 0.03])
    cbar_ax.tick_params(labelsize=fs)
