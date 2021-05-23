# Define the colors and linestyles to use for different simulations
def sim_colors():
    return {'V1M0A0R': 'red', 'V2M0A0R': 'blue', 'V2M0A1R': 'goldenrod',
            'V2M1A0R': 'green', 'V2M1A1R': 'purple'}
def sim_ls():
    return {'1V': '-', '0V': '--'}

# Read in the simulation output files from the file number
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
