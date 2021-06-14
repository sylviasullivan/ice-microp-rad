import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import sys
sys.path.append('/work/bb1018/b380873/tropic_vis/')
from plotting_utilities import sim_colors, sim_ls, file_prefix, sexy_axes2

# Which set of pressure levels to look at?
suffix = '_PL2' # '_PL'

# Which simulation acronyms to calculate for?
acronym = ['0V2M1A1R']
others = ['0V1M0A0R', '1V1M0A0R', '0V2M0A0R', '1V2M0A0R', '0V2M0A1R', '1V2M0A1R',
          '0V2M1A0R', '1V2M1A0R', '1V2M1A1R']

# basedir is the base directory where the nc files are found.
# acronym is how to label the output npy.
# f is the fraction of high cloud coveraged required.
# startindx and endinx are the files over which to iterate.
def meanProfile(basedir, acronym, f, startindx, endindx, fileprefix):
    # How many vertical levels depends on which set we look at
    if suffix == '_PL2':
       c = 120
    elif suffix == '_PL':
       c = 18
    QG = np.zeros((24,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        flx = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix + '.nc')
        clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
        QG[i-startindx] = flx.qg.isel(time=0).where(clch > f).mean(dim={'ncells'})
    print('Saving hourly mean profiles from ' + basedir + '...')
    np.save('../output/QG_' + acronym + suffix + '.npy',QG)
    return QG


for a in acronym:
   QG = meanProfile('/scratch/b/b380873/' + a + '/', a, 0, 1, 25, 'CLCONV_3D')

ll = len(others + acronym)
QG_all = np.zeros((ll, 24, 120))
for indx, o in enumerate(others + acronym):
    QG_all[indx] = np.load('../output/QG_' + o + suffix + '.npy')

# Retrieve the pressure levels
pl = np.loadtxt('../remapping/PMEAN_48-72.txt')
farbe = sim_colors()
stil = sim_ls()


fs = 13
fig = plt.figure()#figsize=(5.5,5.5))
for q, o in zip(QG_all, others + acronym):
    plt.plot(np.nanmean(q,axis=0)*1000, pl/100, color=farbe[o], ls=stil[o[:2]], label=o)
    print(np.nanmean(q*1000,axis=0).max())
    m = np.nanmean(q,axis=0).max()
    i = np.argmax(np.nanmean(q,axis=0))
    plt.plot([0,m*1000],[pl[i]/100,pl[i]/100], lw=0.5, ls='--', color=farbe[o])
plt.plot([0,0], [50,800], lw=0.75, linestyle='--', color='k')
plt.ylabel('Pressure [hPa]', fontsize=fs)
plt.xlabel(r'Domain-mean daily-mean $q_g$ [g kg$^{-1}$]', fontsize=fs)
plt.legend()
sexy_axes2(plt.gca(), fs, True)

fig.savefig('../output/qg-profiles_all.pdf')
plt.show()
