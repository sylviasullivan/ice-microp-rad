import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid, 2mom
#arr = [False, False, False, False]
arr = [False, False, False, False]

def file_prefix(j):
    if len(str(j)) == 1:
       return '000'
    elif len(str(j)) == 2:
       return '00'
    elif len(str(j)) == 3:
       return '0'
    else:
       return 'Inappropriate length of input to file_prefix'

# basedir is the base directory where the nc files are found.
# tag is how to label the output npy.
# f is the fraction of high cloud coveraged required.
# startindx and endinx are the files over which to iterate.
def meanProfile(basedir, tag, f, startindx, endindx, fileprefix):
    # How many vertical levels depends on which set we look at
    if suffix1 == '_PL2':
       c = 120
    elif suffix1 == '_PL':
       c = 18
    WT = np.zeros((2,24,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        w = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
        clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
        WT[0,i-startindx] = w.omega.isel(time=0).where(clch > f). mean(dim={'ncells'})
        WT[1,i-startindx] = w.temp.isel(time=0).where(clch > f). mean(dim={'ncells'})

    print('Saving mean pressure velocity and temperature from ' + basedir + '...')
    np.save('../output/WT_' + tag + '_' + suffix2 + '.npy',WT)
    return WT

if arr[0] == True:
   WT_1mom = meanProfile('/scratch/b/b380873/tropic_run2/','1mom',0,1,24,'WT')
else:
   WT_1mom = np.load('../output/WT_1mom' + suffix2 + '.npy')


if arr[1] == True:
   WT_no2mom = meanProfile('/scratch/b/b380873/tropic_run5_no2mom/','no2mom',0,1,24,'WT')
else:
   WT_no2mom = np.load('../output/WT_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   #WT_novgrid = meanProfile('/scratch/b/b380873/tropic_run5_novgrid/','novgrid',0,1,24,'WT')
   WT_radnovgrid = meanProfile('/scratch/b/b380873/tropic_run7_radnovgrid/','radnovgrid',0,1,24,'WT')
else:
   WT_novgrid = np.load('../output/WT_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   #WT_2mom = meanProfile('/scratch/b/b380873/tropic_run5/','2mom',0,48,72,'WT')
   WT_rad2mom = meanProfile('/scratch/b/b380873/tropic_run7_rad2mom/','rad2mom',0,1,24,'WT')
else:
   WT_2mom = np.load('../output/WT_2mom' + suffix2 + '.npy')


# Retrieve the pressure levels
pl = np.loadtxt('../remapping/PMEAN_48-72.txt')

fs = 13
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(9,6.5))
p = 0.1
#ax[0].plot(np.percentile(WT_1mom[:,0],1,axis=0),pl/100,color='red',label='ICON-1mom')
ax[0].plot((np.percentile(WT_2mom[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='green',label='ICON-2mom')
ax[0].plot((np.percentile(WT_no2mom[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='blue',label='ICON-no2mom')
ax[0].plot((np.percentile(WT_novgrid[:,0],p,axis=0)-np.percentile(WT_1mom[:,0],p,axis=0))*100,\
       pl/100,color='gold',label='ICON-novgrid')
ax[0].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
ax[0].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[0].set_xlabel(r'Domain-mean $\Delta\omega_{99.9}$ '
                  '\n'
                 'from 1-mom [hPa s$^{-1}$]',fontsize=fs)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+4,transform=ax[0].transAxes)
ax[0].legend(loc='upper right')

#ax1 = fig.add_axes([0.33,0.55,0.07,0.25])
#ax1.plot(np.percentile(WT_1mom[:,0],p,axis=0)*100,pl/100,color='red',label='_nolegend_')
#ax1.set_xlabel('$\omega_{99.9}$',fontsize=fs-3)
#ax1.set_yscale('log')
#ax1.set_ylim([50,800])
#ax1.set_yticks([800,500,300,100])
#ax1.set_yticklabels(['800','500','300','100'])
#ax1.invert_yaxis()

ax[0].set_ylim([50,800])
ax[0].set_yscale('log')
ax[0].invert_yaxis()
ax[0].set_yticks([800,500,300,100])
ax[0].set_yticklabels(['800','500','300','100'])
ax[0].tick_params(labelsize=fs)

#ax[1].plot(np.nanmean(WT_1mom[:,0],axis=0),pl/100,color='red',label='ICON-1mom')
ax[1].plot((np.nanmean(WT_2mom[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='green',label='ICON-2mom')
ax[1].plot((np.nanmean(WT_no2mom[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='blue',label='ICON-no2mom')
ax[1].plot((np.nanmean(WT_novgrid[:,0],axis=0)-np.nanmean(WT_1mom[:,0],axis=0))*100,pl/100,\
       color='gold',label='ICON-novgrid')
ax[1].plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
#ax[1].set_ylabel('Pressure [hPa]',fontsize=fs)
ax[1].set_xlabel(r'Domain-mean daily-mean $\Delta\omega$'
                 '\n'
                 'from 1-mom [hPa s$^{-1}$]',fontsize=fs)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+4,transform=ax[1].transAxes)

ax2 = fig.add_axes([0.8,0.55,0.07,0.25])
ax2.plot(np.nanmean(WT_1mom[:,0],axis=0)*100,pl/100,color='red',label='_nolegend_')
ax2.set_xlabel('$<\omega>$',fontsize=fs-3)
ax2.set_yscale('log')
ax2.set_ylim([50,800])
ax2.set_yticks([800,500,300,100])
ax2.set_yticklabels(['800','500','300','100'])
ax2.invert_yaxis()

ax[1].set_ylim([50,800])
ax[1].set_yscale('log')
ax[1].invert_yaxis()
ax[1].set_yticks([800,500,300,100])
ax[1].set_yticklabels(['800','500','300','100'])
ax[1].tick_params(labelsize=fs)

#fig.savefig('../output/wT-profiles2.pdf')
plt.show()
