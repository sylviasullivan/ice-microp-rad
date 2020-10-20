import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid, 2mom, run6 (Atest)
arr = [False, False, False, False, False]
#arr = [True, True, True, True, True]

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
def meanProfile(basedir,tag,f,startindx,endindx,fileprefix):
    # How many vertical levels depends on which set we look at
    if suffix1 == '_PL2':
       c = 120
    elif suffix1 == '_PL':
       c = 18
    QC = np.zeros((24,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        flx = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
        clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
        QC[i-startindx] = flx.qc.isel(time=0).where(clch > f).mean(dim={'ncells'})
    print('Saving hourly mean profiles from ' + basedir + '...')
    np.save('../output/QC_' + tag + suffix2 + '.npy',QC)
    return QC

if arr[0] == True:
   QC_1mom = meanProfile('/scratch/b/b380873/tropic_run2/','1mom',0,1,24,'CLCONV_3D')
else:
   QC_1mom = np.load('../output/QC_1mom' + suffix2 + '.npy')


if arr[1] == True:
   QC_no2mom = meanProfile('/scratch/b/b380873/tropic_run5_no2mom/','no2mom',0,1,24,'Q')
else:
   QC_no2mom = np.load('../output/QC_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   QC_novgrid = meanProfile('/scratch/b/b380873/tropic_run5_novgrid/','novgrid',0,1,24,'Q')
else:
   QC_novgrid = np.load('../output/QC_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   QC_2mom = meanProfile('/scratch/b/b380873/tropic_run5/','2mom',0,49,70,'CLCONV_3D')
else:
   QC_2mom = np.load('../output/QC_2mom' + suffix2 + '.npy')

if arr[4] == True:
   QC_Atest = meanProfile('/scratch/b/b380873/tropic_run6/','Atest',0,1,23,'CLCONV_3D')
else:
   QC_Atest = np.load('../output/QC_Atest' + suffix2 + '.npy')

# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig = plt.figure()#figsize=(5.5,5.5))
plt.plot(np.nanmean(QC_1mom,axis=0)*1000,pl/100,color='red',label='ICON-1mom')
print(np.nanmean(QC_1mom,axis=0).max())
m = np.nanmean(QC_1mom,axis=0).max()
i = np.argmax(np.nanmean(QC_1mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='red')

plt.plot(np.nanmean(QC_2mom,axis=0)*1000,pl/100,color='green',label='ICON-2mom')
print(np.nanmean(QC_2mom,axis=0).max())
m = np.nanmean(QC_2mom,axis=0).max()
i = np.argmax(np.nanmean(QC_2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='green')

plt.plot(np.nanmean(QC_no2mom,axis=0)*1000,pl/100,color='blue',label='ICON-no2mom')
print(np.nanmean(QC_no2mom,axis=0).max())
m = np.nanmean(QC_no2mom,axis=0).max()
i = np.argmax(np.nanmean(QC_no2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='blue')

plt.plot(np.nanmean(QC_novgrid,axis=0)*1000,pl/100,color='gold',label='ICON-novgrid')
print(np.nanmean(QC_novgrid,axis=0).max())
m = np.nanmean(QC_novgrid,axis=0).max()
i = np.argmax(np.nanmean(QC_novgrid,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='gold')

plt.plot(np.nanmean(QC_Atest,axis=0)*1000,pl/100,color='black',label='ICON-Atest')

plt.plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Domain-mean daily-mean $q_c$ [g kg$^{-1}$]',fontsize=fs)
plt.legend()
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.tick_params(labelsize=fs)

fig.savefig('../output/qc-profiles.pdf')
plt.show()
