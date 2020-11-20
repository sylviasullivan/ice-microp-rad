import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid / radnovgrid, 2mom / rad2mom, run6 (Atest) / PDA
arr = [False, False, False, False, False]

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
    QG = np.zeros((24,c))

    for i in np.arange(startindx,endindx):
        print(i)
        prefix = file_prefix(i)
        flx = xr.open_dataset(basedir + fileprefix + '_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
        clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
        QG[i-startindx] = flx.qg.isel(time=0).where(clch > f).mean(dim={'ncells'})
    print('Saving hourly mean profiles from ' + basedir + '...')
    np.save('../output/QG_' + tag + suffix2 + '.npy',QG)
    return QG

if arr[0] == True:
   QG_1mom = meanProfile('/scratch/b/b380873/tropic_run2/','1mom',0,1,24,'CLCONV_3D')
else:
   QG_1mom = np.load('../output/QG_1mom' + suffix2 + '.npy')


if arr[1] == True:
   QG_no2mom = meanProfile('/scratch/b/b380873/tropic_run5_no2mom/','no2mom',0,1,24,'Q')
else:
   QG_no2mom = np.load('../output/QG_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   QG_novgrid = meanProfile('/scratch/b/b380873/tropic_run5_novgrid/','novgrid',0,1,24,'Q')
   QG_radnovgrid = meanProfile('/scratch/b/b380873/tropic_run7_radnovgrid/','radnovgrid',0,1,24,'CLCONV_3D')
else:
   QG_novgrid = np.load('../output/QG_novgrid' + suffix2 + '.npy')
   QG_radnovgrid = np.load('../output/QG_radnovgrid' + suffix2 + '.npy')


if arr[3] == True:
   QG_2mom = meanProfile('/scratch/b/b380873/tropic_run5/','2mom',0,49,73,'CLCONV_3D')
   QG_rad2mom = meanProfile('/scratch/b/b380873/tropic_run7_rad2mom/','rad2mom',0,1,24,'CLCONV_3D')
else:
   QG_2mom = np.load('../output/QG_2mom' + suffix2 + '.npy')
   QG_rad2mom = np.load('../output/QG_rad2mom' + suffix2 + '.npy')

if arr[4] == True:
   QG_Atest = meanProfile('/scratch/b/b380873/tropic_run6/','Atest',0,1,23,'CLCONV_3D')
   QG_PDA = meanProfile('/scratch/b/b380873/tropic_run8_pda/','PDA',0,1,24,'CLCONV_3D')
else:
   QG_Atest = np.load('../output/QG_Atest' + suffix2 + '.npy')
   QG_PDA = np.load('../output/QG_PDA' + suffix2 + '.npy')

# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2


fs = 13
fig = plt.figure()#figsize=(5.5,5.5))
plt.plot(np.nanmean(QG_1mom,axis=0)*1000,pl/100,color='red',label='ICON-1mom')
print(np.nanmean(QG_1mom,axis=0).max())
m = np.nanmean(QG_1mom,axis=0).max()
i = np.argmax(np.nanmean(QG_1mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='red')

plt.plot(np.nanmean(QG_2mom,axis=0)*1000,pl/100,color='green',label='ICON-2mom')
print(np.nanmean(QG_2mom,axis=0).max())
m = np.nanmean(QG_2mom,axis=0).max()
i = np.argmax(np.nanmean(QG_2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='green')

plt.plot(np.nanmean(QG_no2mom,axis=0)*1000,pl/100,color='blue',label='ICON-no2mom')
print(np.nanmean(QG_no2mom,axis=0).max())
m = np.nanmean(QG_no2mom,axis=0).max()
i = np.argmax(np.nanmean(QG_no2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='blue')

plt.plot(np.nanmean(QG_novgrid,axis=0)*1000,pl/100,color='gold',label='ICON-novgrid')
print(np.nanmean(QG_novgrid,axis=0).max())
m = np.nanmean(QG_novgrid,axis=0).max()
i = np.argmax(np.nanmean(QG_novgrid,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='gold')

plt.plot(np.nanmean(QG_Atest,axis=0)*1000,pl/100,color='black',label='ICON-Atest')
plt.plot(np.nanmean(QG_rad2mom,axis=0)*1000,pl/100,color='purple',label='ICON-rad2mom')
plt.plot(np.nanmean(QG_Atest,axis=0)*1000,pl/100,color='purple',linestyle='--',label='ICON-radnovgrid')
plt.plot(np.nanmean(QG_rad2mom,axis=0)*1000,pl/100,color='pink',label='ICON-PDA')

plt.plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Domain-mean daily-mean $q_g$ [g kg$^{-1}$]',fontsize=fs)
plt.legend()
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.tick_params(labelsize=fs)

#fig.savefig('../output/qg-profiles_radsim.pdf')
plt.show()
