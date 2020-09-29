import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# Which set of pressure levels to look at?
suffix1 = '_PL2' # '_PL'
suffix2 = '_PL2'
# Which simulations to look at? Order is 1mom, no2mom, novgrid, 2mom
arr = [False, False, False, False]
#arr = [True, True, True, True]

def file_prefix(j):
    if len(str(j)) == 1:
       return '000'
    elif len(str(j)) == 2:
       return '00'
    elif len(str(j)) == 3:
       return '0'
    else:
       return 'Inappropriate length of input to file_prefix'

# How many vertical levels depends on which set we look at
if suffix1 == '_PL2':
   c = 120
elif suffix1 == '_PL':
   c = 18

QI_1mom = np.zeros((24,c))
QI_2mom = np.zeros((21,c))
QI_no2mom = np.zeros((24,c))
QI_novgrid = np.zeros((24,c))

# Which fraction of high cloud coverage are you requiring?
f = 0

if arr[0] == True:
   for i in np.arange(1,24):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run2/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'QI_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       QI_1mom[i-1] = flx.qi.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 1mom fluxes...')
   np.save('./output/QI_1mom' + suffix2 + '.npy',QI_1mom)
else:
   QI_1mom = np.load('./output/QI_1mom' + suffix2 + '.npy')


if arr[1] == True:
   for i in np.arange(1,24):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_no2mom/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'QI_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       QI_no2mom[i-1] = flx.qi.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving no2mom fluxes...')
   np.save('./output/QI_no2mom' + suffix2 + '.npy',QI_no2mom)
else:
   QI_no2mom = np.load('./output/QI_no2mom' + suffix2 + '.npy')


if arr[2] == True:
   for i in np.arange(1,24):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5_novgrid/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'QI_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       QI_novgrid[i-1] = flx.qi.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving novgrid fluxes...')
   np.save('./output/QI_novgrid' + suffix2 + '.npy',QI_novgrid)
else:
   QI_novgrid = np.load('./output/QI_novgrid' + suffix2 + '.npy')


if arr[3] == True:
   for i in np.arange(49,70):
       print(i)
       basedir = '/scratch/b/b380873/tropic_run5/'
       prefix = file_prefix(i)
       flx = xr.open_dataset(basedir + 'QI_icon_tropic_' + prefix + str(i) + suffix1 + '.nc')
       clch = xr.open_dataset(basedir + 'CLCONV_2D_icon_tropic_' + prefix + str(i) + '.nc').clch.isel(time=0,lev=0)
       QI_2mom[i-49] = flx.qi.isel(time=0).where(clch > f).mean(dim={'ncells'})
   print('Saving 2mom fluxes...')
   np.save('./output/QI_2mom' + suffix2 + '.npy',QI_2mom)
else:
   QI_2mom = np.load('./output/QI_2mom' + suffix2 + '.npy')


# Retrieve the pressure levels
basedir = '/scratch/b/b380873/tropic_run5/'
flx = xr.open_dataset(basedir + 'FLX_icon_tropic_0049' + suffix1 + '.nc')
pl = flx.plev_2

fs = 13
fig = plt.figure(figsize=(5.5,5.5))
plt.plot(np.nanmean(QI_1mom,axis=0)*1000,pl/100,color='red',label='ICON-1mom')
print(np.nanmean(QI_1mom,axis=0).max())
m = np.nanmean(QI_1mom,axis=0).max()
i = np.argmax(np.nanmean(QI_1mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='red')

plt.plot(np.nanmean(QI_2mom,axis=0)*1000,pl/100,color='green',label='ICON-2mom')
print(np.nanmean(QI_2mom,axis=0).max())
m = np.nanmean(QI_2mom,axis=0).max()
i = np.argmax(np.nanmean(QI_2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='green')

plt.plot(np.nanmean(QI_no2mom,axis=0)*1000,pl/100,color='blue',label='ICON-no2mom')
print(np.nanmean(QI_no2mom,axis=0).max())
m = np.nanmean(QI_no2mom,axis=0).max()
i = np.argmax(np.nanmean(QI_no2mom,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='blue')

plt.plot(np.nanmean(QI_novgrid,axis=0)*1000,pl/100,color='gold',label='ICON-novgrid')
print(np.nanmean(QI_novgrid,axis=0).max())
m = np.nanmean(QI_novgrid,axis=0).max()
i = np.argmax(np.nanmean(QI_novgrid,axis=0))
plt.plot([0,m*1000],[pl[i]/100,pl[i]/100],lw=0.5,ls='--',color='gold')

plt.plot([0,0],[50,800],lw=0.75,linestyle='--',color='k')
plt.ylabel('Pressure [hPa]',fontsize=fs)
plt.xlabel(r'Domain-mean daily-mean $q_i$ [g kg$^{-1}$]',fontsize=fs)
plt.legend()
plt.ylim([50,800])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.gca().set_yticks([800,500,300,100])
plt.gca().set_yticklabels(['800','500','300','100'])
plt.tick_params(labelsize=fs)

fig.savefig('./output/qi-profiles2.pdf')
plt.show()
