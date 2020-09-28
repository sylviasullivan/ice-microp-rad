import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

basedir = '/work/bb1131/b380873/'
vgrid = xr.open_dataset(basedir + 'vgrid_icon-grid_tropic_55e115e5s40n.nc')
vct_a = vgrid['vct_a'].values
diff1 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff1[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff1[i,1] = d[0]

vgrid = xr.open_dataset(basedir + 'vgrid_icon-grid_tropic_55e115e5s40n_tropic_TTL.nc')
vct_a = vgrid['vct_a'].values
diff = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff[i,1] = d[0]

vgrid = xr.open_dataset(basedir + 'vgrid_icon-grid_tropic_55e115e5s40n_tropic_TTL2.nc')
vct_a = vgrid['vct_a'].values
diff2 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff2[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff2[i,1] = d[0]

vgrid = xr.open_dataset(basedir + 'vgrid_icon-grid_tropic_55e115e5s40n_tropic_TTL3.nc')
vct_a = vgrid['vct_a'].values
diff3 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff3[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff3[i,1] = d[0]

vgrid = xr.open_dataset(basedir + 'vgrid_max_lay_thckn_300.nc')
vct_a = vgrid['vct_a'].values
diff4 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff4[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff4[i,1] = d[0]

vgrid = xr.open_dataset('/scratch/b/b380873/tropic_run5/vgrid_DOM01.nc') #vgrid_max_lay_thckn_200.nc')
vct_a = vgrid['vct_a'].values
diff5 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff5[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff5[i,1] = d[0]

fs = 14
fig = plt.figure()
plt.plot(diff1[:,0],diff1[:,1],color='red',label='default vgrid')
plt.plot(diff[:,0],diff[:,1],color='blue',label='TTL vgrid')
plt.plot(diff2[:,0],diff2[:,1],color='green',label='TTL vgrid2')
plt.plot(diff3[:,0],diff3[:,1],color='black',label='TTL vgrid3')
plt.plot(diff4[:,0],diff4[:,1],color='red',linestyle='--',label='htop_thcknlimit = 300')
plt.plot(diff5[:,0],diff5[:,1],color='green',linestyle='--',label='htop_thcknlimit = 200')
plt.xlabel('Layer thickness [m]',fontsize=fs)
plt.ylabel('Altitude [m]',fontsize=fs)
plt.legend()
#fig.savefig('vgrid_spacings.pdf',bbox_inches='tight')
plt.show()
