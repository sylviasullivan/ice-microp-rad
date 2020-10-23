import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

basedir = '/work/bb1018/b380873/'
vgrid = xr.open_dataset(basedir + 'vgrid_icon-grid_tropic_55e115e5s40n.nc')
vct_a0 = vgrid['vct_a'].values
diff1 = np.zeros((vct_a0.shape[0]-1,2))
for i in np.arange(vct_a0.shape[0]-1):
    d = vct_a0[i] - vct_a0[i+1]
    diff1[i,0] = d[0]
    d = (vct_a0[i] + vct_a0[i+1])/2
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

vgrid = xr.open_dataset(basedir + 'vgrid_max_lay_thckn_200.nc')
vct_a = vgrid['vct_a'].values
diff5 = np.zeros((vct_a.shape[0]-1,2))
for i in np.arange(vct_a.shape[0]-1):
    d = vct_a[i] - vct_a[i+1]
    diff5[i,0] = d[0]
    d = (vct_a[i] + vct_a[i+1])/2
    diff5[i,1] = d[0]

fs = 14
fig = plt.figure(figsize=(5.5,4.5))
#plt.subplot(1,2,1)
#plt.plot(diff1[:,0],diff1[:,1],color='red',label='default vgrid')
#plt.plot(diff[:,0],diff[:,1],color='blue',label='TTL vgrid')
#plt.plot(diff2[:,0],diff2[:,1],color='green',label='TTL vgrid2')
#plt.plot(diff3[:,0],diff3[:,1],color='black',label='TTL vgrid3')
#plt.plot(diff4[:,0],diff4[:,1],color='red',linestyle='--',label='htop_thcknlimit = 300')
#plt.plot(diff5[:,0],diff5[:,1],color='green',linestyle='--',label='htop_thcknlimit = 200')
#plt.xlabel('Layer thickness [m]',fontsize=fs)
plt.ylabel('Altitude [km]',fontsize=fs)
#plt.legend()

#plt.subplot(1,2,2)
plt.plot(np.arange(120,-1,-1),vct_a[:,0]/1000,'gx',label='higher res - upper troposphere')
plt.plot(np.arange(75,-1,-1),vct_a0[:,0]/1000,'rx',label='default')
plt.xlabel('Level count',fontsize=fs)
plt.ylim(0,21)
plt.legend()
#fig.savefig('tropic_vgrid.pdf',bbox_inches='tight')
plt.show()
