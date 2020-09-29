import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

# I'm not sure there is a cleaner way to do this. At least not offhand as each list has a different length.
# All cloud
basedir = '/work/bb1131/b380873/tropic_run2_output/'
tqi_fi = xr.open_dataset(basedir + 'tqi_2017080800-2017080806.nc')
tqi_all_vals = tqi_fi.tqi.values
tqi_all_vals = np.reshape(tqi_all_vals,(12*1800*4600))*1000
tqi_all_vals = tqi_all_vals[tqi_all_vals > 10**(-5)]
print(tqi_all_vals.shape)
print(np.nanmin(tqi_all_vals),np.nanmedian(tqi_all_vals),np.nanmean(tqi_all_vals),np.nanmax(tqi_all_vals))

# P > 0
tqi_fi = xr.open_dataset(basedir + 'tqi_pgt0_2017080800-2017080806.nc')
tqi_pgt0_all_vals = tqi_fi.tqi.values
tqi_pgt0_all_vals = np.reshape(tqi_pgt0_all_vals,(13*1800*4600))*1000
tqi_pgt0_real_vals = tqi_pgt0_all_vals[~np.isnan(tqi_pgt0_all_vals)]
tqi_pgt0_real_vals = tqi_pgt0_real_vals[tqi_pgt0_real_vals > 10**(-5)]
print(tqi_pgt0_real_vals.shape)
print(np.nanmin(tqi_pgt0_real_vals),np.nanmedian(tqi_pgt0_real_vals),np.nanmean(tqi_pgt0_real_vals),np.nanmax(tqi_pgt0_real_vals))

# P > 1
tqi_fi = xr.open_dataset(basedir + 'tqi_pgt1_2017080800-2017080806.nc')
tqi_pgt25_all_vals = tqi_fi.tqi.values
tqi_pgt25_all_vals = np.reshape(tqi_pgt25_all_vals,(13*1800*4600))*1000
tqi_pgt25_real_vals = tqi_pgt25_all_vals[~np.isnan(tqi_pgt25_all_vals)]
tqi_pgt25_real_vals = tqi_pgt25_real_vals[tqi_pgt25_real_vals > 10**(-5)]
print(tqi_pgt25_real_vals.shape)
print(np.nanmin(tqi_pgt25_real_vals),np.nanmedian(tqi_pgt25_real_vals),np.nanmean(tqi_pgt25_real_vals),np.nanmax(tqi_pgt25_real_vals))

# P > 5
tqi_fi = xr.open_dataset(basedir + 'tqi_pgt5_2017080800-2017080806.nc')
tqi_pgt5_all_vals = tqi_fi.tqi.values
tqi_pgt5_all_vals = np.reshape(tqi_pgt5_all_vals,(13*1800*4600))*1000
tqi_pgt5_real_vals = tqi_pgt5_all_vals[~np.isnan(tqi_pgt5_all_vals)]
tqi_pgt5_real_vals = tqi_pgt5_real_vals[tqi_pgt5_real_vals > 10**(-5)]
print(tqi_pgt5_real_vals.shape)
print(np.nanmin(tqi_pgt5_real_vals),np.nanmedian(tqi_pgt5_real_vals),np.nanmean(tqi_pgt5_real_vals),np.nanmax(tqi_pgt5_real_vals))

# P > 10
tqi_fi = xr.open_dataset(basedir + 'tqi_pgt10_2017080800-2017080806.nc')
tqi_pgt10_all_vals = tqi_fi.tqi.values   # 13 time steps, 1800 lat, 4600 lon
tqi_pgt10_all_vals = np.reshape(tqi_pgt10_all_vals,(13*1800*4600))*1000
tqi_pgt10_real_vals = tqi_pgt10_all_vals[~np.isnan(tqi_pgt10_all_vals)]
tqi_pgt10_real_vals = tqi_pgt10_real_vals[tqi_pgt10_real_vals > 10**(-5)]
print(tqi_pgt10_real_vals.shape)
print(np.nanmin(tqi_pgt10_real_vals),np.nanmedian(tqi_pgt10_real_vals),np.nanmean(tqi_pgt10_real_vals),np.nanmax(tqi_pgt10_real_vals))

# P > 20
tqi_fi = xr.open_dataset(basedir + 'tqi_pgt20_2017080800-2017080806.nc')
tqi_pgt20_all_vals = tqi_fi.tqi.values
tqi_pgt20_all_vals = np.reshape(tqi_pgt20_all_vals,(13*1800*4600))*1000
tqi_pgt20_real_vals = tqi_pgt20_all_vals[~np.isnan(tqi_pgt20_all_vals)]
tqi_pgt20_real_vals = tqi_pgt20_real_vals[tqi_pgt20_real_vals > 10**(-5)]
print(tqi_pgt20_real_vals.shape)
print(np.nanmin(tqi_pgt20_real_vals),np.nanmedian(tqi_pgt20_real_vals),np.nanmean(tqi_pgt20_real_vals),np.nanmax(tqi_pgt20_real_vals))

fs = 13
fig = plt.figure()
d = -3.5; u = 4; n = 30
#wgts = np.ones_like(tqi_all_vals)/len(tqi_all_vals)*100
h1, bin_edges = np.histogram(tqi_all_vals,bins=np.logspace(d,u,n)) #,weights=wgts
#wgts = np.ones_like(tqi_pgt0_all_vals)/len(tqi_pgt0_all_vals)*100
h2, _ = np.histogram(tqi_pgt0_real_vals,bins=np.logspace(d,u,n)) #,weights=wgts
h3, _ = np.histogram(tqi_pgt25_real_vals,bins=np.logspace(d,u,n))
h4, _ = np.histogram(tqi_pgt5_real_vals,bins=np.logspace(d,u,n))
h5, _ = np.histogram(tqi_pgt10_real_vals,bins=np.logspace(d,u,n))
h6, _ = np.histogram(tqi_pgt20_real_vals,bins=np.logspace(d,u,n))

#bins=np.logspace(np.floor(tqi_all_vals.min()),np.ceil(tqi_all_vals.max()))
#h1, bin_edges = np.histogram(tqi_all_vals*1000,density=True,bins=bins) #np.logspace(-3.5,3.5,30))
#h2, _ = np.histogram(tqi_pgt0_all_vals*1000,density=True,bins=bins) #np.logspace(-3.5,3.5,30))


bin_center = np.zeros((bin_edges.shape[0]-1,))
for i,b in enumerate(bin_edges):
    if i != len(bin_edges)-1:
       bin_center[i] = (b + bin_edges[i+1])/2
    else:
       break

plt.step(bin_center,h1/len(tqi_all_vals)*100,color='r',label='All cloud',linewidth=0.75)
plt.plot([np.nanmean(tqi_all_vals),np.nanmean(tqi_all_vals)],[0,8.5],color='r',linewidth=1.25)
plt.step(bin_center,h2/len(tqi_pgt0_real_vals)*100,color='orange',label=r'$P >$ 0 mm h$^{-1}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_pgt0_real_vals),np.nanmean(tqi_pgt0_real_vals)],[0,8.5],color='orange',linewidth=1.25)
plt.step(bin_center,h3/len(tqi_pgt25_real_vals)*100,color='gold',label=r'$P >$ 1 mm h$^{-1}$',linewidth=0.75)  # 2.5
plt.plot([np.nanmean(tqi_pgt25_real_vals),np.nanmean(tqi_pgt25_real_vals)],[0,8.5],color='gold',linewidth=1.25)
plt.step(bin_center,h4/len(tqi_pgt5_real_vals)*100,color='green',label=r'$P >$ 5 mm h$^{-1}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_pgt5_real_vals),np.nanmean(tqi_pgt5_real_vals)],[0,8.5],color='green',linewidth=1.25)
plt.step(bin_center,h5/len(tqi_pgt10_real_vals)*100,color='blue',label=r'$P >$ 10 mm h$^{-1}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_pgt10_real_vals),np.nanmean(tqi_pgt10_real_vals)],[0,8.5],color='blue',linewidth=1.25)
plt.step(bin_center,h6/len(tqi_pgt20_real_vals)*100,color='purple',label=r'$P >$ 20 mm h$^{-1}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_pgt20_real_vals),np.nanmean(tqi_pgt20_real_vals)],[0,8.5],color='purple',linewidth=1.25)
plt.legend()
#plt.title('(a) 2017080800-2017080806 IWP dist')
plt.ylim([0,9]); plt.ylabel('Probability [%]',fontsize=fs)
plt.xlim([10**(d),10**(3.5)]); plt.xlabel(r'IWP [kg m$^{-2}$]',fontsize=fs)
plt.gca().set_xscale('log')
#fig.savefig('IWPdist-log-precip2.pdf',bbox_inches='tight')
plt.show()
