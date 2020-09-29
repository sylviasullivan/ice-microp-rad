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

# OLR > 240, below increasingly stringent criteria for OLR
tqi_fi = xr.open_dataset(basedir + 'tqi_olrgt240_2017080800-2017080806.nc')
tqi_1_all_vals = tqi_fi.tqi.values
tqi_1_all_vals = np.reshape(tqi_1_all_vals,(13*1800*4600))*1000
tqi_1_real_vals = tqi_1_all_vals[~np.isnan(tqi_1_all_vals)]
tqi_1_real_vals = tqi_1_real_vals[tqi_1_real_vals > 10**(-5)]
print(tqi_1_real_vals.shape)
print(np.nanmin(tqi_1_real_vals),np.nanmedian(tqi_1_real_vals),np.nanmean(tqi_1_real_vals),np.nanmax(tqi_1_real_vals))

# OLR > 220
tqi_fi = xr.open_dataset(basedir + 'tqi_olrgt220_2017080800-2017080806.nc')
tqi_2_all_vals = tqi_fi.tqi.values
tqi_2_all_vals = np.reshape(tqi_2_all_vals,(13*1800*4600))*1000
tqi_2_real_vals = tqi_2_all_vals[~np.isnan(tqi_2_all_vals)]
tqi_2_real_vals = tqi_2_real_vals[tqi_2_real_vals > 10**(-5)]
print(tqi_2_real_vals.shape)
print(np.nanmin(tqi_2_real_vals),np.nanmedian(tqi_2_real_vals),np.nanmean(tqi_2_real_vals),np.nanmax(tqi_2_real_vals))

# OLR > 200
tqi_fi = xr.open_dataset(basedir + 'tqi_olrgt200_2017080800-2017080806.nc')
tqi_3_all_vals = tqi_fi.tqi.values
tqi_3_all_vals = np.reshape(tqi_3_all_vals,(13*1800*4600))*1000
tqi_3_real_vals = tqi_3_all_vals[~np.isnan(tqi_3_all_vals)]
tqi_3_real_vals = tqi_3_real_vals[tqi_3_real_vals > 10**(-5)]
print(tqi_3_real_vals.shape)
print(np.nanmin(tqi_3_real_vals),np.nanmedian(tqi_3_real_vals),np.nanmean(tqi_3_real_vals),np.nanmax(tqi_3_real_vals))

# OLR > 180
tqi_fi = xr.open_dataset(basedir + 'tqi_olrgt180_2017080800-2017080806.nc')
tqi_4_all_vals = tqi_fi.tqi.values   # 13 time steps, 1800 lat, 4600 lon
tqi_4_all_vals = np.reshape(tqi_4_all_vals,(13*1800*4600))*1000
tqi_4_real_vals = tqi_4_all_vals[~np.isnan(tqi_4_all_vals)]
tqi_4_real_vals = tqi_4_real_vals[tqi_4_real_vals > 10**(-5)]
print(tqi_4_real_vals.shape)
print(np.nanmin(tqi_4_real_vals),np.nanmedian(tqi_4_real_vals),np.nanmean(tqi_4_real_vals),np.nanmax(tqi_4_real_vals))

# OLR > 160
tqi_fi = xr.open_dataset(basedir + 'tqi_olrgt160_2017080800-2017080806.nc')
tqi_5_all_vals = tqi_fi.tqi.values
tqi_5_all_vals = np.reshape(tqi_5_all_vals,(13*1800*4600))*1000
tqi_5_real_vals = tqi_5_all_vals[~np.isnan(tqi_5_all_vals)]
tqi_5_real_vals = tqi_5_real_vals[tqi_5_real_vals > 10**(-5)]
print(tqi_5_real_vals.shape)
print(np.nanmin(tqi_5_real_vals),np.nanmedian(tqi_5_real_vals),np.nanmean(tqi_5_real_vals),np.nanmax(tqi_5_real_vals))

fs = 13
fig = plt.figure(figsize=(6,5))
d = -3.5; u = 4; n = 30
#wgts = np.ones_like(tqi_all_vals)/len(tqi_all_vals)*100
h1, bin_edges = np.histogram(tqi_all_vals,bins=np.logspace(d,u,n)) #,weights=wgts
#wgts = np.ones_like(tqi_1_all_vals)/len(tqi_1_all_vals)*100
h2, _ = np.histogram(tqi_1_real_vals,bins=np.logspace(d,u,n)) #,weights=wgts
h3, _ = np.histogram(tqi_2_real_vals,bins=np.logspace(d,u,n))
h4, _ = np.histogram(tqi_3_real_vals,bins=np.logspace(d,u,n))
h5, _ = np.histogram(tqi_4_real_vals,bins=np.logspace(d,u,n))
h6, _ = np.histogram(tqi_5_real_vals,bins=np.logspace(d,u,n))

#bins=np.logspace(np.floor(tqi_all_vals.min()),np.ceil(tqi_all_vals.max()))
#h1, bin_edges = np.histogram(tqi_all_vals*1000,density=True,bins=bins) #np.logspace(-3.5,3.5,30))
#h2, _ = np.histogram(tqi_1_all_vals*1000,density=True,bins=bins) #np.logspace(-3.5,3.5,30))


bin_center = np.zeros((bin_edges.shape[0]-1,))
for i,b in enumerate(bin_edges):
    if i != len(bin_edges)-1:
       bin_center[i] = (b + bin_edges[i+1])/2
    else:
       break

print(bin_center,h1,len(tqi_all_vals))
#plt.step(bin_center,h1/len(tqi_all_vals)*100,color='r',label='All cloud',linewidth=0.75)
#plt.plot([np.nanmean(tqi_all_vals),np.nanmean(tqi_all_vals)],[0,8.5],color='r',linewidth=1.25)
plt.step(bin_center,h2/float(len(tqi_1_real_vals))*100.,color='red',\
         label=r'|OLR| $\leq$ 240 W m$^{-2}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_1_real_vals),np.nanmean(tqi_1_real_vals)],[0,8.5],color='red',linewidth=1.25)
plt.step(bin_center,h3/float(len(tqi_2_real_vals))*100.,color='gold',\
         label=r'|OLR| $\leq$ 220 W m$^{-2}$',linewidth=0.75)  # 2.5
plt.plot([np.nanmean(tqi_2_real_vals),np.nanmean(tqi_2_real_vals)],[0,8.5],color='gold',linewidth=1.25)
plt.step(bin_center,h4/float(len(tqi_3_real_vals))*100.,color='green',\
         label=r'|OLR| $\leq$ 200 W m$^{-2}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_3_real_vals),np.nanmean(tqi_3_real_vals)],[0,8.5],color='green',linewidth=1.25)
plt.step(bin_center,h5/float(len(tqi_4_real_vals))*100.,color='blue',\
         label=r'|OLR| $\leq$ 180 W m$^{-2}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_4_real_vals),np.nanmean(tqi_4_real_vals)],[0,8.5],color='blue',linewidth=1.25)
plt.step(bin_center,h6/float(len(tqi_5_real_vals))*100.,color='purple',\
         label=r'|OLR| $\leq$ 160 W m$^{-2}$',linewidth=0.75)
plt.plot([np.nanmean(tqi_5_real_vals),np.nanmean(tqi_5_real_vals)],[0,8.5],color='purple',linewidth=1.25)
#plt.title('(a) 2017080800-2017080806 IWP dist')
plt.ylim([0,12]); plt.ylabel('Probability [%]',fontsize=fs)
plt.xlim([10**(d),10**(3.)]); plt.xlabel(r'IWP [kg m$^{-2}$]',fontsize=fs)
plt.gca().set_xscale('log')
plt.legend(loc='upper left',fontsize=10)
fig.savefig('IWPdist-log-olr.pdf',bbox_inches='tight')
plt.show()
