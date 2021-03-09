import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

fs =14
basedir = '/work/bb1018/b380873/model_output/ICON/'
TQI_fi = xr.open_dataset(basedir + 'TQI_ALL_0.025deg_tropic_run2.nc')
TQI = TQI_fi.tqi.sel(lon=slice(55,115)).values.flatten()
TQI = TQI[(TQI >= 10**(-7)) & (TQI <= 10)]
print(np.nanmean(TQI))

fig = plt.figure(figsize=(7,5))
d = -4; u = 4; b = 40
wgts = np.ones_like(TQI)/float(len(TQI))*100
sns.distplot(TQI*1000,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   label=r'1-mom',bins=np.logspace(d,u,b),color='r',hist_kws={'weights':wgts})

# Make sure the times correspond between datasets and that the lons
# correspond to the smaller domain.
basedir = '/work/bb1018/b380873/model_output/ICON/'
TQI_fi = xr.open_dataset(basedir + 'TQI_120-141_0.025deg_tropic_run5.nc')
TQI = TQI_fi.tqi.values.flatten()
TQI = TQI[(TQI >= 10**(-7)) & (TQI <= 10)]
print(np.nanmean(TQI))

d = -4; u = 4; b = 40
wgts = np.ones_like(TQI)/float(len(TQI))*100
sns.distplot(TQI*1000,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   label=r'2-mom',bins=np.logspace(d,u,b),color='b',hist_kws={'weights':wgts})

plt.xscale('log')
plt.xlim([10**(-4),10**4])
plt.xlabel(r'Ice water path [g m$^{-2}$]',fontsize=fs)
plt.ylabel('Probability [%]',fontsize=fs)
plt.legend()
plt.gca().tick_params('both',labelsize=fs-2)

#fig.savefig('./output/IWPdist.pdf',bbox_inches='tight')
plt.show()

