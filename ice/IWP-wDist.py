import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

fs =14
basedir = '/work/bb1131/b380873/tropic_run2_output/'
TQI_fi = xr.open_dataset(basedir + 'TQI_ALL_0.025.nc')
W_fi = xr.open_dataset(basedir + 'W_ALL-0051-0061_0.025.nc')
OLR_fi = xr.open_dataset(basedir + 'OLR_120-132_0.025deg.nc')

# Extract the overlapping times between the datasets.
t1 = W_fi.time
t2 = TQI_fi.time
t3 = OLR_fi.time

# Make sure the times correspond between datasets and that the lons
# correspond to the smaller domain.
TQI = TQI_fi.tqi.where(t1 == t3)
TQI = TQI.sel(lon=slice(55,115)).values.flatten()
W250 = W_fi.w.where(t1 == t3)
W250 = W250.sel(height_2=53).values.flatten()
W500 = W_fi.w.where(t1 == t3)
W500 = W500.sel(height_2=63).values.flatten()
print(W250.min(),W250.mean(),W250.max())
print(W500.min(),W500.mean(),W500.max())

fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(13,5))
d = -3.5; u = 3.5; b = 40
wgts = np.ones_like(TQI)/float(len(TQI))*100
sns.distplot(TQI*1000,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[0],label=r'1-mom',bins=np.logspace(d,u,b),color='r',hist_kws={'weights':wgts})
del TQI

d = -5; u = 15; b = 100
wgts = np.ones_like(W250)/float(len(W250))*100
sns.distplot(W250,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[1],label=r'1-mom',bins=np.linspace(d,u,b),color='r',hist_kws={'weights':wgts})
wgts = np.ones_like(W500)/float(len(W500))*100
sns.distplot(W500,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[2],label=r'1-mom',bins=np.linspace(d,u,b),color='r',hist_kws={'weights':wgts})

# Make sure the times correspond between datasets and that the lons
# correspond to the smaller domain.
basedir = '/work/bb1131/b380873/tropic_run5_output/'
OLR_fi = xr.open_dataset(basedir + 'OLR_120-141_0.025deg.nc')
TQI_fi = xr.open_dataset(basedir + 'TQI_120-141_0.025deg.nc')
W_fi = xr.open_dataset(basedir + 'W_60-70_0.025deg_plev.nc')

# Make sure the times correspond (must be hourly given 3D w field)
OLR = OLR_fi.sel(time=W_fi.time).thb_t
TQI = TQI_fi.sel(time=W_fi.time).tqi.values.flatten()
W250 = W_fi.sel(time=W_fi.time,plev_2=25000).w.values.flatten()
W500 = W_fi.sel(time=W_fi.time,plev_2=50000).w.values.flatten()
W500 = W500[~np.isnan(W500)]
print(W250.min(),W250.mean(),W250.max())
print(W500.min(),W500.mean(),W500.max())

d = -3.5; u = 4; b = 40
wgts = np.ones_like(TQI)/float(len(TQI))*100
sns.distplot(TQI*1000,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[0],label=r'2-mom',bins=np.logspace(d,u,b),color='b',hist_kws={'weights':wgts})
del TQI

d = -5; u = 15; b = 100
wgts = np.ones_like(W250)/float(len(W250))*100
sns.distplot(W250,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[1],label=r'2-mom',bins=np.linspace(d,u,b),color='b',hist_kws={'weights':wgts})
wgts = np.ones_like(W500)/float(len(W500))*100
sns.distplot(W500,kde=False,hist=True,kde_kws={'shade':True,'linewidth':3},\
   ax=ax[2],label=r'2-mom',bins=np.linspace(d,u,b),color='b',hist_kws={'weights':wgts})


ax[0].set_xscale('log')
ax[0].set_xlim([10**(-3),10**4])
ax[0].set_xlabel(r'Ice water path [g m$^{-2}$]',fontsize=fs)
ax[0].set_ylabel('Probability [%]',fontsize=fs)
ax[0].legend()
ax[0].tick_params('both',labelsize=fs-2)
ax[0].text(0.05,0.92,'(a)',weight='bold',fontsize=fs+2,transform=ax[0].transAxes)

ax[1].set_yscale('log')
ax[1].set_xlim([-3,10])
ax[1].set_ylim([10**(-3),100])
ax[1].set_xlabel(r'$w_{250}$ [m s$^{-1}$]',fontsize=fs)
ax[1].tick_params('both',labelsize=fs-2)
ax[1].text(0.05,0.92,'(b)',weight='bold',fontsize=fs+2,transform=ax[1].transAxes)

ax[2].set_yscale('log')
ax[2].set_xlim([-3,10])
ax[2].set_ylim([10**(-3),100])
ax[2].set_xlabel(r'$w_{500}$ [m s$^{-1}$]',fontsize=fs)
ax[2].tick_params('both',labelsize=fs-2)
ax[2].text(0.05,0.92,'(c)',weight='bold',fontsize=fs+2,transform=ax[2].transAxes)

#fig.savefig('output/IWP-wDist.pdf',bbox_inches='tight')
plt.show()

