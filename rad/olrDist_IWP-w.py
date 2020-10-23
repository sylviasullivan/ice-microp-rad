import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
import sys

basedir = '/work/bb1018/b380873/tropic_run2_output/'
OLR_fi = xr.open_dataset(basedir + 'OLR_120-132_0.025deg.nc')
TQI_fi = xr.open_dataset(basedir + 'TQI_ALL_0.025.nc')
W_fi = xr.open_dataset(basedir + 'W_ALL-0051-0061_0.025.nc')

# Extract the overlapping times between the datasets.
t1 = OLR_fi.time
t2 = W_fi.time
t3 = TQI_fi.time
print(t1)
print(t2)
print(t3)
sys.exit()

# Make sure the times correspond between datasets and that the lons
# correspond to the smaller domain.
OLR = OLR_fi.lwflxall.where(t1 == t2)
OLR = OLR.sel(lon=slice(55,115))
TQI = TQI_fi.tqi.where(t2 == t3)
TQI = TQI.sel(lon=slice(55,115))
W250 = W_fi.w.where(t1 == t2)
W250 = W250.sel(height_2=53)
W500 = W_fi.w.where(t1 == t2)
W500 = W500.sel(height_2=64)
print(W250.min(),W250.max(),W500.min(),W500.max())

# Filter OLR into ranges of ice water path (TQI).
OLR_1 = OLR.where((TQI > 0) & (TQI < 10**(-4)),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((TQI >= 10**(-4)) & (TQI < 10**(-3)),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((TQI >= 10**(-3)) & (TQI < 10**(-2)),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((TQI >= 10**(-2)) & (TQI < 10**(-1)),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((TQI >= 10**(-1)),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

fig,ax = plt.subplots(nrows=2,ncols=3,figsize=(12,8))
d = 80; u = 360; b = 50; fs = 13
farbe = ['red','orange','yellow','green','blue','purple']
sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},label=r'Clear sky',ax=ax[0,0])
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,0],label=r'IWP $\in$ [10$^{-4}$,10$^{-3}$] g m$^{-2}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,0],label=r'IWP $\in$ [10$^{-3}$,10$^{-2}$] g m$^{-2}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,0],label=r'IWP $\in$ [10$^{-2}$,10$^{-1}$] g m$^{-2}$')
sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,0],label=r'IWP $>$ 10$^{-1}$ g m$^{-2}$')
ax[0,0].set_xlim([d-20,u])
ax[0,0].set_ylabel('Density',fontsize=fs)
ax[0,0].legend(loc='upper left',fontsize=8)

# Filter the OLR into ranges of vertical velocity at 250 hPa (W250).
OLR_1 = OLR.where((W250 > 0) & (W250 < 0.1),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((W250 >= 0.1) & (W250 < 1),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((W250 >= 1) & (W250 < 5),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((W250 >= 5),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((W250 < 0),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,1],label=r'$w_{250}$ $<$ 0 m s$^{-1}$')
sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},label=r'$w_{250}$ $\in$ [0,0.1] m s$^{-1}$',ax=ax[0,1])
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,1],label=r'$w_{250}$ $\in$ [0.1,1] m s$^{-1}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,1],label=r'$w_{250}$ $\in$ [1,5] m s$^{-1}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,1],label=r'$w_{250}$ $>$ 5 m s$^{-1}$')
ax[0,1].set_xlim([d-20,u])
ax[0,1].legend(loc='upper left',fontsize=8)

# Filter the OLR into ranges of vertical velocity at 500 hPa (W500).
OLR_1 = OLR.where((W500 > 0) & (W500 < 0.1),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((W500 >= 0.1) & (W500 < 1),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((W500 >= 1) & (W500 < 5),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((W500 >= 5),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((W500 < 0),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,2],label=r'$w_{500}$ $<$ 0 m s$^{-1}$')
sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},label=r'$w_{500}$ $\in$ [0,0.1] m s$^{-1}$',ax=ax[0,2])
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,2],label=r'$w_{500}$ $\in$ [0.1,1] m s$^{-1}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,2],label=r'$w_{500}$ $\in$ [1,5] m s$^{-1}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},ax=ax[0,2],label=r'$w_{500}$ $>$ 5 m s$^{-1}$')
ax[0,2].set_xlim([d-20,u])
ax[0,2].legend(loc='upper left',fontsize=8)

basedir = '/work/bb1018/b380873/tropic_run5_output/'
OLR_fi = xr.open_dataset(basedir + 'OLR_120-141_0.025deg.nc')
TQI_fi = xr.open_dataset(basedir + 'TQI_120-141_0.025deg.nc')
W_fi = xr.open_dataset(basedir + 'W_60-70_0.025deg_plev.nc')

# Make sure the times correspond (must be hourly given 3D w field)
OLR = OLR_fi.sel(time=W_fi.time).thb_t
TQI = TQI_fi.sel(time=W_fi.time).tqi
W250 = W_fi.sel(time=W_fi.time,plev_2=25000).w
W500 = W_fi.sel(time=W_fi.time,plev_2=50000).w
print(W250.min(),W250.max(),W500.min(),W500.max())

# Filter the OLR into ranges of ice water path (TQI). Flatten them to just look at statistics / density.
OLR_1 = OLR.where((TQI > 0) & (TQI < 10**(-4)),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((TQI >= 10**(-4)) & (TQI < 10**(-3)),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((TQI >= 10**(-3)) & (TQI < 10**(-2)),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((TQI >= 10**(-2)) & (TQI < 10**(-1)),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((TQI >= 10**(-1)),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,0])#label=r'Clear sky',ax=ax[1,0]
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,0])#,label=r'IWP $\in$ [10$^{-4}$,10$^{-3}$] g m$^{-2}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,0])#,label=r'IWP $\in$ [10$^{-3}$,10$^{-2}$] g m$^{-2}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,0])#,label=r'IWP $\in$ [10$^{-2}$,10$^{-1}$] g m$^{-2}$')
sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,0])#,label=r'IWP $>$ 10$^{-1}$ g m$^{-2}$')
ax[1,0].set_xlabel(r'OLR [W m$^{-2}$]',fontsize=fs)
ax[1,0].set_xlim([d-20,u])
ax[1,0].set_ylabel('Density',fontsize=fs)

# Filter the OLR into ranges of vertical velocity at 250 hPa (W250).
OLR_1 = OLR.where((W250 > 0) & (W250 < 0.1),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((W250 >= 0.1) & (W250 < 1),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((W250 >= 1) & (W250 < 5),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((W250 >= 5),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((W250 < 0),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,1])#,label=r'$w_{250}$ $<$ 0 m s$^{-1}$')
sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,1])#label=r'$w_{250}$ $\in$ [0,0.1] m s$^{-1}$')
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,1])#,label=r'$w_{250}$ $\in$ [0.1,1] m s$^{-1}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,1])#,label=r'$w_{250}$ $\in$ [1,5] m s$^{-1}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,1])#,label=r'$w_{250}$ $>$ 5 m s$^{-1}$')
ax[1,1].set_xlim([d-20,u])
ax[1,1].set_xlabel(r'OLR [W m$^{-2}$]',fontsize=fs)

# Filter the OLR into ranges of vertical velocity at 500 hPa (W500).
OLR_1 = OLR.where((W500 > 0) & (W500 < 0.1),drop=True).values.flatten()
OLR_1 = -1*OLR_1[~np.isnan(OLR_1)]
OLR_2 = OLR.where((W500 >= 0.1) & (W500 < 1),drop=True).values.flatten()
OLR_2 = -1*OLR_2[~np.isnan(OLR_2)]
OLR_3 = OLR.where((W500 >= 1) & (W500 < 5),drop=True).values.flatten()
OLR_3 = -1*OLR_3[~np.isnan(OLR_3)]
OLR_4 = OLR.where((W500 >= 5),drop=True).values.flatten()
OLR_4 = -1*OLR_4[~np.isnan(OLR_4)]
OLR_5 = OLR.where((W500 < 0),drop=True).values.flatten()
OLR_5 = -1*OLR_5[~np.isnan(OLR_5)]

sns.distplot(OLR_5,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[4],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,2])#,label=r'$w_{500}$ $<$ 0 m s$^{-1}$')
sns.distplot(OLR_1,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,2])#label=r'$w_{500}$ $\in$ [0,0.1] m s$^{-1}$')
sns.distplot(OLR_2,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,2])#,label=r'$w_{500}$ $\in$ [0.1,1] m s$^{-1}$')
sns.distplot(OLR_3,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,2])#,label=r'$w_{500}$ $\in$ [1,5] m s$^{-1}$')
sns.distplot(OLR_4,bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
    kde_kws={'shade':True,'linewidth':3},ax=ax[1,2])#,label=r'$w_{500}$ $>$ 5 m s$^{-1}$')
ax[1,2].set_xlim([d-20,u])
ax[1,2].set_xlabel(r'OLR [W m$^{-2}$]',fontsize=fs)

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
#fig.savefig('olrDist_IWP-w.pdf',bbox_inches='tight')
plt.show()
