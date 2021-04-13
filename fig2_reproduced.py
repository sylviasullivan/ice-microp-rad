import matplotlib.pyplot as plt
from datetime import datetime
from functools import reduce
import numpy.ma as ma
import scipy.interpolate
import numpy as np
import xarray as xr
import sys,time,os

# Function to calculate the saturation vapor pressure over ice.
def satVapP_ice(t_in):
    a1 = 9.550426
    a2 = -5723.265
    a3 = 3.53068
    a4 = -0.00728332
    psatI = a1 + a2/t_in + a3*np.log(t_in) + a4*t_in
    psatI = np.exp(psatI)
    return psatI

daten = xr.open_dataset('obs/stratoclim2017.geophysika.0808_1.master.ci_eval.nc')
zeit = daten.time

# In Figure 2 of Lee et al., only values between 6:20 and 6:48 UTC are used.
# Vertical profile of water vapour [ppmv] for first panel.
# Filter for lack of altitude data, recorded as -99999.
alt = daten['BEST:ALT'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
h2o_flash = daten['BEST:H2O_gas'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
h2o_fish = daten['BEST:H2O_enh'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
iwc = daten['BEST:IWC'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
temp = daten['BEST:TEMP'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
theta = daten['BEST:THETA'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
rhice_flash = daten['BEST:RH_ice_gas'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))
rhice_fish = daten['BEST:RH_ice_enh'].sel(time=slice(datetime(2017,8,8,6,20),datetime(2017,8,8,6,48)))

sys.exit()
# Different altitude variables extracted with different filters.
alt1 = alt.where((alt > 0) & (h2o_flash > 0) & (h2o_fish > 0))
h2o_flash = h2o_flash.where((alt > 0) & (h2o_flash > 0) & (h2o_fish > 0))
h2o_fish = h2o_fish.where((alt > 0) & (h2o_flash > 0) & (h2o_fish > 0))

alt2 = alt.where((alt > 0) & (iwc > 0))
iwc = iwc.where((alt > 0) & (iwc > 0))

alt3 = alt.where((temp > 0) & (theta > 0))
temp = temp.where((temp > 0) & (theta > 0))
theta = theta.where((temp > 0) & (theta > 0))

alt4 = alt.where((rhice_flash > 0) & (rhice_fish > 0))
rhice_flash = rhice_flash.where((rhice_flash > 0) & (rhice_fish > 0))
rhice_fish = rhice_fish.where((rhice_flash > 0) & (rhice_fish > 0))

# If binning between <u> and <d> with <n> bins, which elements go in which bin?
# Make a multidimensional list of alt and h2o values in each.
u = 14000
d = 22000
n = 100
i1 = np.digitize(alt1,bins=np.linspace(u,d,n))
i2 = np.digitize(alt2,bins=np.linspace(u,d,n))

# All elements
alt1_list = [[] for i in np.arange(n)]
h2o_flash_list = [[] for i in np.arange(n)]
h2o_fish_list = [[] for i in np.arange(n)]
alt2_list = [[] for i in np.arange(n)]
iwc_list = [[] for i in np.arange(n)]

# Their mean
h2o_flash_m = np.zeros((n,)); h2o_flash_m[:] = np.nan
h2o_fish_m = np.zeros((n,)); h2o_fish_m[:] = np.nan
iwc_m = np.zeros((n,)); iwc_m[:] = np.nan

for elem_idx, group_idx in enumerate(i1):
    alt1_list[group_idx-1].append(alt1[elem_idx])
    h2o_flash_list[group_idx-1].append(h2o_flash[elem_idx])
    h2o_fish_list[group_idx-1].append(h2o_fish[elem_idx])
    alt2_list[group_idx-1].append(alt2[elem_idx])
    iwc_list[group_idx-1].append(iwc[elem_idx])

# Take the mean across the lists if there are at least 5 elements.
for i in np.arange(n):
    if (len(h2o_flash_list[i]) > 5):
       h2o_flash_m[i] = sum(h2o_flash_list[i]) / len(h2o_flash_list[i])
    if (len(h2o_fish_list[i]) > 5):
       h2o_fish_m[i] = sum(h2o_fish_list[i]) / len(h2o_fish_list[i])
    if (len(iwc_list[i]) > 5):
       iwc_m[i] = sum(iwc_list[i]) / len(iwc_list[i])

# Read in the simulation values.
# Filter for values between 25-26.5 N lat and 85-85.5 E lon and store only positive values.
basedir = '/work/bb1018/b380873/tropic_run5_output/'
sims = xr.open_dataset(basedir + 'QV_QI_60-70_0.025deg_HL.nc')
qi_s = sims.qi.sel(time=datetime(2017,8,8,6,0),lat=slice(25,26.5),lon=slice(85,85.5))
qv_s = sims.qv.sel(time=datetime(2017,8,8,6,0),lat=slice(25,26.5),lon=slice(85,85.5))
qi_s = qi_s.where(qi_s > 0)
qv_s = qv_s.where(qv_s > 0) # (100, 60, 20)
z_s = sims.height
sims = xr.open_dataset(basedir + 'UVT_60-70_0.025deg_HL.nc')
T_s = sims.temp.sel(time=datetime(2017,8,8,6,0),lat=slice(25,26.5),lon=slice(85,85.5))
P_s = sims.air_pressure.sel(time=datetime(2017,8,8,6,0),lat=slice(25,26.5),lon=slice(85,85.5))

# Multiply qv_s by 10^6 to translate kg kg-1 to ppmv and take its lat-lon mean.
mw_dryair = 28.97*1000    # kg air (mol air)-1
mw_watvap = 18.02*1000    # kg wv (mol wv)-1
conv = mw_dryair / mw_watvap
qv_s_ppmv = qv_s.mean({'lat','lon'}) * conv * 10**6
qi_s_ppmv = qi_s.mean({'lat','lon'}) * conv * 10**6
T_s_profile = T_s.mean({'lat','lon'})
P_s_profile = P_s.mean({'lat','lon'})
# Poisson constant assumed to be 2/7, the ratio of the gas constant to specific cp
# at constant pressure for an ideal diatomic gas. P_0 is a reference pressure.
kappa = 0.285714
P_0 = 100000
theta_s_profile = (T_s * (P_0 / P_s) ** kappa).mean({'lat','lon'})
# Use Murphy & Koop equation for p_sat,ice. eps = ratio of water vapor and dry air
# molar masses.
e = qv_s / (1 - qv_s)
eps = 0.622
e_si = eps * satVapP_ice(T_s) / (P_s - satVapP_ice(T_s))
rhi_s_profile = (e / e_si).mean({'lat','lon'}) * 100

fs = 12
fig, ax = plt.subplots(nrows=1,ncols=4,figsize=(13,5.5))
ax[0].grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
ax[0].scatter(h2o_flash,alt1/1000,color='k',s=10,marker='*',label='FLASH')
ax[0].scatter(h2o_fish,alt1/1000,color='gray',s=10,marker='d',label='FISH')
ax[0].scatter(qv_s_ppmv,z_s/1000,marker='o',color='r',s=15,label='ICON')
ax[0].plot(h2o_flash_m,np.linspace(u,d,n)/1000,linewidth=1.25,color='k')
ax[0].plot(h2o_fish_m,np.linspace(u,d,n)/1000,linewidth=1.25,color='gray')
ax[0].plot([4,4],[14,21],linewidth=0.75,linestyle='--',color='k')
ax[0].set_xlim([2,10])
ax[0].set_ylim([15,19])
ax[0].set_xlabel('Water vapour [ppmv]',fontsize=fs)
ax[0].set_ylabel('Altitude [km]',fontsize=fs)
ax[0].legend(loc='upper right')
plt.tick_params(labelsize=fs)

ax[1].grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
ax[1].scatter(iwc,alt2/1000,color='k',s=10,marker='*')
ax[1].scatter(qi_s_ppmv,z_s/1000,marker='o',color='r',s=15)
ax[1].plot(iwc_m,np.linspace(u,d,n)/1000,linewidth=1.25,color='k')
ax[1].set_ylim([15,19])
ax[1].set_xlim([0,2])
ax[1].set_xlabel('Ice [eq.ppmv]',fontsize=fs)
plt.gca().set_yticklabels([])

ax[2].grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
ax2 = ax[2].twiny()
ax[2].scatter(temp-273,alt3/1000,color='k',s=10,marker='*',label=r'in-situ $T$')
ax2.scatter(theta,alt3/1000,color='b',s=10,marker='*',label=r'in-situ $\theta$')
ax[2].scatter(T_s_profile-273,z_s/1000,marker='o',color='r',s=15,label='ICON $T$')
ax2.scatter(theta_s_profile,z_s/1000,marker='o',color='deepskyblue',s=15,label=r'ICON $\theta$')
ax[2].set_ylim([15,19])
ax[2].set_xlim([-85,-65]); ax2.set_xlim([350,430])
ax[2].set_xlabel('Temperature [deg C]',fontsize=fs)
ax2.set_xlabel('Potential temperature [K]',fontsize=fs,color='b')
ax2.spines['top'].set_color('b'); ax2.tick_params(axis='x',colors='b')
ax[2].legend(loc='center right')
ax2.legend(loc='lower right')
plt.gca().set_yticklabels([])

ax[3].grid(b=True,which='both',axis='both',color='gray',linewidth=0.5)
ax[3].scatter(rhice_fish,alt4/1000,color='k',s=10,marker='*',label='FISH')
ax[3].scatter(rhice_flash,alt4/1000,color='gray',s=10,marker='*',label='FLASH')
ax[3].scatter(rhi_s_profile,z_s/1000,marker='o',color='r',s=15)
ax[3].set_ylim([15,19])
ax[3].set_xlim([40,150])
ax[3].set_xlabel('RHice [%]',fontsize=fs)
ax[3].legend(loc='upper right')
plt.gca().set_yticklabels([])

plt.tight_layout()
#plt.savefig('StratoClim_ICON_profiles.pdf',dpi=200,bbox_inches='tight')
plt.show()
