import matplotlib.dates as mdates
from datetime import datetime
from datetime import timedelta
import time
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

basedir = '/work/bb1131/b380873/traj_output/test2h/'
clams = xr.open_dataset(basedir + 'cirrus_tst00000751_p001_trim_clams.nc')
icon = xr.open_dataset(basedir + 'traj_tst00000751_p001_trim.nc')
clams.rename({'T':'temp'},inplace=True)

clams_t = clams.time.values
clams_t = np.array([t.astype('int') for t in clams_t])/1000000000 - 946684800
icon_t = icon.rtime.values

clams_qi = clams.IWC_hom + clams.IWC_het + clams.IWC_pre
# kg kg-1, ice mass mixing ratio from homog. + heterog. nucleation and preexisting ice
# Dims = time 120, NPARTS 5433

icon_qi = icon.qi
# kg kg-1, ice mass mixing ratio
# Dims = time_step 300, traj_id 5433

clams_temp = clams.temp.values
icon_temp = icon.t.values

# How many of the trajectories do you want to visualize in the background?
n = 500
colors1 = plt.cm.Reds(np.linspace(0,1,n))
colors2 = plt.cm.Blues(np.linspace(0,1,n))

icon_qi_cubic = np.zeros((clams_qi.shape[1],clams_qi.shape[0]))
for i in np.arange(clams_qi.shape[1]):
    cubic_interp = interp1d(icon_t,icon_qi[:,i].values,kind='cubic')
    # Factor of 10**6 converts to mg per kg
    icon_qi_cubic[i] = cubic_interp(clams_t)*10**6

#print(clams_qi.max(),clams_qi.min())
#print(icon_qi.max(),icon_qi.min())

fs = 15
fig = plt.figure(figsize=(11,5.5))
dates = [datetime(2017,8,5,16,0) + timedelta(seconds=t) for t in clams_t]

plt.gca().xaxis_date()
for i in np.arange(n):
    plt.plot(dates,icon_qi_cubic[i],color=colors1[i],alpha=0.25,lw=0.25)
    plt.plot(dates,clams_qi[:,i]*10**6,color=colors2[i],alpha=0.25,lw=0.25)
    # Factor of 10**6 converts to mg per kg
plt.plot(dates,np.nanmean(icon_qi_cubic,axis=0),color='red',lw=1.25,label='ICON')
plt.plot(dates,np.nanmean(clams_qi*10**6,axis=1),color='blue',lw=1.25,label='CLaMS-ice')
plt.gca().set_ylim([10**(-2),20])
plt.gca().set_ylabel(r'Ice mass mixing ratio [mg kg$^{-1}$]',fontsize=fs)
plt.gca().set_xlabel('Trajectory time',fontsize=fs)
years_fmt = mdates.DateFormatter('%Y-%m-%d %H:%M')
plt.gca().xaxis.set_major_formatter(years_fmt)
fig.autofmt_xdate()
plt.legend(loc='upper left')
plt.savefig('../output/CLAMS-ICON-traj_time.png')
plt.show()
sys.exit()

# Classify the trajectories by their mean temperature.
meant = np.zeros((icon_qi_cubic.shape[0],))
for i in np.arange(icon_qi_cubic.shape[0]):
    t1 = np.nanmean(icon_temp[:,i])
    t2 = np.nanmean(clams_temp[:,i])
    if((238 <= t1 < 225) & (238 <= t2 < 225)):
       meant[i] = 1
    elif((225 <= t1 < 210) & (225 <= t2 < 210)):
       meant[i] = 2
    elif((t1 <= 210) & (t2 <= 210)):
       meant[i] = 3

ii = np.argwhere(meant == 3)
ax[1].xaxis_date()
ax[1].plot(dates,np.nanmean(icon_qi_cubic[ii[:,0]],axis=0),color='red',lw=1.25)
ax[1].plot(dates,np.nanmean(clams_qi[:,ii[:,0]],axis=1),color='blue',lw=1.25)
print(np.nanmean(icon_qi_cubic[ii[:,0]],axis=0))
print(np.nanmean(clams_qi[:,ii[:,0]],axis=1))
#plt.gca().set_ylim([10**(-2),10**2])
#plt.gca().set_yscale('log')
#fig.autofmt_xdate()

plt.show()

