import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import sys

def prefix(i):
    if i < 10:
       prefix = '00'
    else:
       prefix = '0'
    return prefix

basedir = '/work/bb1018/b380873/tropic_vis/obs/POSIDON/'
posidon_winds = pd.read_csv(basedir + 'posidon_vertical_wind2.dat',sep='\s+')

# Filter out the w values only for z > 15 km
posidon_winds = posidon_winds[(posidon_winds['z'] >= 15)]
w_obs = posidon_winds['w']

# Calculate the normalized / relative frequency of these w values
wgts = np.ones_like(w_obs)/len(w_obs)
h_obs, bin_edges = np.histogram(w_obs, bins=np.linspace(-3,3,100), weights=wgts) # density=True)

basedir = '/work/bb1018/b380873/tropic_vis/obs/ATTREX/'
attrex_winds = pd.read_csv(basedir + 'attrex3_vertical_wind2.dat',sep='\s+')

# Filter out the w values only for z > 15 km
attrex_winds = attrex_winds[(attrex_winds['z'] >= 15)]
w_obs = attrex_winds['w']

# Calculate the normalized / relative frequency of these w values
wgts = np.ones_like(w_obs)/len(w_obs)
h_obs2, bin_edges = np.histogram(w_obs, bins=np.linspace(-3,3,100), weights=wgts) # density=True)

# Read in the vertical velocities from all 30 trajectories
#basedir = '/work/bb1018/b380873/traj_output/full51h_fast/51h_trim/'
#fi = basedir + 'traj_tst00000450_p001_trim.nc'
make_h_sim = False
if make_h_sim == True:
   h_sim = np.zeros((30,99))
   for i in np.arange(1,31):
       print(i)
       basedir = '/scratch/b/b380873/traj_full51h_fast/'
       fi = basedir + 'traj_tst00000450_p' + prefix(i) + str(i) + '_trim.nc'
       z = xr.open_dataset(fi)['alt']
       #w_sim = xr.open_dataset(fi)['w_v'].where((z > 10000) & (z <= 15000)).values
       w_sim = xr.open_dataset(fi)['w_v'].where((z > 15000)).values
       w_sim = w_sim[~np.isnan(w_sim)]
       print(np.nanmean(w_sim),np.nanstd(w_sim))

       # Calculate the normalized / relative frequency of these w values
       wgts = np.ones_like(w_sim)/len(w_sim)
       h_sim[i-1], _ = np.histogram(w_sim, bins=np.linspace(-3,3,100), weights=wgts)
   np.save('../output/traj_w_histogram_sim_extract.npy',h_sim)
else:
   h_sim = np.load('../output/traj_w_histogram_sim.npy')  #'traj_w_histogram_sim_extract.npy'
   h_sim2 = np.load('../output/traj_w_histogram_sim-10-15.npy')

fs = 16
fig = plt.figure(figsize=(7,6))
plt.plot((bin_edges[1:] + bin_edges[:-1])/2., h_obs, color='k')#, lineweight=1.2)
plt.scatter((bin_edges[1:] + bin_edges[:-1])/2., h_obs, marker='x', color='k')
plt.plot((bin_edges[1:] + bin_edges[:-1])/2., h_obs2, color='gray',label='ATTREX')
plt.scatter((bin_edges[1:] + bin_edges[:-1])/2., h_obs2, marker='x', color='gray')
plt.plot((bin_edges[1:] + bin_edges[:-1])/2., h_sim.T, color='r')
plt.plot((bin_edges[1:] + bin_edges[:-1])/2., h_sim2.T,color='b')
#plt.scatter((bin_edges[1:] + bin_edges[:-1])/2., h_sim, marker='x', color='r')
plt.gca().set_yscale('log')
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.gca().tick_params('both',labelsize=fs)
plt.ylim([0.0001, 0.3])
plt.ylabel('Relative frequency',fontsize=fs)
plt.xlim([-3, 3])
plt.xlabel(r'Vertical velocity [m s$^{-1}$]',fontsize=fs)
#fig.savefig('../output/traj_POSIDON_w_histogram_small-range.pdf',bbox_inches='tight')
#plt.show()

fig2 = plt.figure(figsize=(7,6))
plt.plot((bin_edges[1:] + bin_edges[:-1])/2., h_sim.T)#, color='r')
#plt.scatter((bin_edges[1:] + bin_edges[:-1])/2., h_sim, marker='x', color='r')
plt.gca().set_yscale('log')
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.gca().tick_params('both',labelsize=fs)
plt.ylim([0.0001, 0.4])
plt.ylabel('Relative frequency',fontsize=fs)
plt.xlim([-1.5, 1.5])
plt.xlabel(r'Vertical velocity [m s$^{-1}$]',fontsize=fs)
#fig2.savefig('../output/traj_w_histogram_small-range.pdf',bbox_inches='tight')
plt.show()


