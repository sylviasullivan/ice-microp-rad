import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys, time
from datetime import datetime

# The CERES and ERA5 block can be read in and arranged from OLRDist2.py
# Read in the ICON no2mom simulation.
#basedir = '/work/bb1131/b380873/tropic_run5_output/no2mom/'
#olr_file = basedir + 'OLR_120-141_0.025deg.nc'
#olr_data = xr.open_dataset(olr_file)
#olr_no2mom = np.abs(olr_data.lwflxall.values)
#zeit_no2mom = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
#exttime = datetime(2017,8,8,0,0)
#ii = np.argwhere(zeit_no2mom > np.datetime64(exttime))[0,0]
#olr_no2mom_sub = np.reshape(olr_no2mom[ii],(olr_no2mom.shape[1]*olr_no2mom.shape[2],))
#exttime = datetime(2017,8,8,6,0)
#ii = np.argwhere(zeit_no2mom > np.datetime64(exttime))[0,0]
#olr_no2mom_sub1 = np.reshape(olr_no2mom[ii],(olr_no2mom.shape[1]*olr_no2mom.shape[2],))

# Read in the CERES data.
basedir = '/work/bb1131/b380873/tropic_vis/obs/CERES/'
olr_file = basedir + 'CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20170801-20170831_full-domain.nc'
olr_data = xr.open_dataset(olr_file)
olr_ceres = np.abs(olr_data.toa_lw_all_1h.values)
lons_ceres = olr_data.lon
lats_ceres = olr_data.lat
zeit_ceres = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub = np.reshape(olr_ceres[ii],(olr_ceres.shape[1]*olr_ceres.shape[2],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub1 = np.reshape(olr_ceres[ii],(olr_ceres.shape[1]*olr_ceres.shape[2],))

# Read in the ERA5 values.
basedir = '/work/bb1131/b380873/tropic_vis/obs/ERA5/'
olr_file = basedir + 'ERA5_OLR_1deg[55-170]-20170805-20170809.nc'
olr_data = xr.open_dataset(olr_file)
olr_era5 = np.abs(olr_data.mtnlwrf.values)
zeit_era5 = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub1 = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))

# Read in the ICON simulation.
basedir = '/work/bb1131/b380873/tropic_run2_output/'
olr_file = basedir + 'OLR_TOA_all_1deg.nc'
olr_data = xr.open_dataset(olr_file)
olr_icon = np.abs(olr_data.lwflxall.values)
zeit_icon = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub1 = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))

print('ceres   ERA5    ICON')
print('Shape: ' + str(olr_ceres_sub.shape) + ' ' + str(olr_era5_sub.shape) + ' ' + str(olr_icon_sub.shape))
print(np.nanmean(olr_ceres_sub),np.nanmean(olr_era5_sub),np.nanmean(olr_icon_sub))
print(np.nanmedian(olr_ceres_sub),np.nanmedian(olr_era5_sub),np.nanmedian(olr_icon_sub))
print(np.nanstd(olr_ceres_sub),np.nanstd(olr_era5_sub),np.nanstd(olr_icon_sub))

def kl_divergence(p,q):
    return np.nansum(np.where(((p != 0) & (q != 0)), p*np.log2(p / q),0))

titre = ['ceres TOA OLR: 2017-0808T00:30:14', 'ERA5 TOA OLR: 2017-08-08T00:00:00', \
         'ICON TOA OLR: 2017-08-08T00:00:00', 'ceres TOA OLR: 2017-0808T06:30:14', \
         'ERA5 TOA OLR: 2017-08-08T06:00:00', 'ICON TOA OLR: 2017-08-08T06:00:00']
olr = [olr_ceres_sub, olr_era5_sub, olr_icon_sub, olr_ceres_sub1, olr_era5_sub1, olr_icon_sub1]
farbe = ['gold','blue','red','gold','blue','red']
let = ['(a)','(c)','(e)','(b)','(d)','(f)']
# KL divergence values
#KL = [
fs = 13

fig, ax = plt.subplots(3,2,figsize=(7,9))
#fig.subplots_adjust(hspace=0.1, wspace=0, top=0.95, left=0.1) #top=0.925
h = np.zeros((49,6))

c = 0
for j in np.arange(2):
    for i in np.arange(3):
        ax[i,j].set_title(titre[c],fontsize=fs-2)
        wgts = np.ones_like(olr[c])/float(len(olr[c]))*100
        ax[i,j].hist(olr[c],bins=np.linspace(80,360,50),color=farbe[c],edgecolor='k',\
                weights=wgts)
        h[:,c],bar_edges = np.histogram(olr[c],bins=np.linspace(80,360,50),weights=wgts)
        if i == 2:
           ax[i,j].set_xlabel('OLR [W m$^{-2}$]',fontsize=fs)
        if j == 0:
           ax[i,j].set_ylabel('Probability',fontsize=fs) 
        if j == 0 and i == 1:
           ax[i,j].text(0.05,0.8,'KL(CERES | ERA5) = 3.99',transform=ax[i,j].transAxes)
        if j == 1 and i == 1:
           ax[i,j].text(0.05,0.8,'KL(CERES | ERA5) = 5.30',transform=ax[i,j].transAxes)
        if j == 0 and i == 2:
           ax[i,j].text(0.05,0.78,'KL(CERES | ICON) = '
                                 '\n'
                                 '10.04',transform=ax[i,j].transAxes)
           ax[i,j].text(0.05,0.63,'KL(ERA5 | ICON) = '
                                  '\n'
                                  '9.15',transform=ax[i,j].transAxes)
        if j == 1 and i == 2:
           ax[i,j].text(0.05,0.78,'KL(CERES | ICON) = '
                                 '\n'
                                 '14.37',transform=ax[i,j].transAxes)
           ax[i,j].text(0.05,0.63,'KL(ERA5 | ICON) = '
                                 '\n'
                                 '9.83',transform=ax[i,j].transAxes)
        ax[i,j].set_ylim([0,15])
        ax[i,j].text(0.05,0.92,let[c],fontsize=fs,fontweight='bold',transform=ax[i,j].transAxes)
        c = c + 1

print(h[:,2])
#print(np.stack(bar_edges,h[:,2]))
print('ERA5 from CERES (0h): ' + str(kl_divergence(h[:,0],h[:,1])))
print('CERES from ERA5 (0h): ' + str(kl_divergence(h[:,1],h[:,0])))
print('ICON from CERES (0h): ' + str(kl_divergence(h[:,0],h[:,2])))
print('CERES from ICON (0h): ' + str(kl_divergence(h[:,2],h[:,0])))
print('ICON from ERA5 (0h): ' + str(kl_divergence(h[:,1],h[:,2])))
print('ERA5 from ICON (0h): ' + str(kl_divergence(h[:,2],h[:,1])))
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('ERA5 from CERES (6h): ' + str(kl_divergence(h[:,3],h[:,4])))
print('CERES from ERA5 (6h): ' + str(kl_divergence(h[:,4],h[:,3])))
print('ICON from CERES (6h): ' + str(kl_divergence(h[:,3],h[:,5])))
print('CERES from ICON (6h): ' + str(kl_divergence(h[:,5],h[:,3])))
print('ICON from ERA5 (6h): ' + str(kl_divergence(h[:,4],h[:,5])))
print('ERA5 from ICON (6h): ' + str(kl_divergence(h[:,5],h[:,4])))



fig.tight_layout(w_pad=0.2)
#fig.savefig('olr-distribution[1deg].pdf',bbox_inches='tight')
plt.show()
