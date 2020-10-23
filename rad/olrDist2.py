import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys, time
from datetime import datetime
import seaborn as sns

# Read in the CERES data.
basedir = '/work/bb1018/b380873/tropic_vis/obs/CERES/'
olr_file = basedir + 'CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20170801-20170831_full-domain.nc'
olr_data = xr.open_dataset(olr_file)
print(olr_data.dims)
olr_ceres = np.abs(olr_data.toa_lw_all_1h.values)[:,:61]
zeit_ceres = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub = np.reshape(olr_ceres[ii],(olr_ceres.shape[1]*olr_ceres.shape[2],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_ceres > np.datetime64(exttime))[0,0]
olr_ceres_sub1 = np.reshape(olr_ceres[ii],(olr_ceres.shape[1]*olr_ceres.shape[2],))

# Read in the ERA5 values.
basedir = '/work/bb1018/b380873/tropic_vis/obs/ERA5/'
olr_file = basedir + 'ERA5_OLR_1deg[55-170]-20170805-20170809.nc'
olr_data = xr.open_dataset(olr_file)
#olr_era5 = np.abs(olr_data.ttr.values)[:,:61]/3600
olr_era5 = np.abs(olr_data.mtnlwrf.values)[:,:61]
zeit_era5 = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
olr_era5_sub1 = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))

# Read in the ICON-1mom simulation.
basedir = '/work/bb1018/b380873/tropic_run2_output/'
olr_file = basedir + 'OLR_TOA_all.nc' # 1deg
olr_data = xr.open_dataset(olr_file)
olr_icon = np.abs(olr_data.lwflxall.values)[:,:61]
zeit_icon = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub1 = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))

# Read in the ICON-2mom simulation.
basedir = '/work/bb1018/b380873/tropic_run5_output/'
olr_file2 = basedir + 'OLR_120-141_0.025deg.nc' # 1deg
olr_data2 = xr.open_dataset(olr_file2)
olr_icon2 = np.abs(olr_data2.thb_t.values)
zeit_icon2 = olr_data2.time.values

# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon2 >= np.datetime64(exttime))[0,0]
olr_icon_sub2 = np.reshape(olr_icon2[ii],(olr_icon2.shape[1]*olr_icon2.shape[2],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon2 >= np.datetime64(exttime))[0,0]
olr_icon_sub12 = np.reshape(olr_icon2[ii],(olr_icon2.shape[1]*olr_icon2.shape[2],))

print('CERES   ERA5    ICON-1mom   ICON-2mom')
print('Shape: ' + str(olr_ceres_sub.shape) + ' ' + str(olr_era5_sub.shape) + ' ' + str(olr_icon_sub.shape) +
      ' ' + str(olr_icon_sub2.shape))
print(np.nanmean(olr_ceres_sub),np.nanmean(olr_era5_sub),np.nanmean(olr_icon_sub),np.nanmean(olr_icon_sub2))
print(np.nanmedian(olr_ceres_sub),np.nanmedian(olr_era5_sub),np.nanmedian(olr_icon_sub),np.nanmedian(olr_icon_sub2))
print(np.nanstd(olr_ceres_sub),np.nanstd(olr_era5_sub),np.nanstd(olr_icon_sub),np.nanstd(olr_icon_sub2))

def kl_divergence(p,q):
    return np.nansum(np.where(((p != 0) & (q != 0)), p*np.log2(p / q),0))

titre = ['2017-08-08 00:00', '2017-08-08 06:00']
olr = [olr_ceres_sub, olr_ceres_sub1, olr_era5_sub, olr_era5_sub1, olr_icon_sub, olr_icon_sub1, olr_icon_sub2, olr_icon_sub12]
farbe = ['gold','blue','red','green']
let = ['a','b']
fs = 13

fig, ax = plt.subplots(1,2,figsize=(9,5))
#fig.subplots_adjust(hspace=0.1, wspace=0, top=0.95, left=0.1) #top=0.925
h = np.zeros((49,4))

for j in np.arange(2):
    ax[j].set_title(titre[j],fontsize=fs-2)
    d = 80; u = 360; b = 50
    sns.distplot(olr[j],bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[0],
        kde_kws={'shade':True,'linewidth':3},ax=ax[j],label='CERES')
    sns.distplot(olr[j+2],bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[1],
        kde_kws={'shade':True,'linewidth':3},ax=ax[j],label='ERA5')
    sns.distplot(olr[j+4],bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[2],
        kde_kws={'shade':True,'linewidth':3},ax=ax[j],label='ICON-1mom')
    sns.distplot(olr[j+6],bins=np.linspace(d,u,b),kde=True,hist=False,color=farbe[3],
        kde_kws={'shade':True,'linewidth':3},ax=ax[j],label='ICON-2mom')
    #ax[j].hist(olr[c],bins=np.linspace(80,360,50),color=farbe[c],edgecolor='k',\
    #        weights=wgts)
    if j == 0:
       wgts = np.ones_like(olr[j])/float(len(olr[j]))*100
       h[:,j],bar_edges = np.histogram(olr[j],bins=np.linspace(80,360,50),weights=wgts)
       wgts = np.ones_like(olr[j+2])/float(len(olr[j+2]))*100
       h[:,j+1],bar_edges = np.histogram(olr[j+2],bins=np.linspace(80,360,50),weights=wgts)
       wgts = np.ones_like(olr[j+4])/float(len(olr[j+4]))*100
       h[:,j+2],bar_edges = np.histogram(olr[j+4],bins=np.linspace(80,360,50),weights=wgts)
       wgts = np.ones_like(olr[j+6])/float(len(olr[j+6]))*100
       h[:,j+3],bar_edges = np.histogram(olr[j+6],bins=np.linspace(80,360,50),weights=wgts)
    ax[j].set_xlabel('OLR [W m$^{-2}$]',fontsize=fs)
    ax[j].set_xlim([80,360])
    if j == 0:
       ax[j].set_ylabel('Probability density',fontsize=fs)
    #if j == 0 and i == 1:
    #   ax[i,j].text(0.05,0.8,'KL(CERES | ERA5) = 3.99',transform=ax[i,j].transAxes)
    #if j == 1 and i == 1:
    #   ax[i,j].text(0.05,0.8,'KL(CERES | ERA5) = 5.30',transform=ax[i,j].transAxes)
    #if j == 0 and i == 2:
    #   ax[i,j].text(0.05,0.78,'KL(CERES | ICON) = '
    #                 '\n'
    #                             '10.04',transform=ax[i,j].transAxes)
    #       ax[i,j].text(0.05,0.63,'KL(ERA5 | ICON) = '
    #                              '\n'
    #                              '9.15',transform=ax[i,j].transAxes)
    #    if j == 1 and i == 2:
    #       ax[i,j].text(0.05,0.78,'KL(CERES | ICON) = '
    #                             '\n'
    #                             '14.37',transform=ax[i,j].transAxes)
    #       ax[i,j].text(0.05,0.63,'KL(ERA5 | ICON) = '
    #                             '\n'
    #                             '9.83',transform=ax[i,j].transAxes)
    ax[j].text(0.05,0.92,let[j],fontsize=fs,fontweight='bold',transform=ax[j].transAxes)
    if j == 1:
       ax[j].legend(loc='center left')
    else:
       ax[j].get_legend().remove()

#print(h[:,2])
#print(np.stack(bar_edges,h[:,2]))
print('ERA5 from CERES (0h): ' + str(kl_divergence(h[:,0],h[:,1])))
#print('CERES from ERA5 (0h): ' + str(kl_divergence(h[:,1],h[:,0])))
print('ICON 1mom from CERES (0h): ' + str(kl_divergence(h[:,0],h[:,2])))
print('ICON 2mom from CERES (0h): ' + str(kl_divergence(h[:,0],h[:,3])))
#print('CERES from ICON (0h): ' + str(kl_divergence(h[:,2],h[:,0])))
#print('ICON from ERA5 (0h): ' + str(kl_divergence(h[:,1],h[:,2])))
#print('ERA5 from ICON (0h): ' + str(kl_divergence(h[:,2],h[:,1])))
#print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
#print('ERA5 from CERES (6h): ' + str(kl_divergence(h[:,3],h[:,4])))
#print('CERES from ERA5 (6h): ' + str(kl_divergence(h[:,4],h[:,3])))
#print('ICON from CERES (6h): ' + str(kl_divergence(h[:,3],h[:,5])))
#print('CERES from ICON (6h): ' + str(kl_divergence(h[:,5],h[:,3])))
#print('ICON from ERA5 (6h): ' + str(kl_divergence(h[:,4],h[:,5])))
#print('ERA5 from ICON (6h): ' + str(kl_divergence(h[:,5],h[:,4])))



fig.tight_layout(w_pad=0.2)
#fig.savefig('../output/olr-distribution_115e.pdf',bbox_inches='tight')
plt.show()
