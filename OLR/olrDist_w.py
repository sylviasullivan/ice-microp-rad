import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import sys, time
from datetime import datetime

# Read in the ERA5 values.
#basedir = '/work/bb1131/b380873/tropic_vis/obs/'
#olr_file = basedir + 'ERA5_OLR_1deg[55-170]-20170805-20170809.nc'
#olr_data = xr.open_dataset(olr_file)
#olr_era5 = np.abs(olr_data.mtnlwrf.values)
#zeit_era5 = olr_data.time.values

# Extract 8 August 2017 at midnight / 6 am
#exttime = datetime(2017,8,8,0,0)
#ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
#olr_era5_sub = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))
#exttime = datetime(2017,8,8,6,0)
#ii = np.argwhere(zeit_era5 >= np.datetime64(exttime))[0,0]
#olr_era5_sub1 = np.reshape(olr_era5[ii],(olr_era5.shape[1]*olr_era5.shape[2],))

# Read in the ICON simulation.
basedir = '/work/bb1131/b380873/tropic_run2_output/'
olr_file = basedir + 'olr_ALL-0102-0122.nc'
olr_data = xr.open_dataset(olr_file)
olr_icon = np.abs(olr_data.lwflxall.values)[0]   # dims = (11,1800,4600)

w_file = basedir + 'w500_ALL-0051-0061.nc'
w_data = xr.open_dataset(w_file)
w_icon = w_data.w.values[0]                      # dims = (11,1800,4600)

# Generate a field of convective and quiescent OLR values.
olr_icon_conv = np.zeros((3,11,1800,4600))
wlim = [0.1, 1, 10]
for i in np.arange(3):
    olr_icon_conv[i] = np.where(w_icon > wlim[i], olr_icon, np.nan)
olr_icon_situ = np.where(w_icon <= 0, olr_icon, np.nan)

fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(7,7))
u = 360; d = 80; b = 50
fs = 13
y1 = 0; y2 = 15
let = ['(a)','(b)','(c)']

# Filter NaN instances from the convective OLR values. Normalize them by count.
a = np.asarray([[0,0],[0,1],[1,0]])
for i in np.arange(3):
    olr_icon_conv_nonnan = olr_icon_conv[i]
    olr_icon_conv_nonnan = olr_icon_conv_nonnan[~np.isnan(olr_icon_conv_nonnan)]
    wgts = np.ones_like(olr_icon_conv_nonnan)/float(len(olr_icon_conv_nonnan))*100
    ax[a[i,0],a[i,1]].hist(olr_icon_conv_nonnan,bins=np.linspace(d,u,b),facecolor='red',alpha=0.6,weights=wgts)
    if i != 1:
       ax[a[i,0],a[i,1]].set_ylabel('Probability',fontsize=fs)
    ax[a[i,0],a[i,1]].set_xlabel(r'OLR (w > ' + str(wlim[i]) + ' m s$^{-1}$)',fontsize=fs)
    ax[a[i,0],a[i,1]].set_ylim([y1,y2])
    ax[a[i,0],a[i,1]].text(0.05,0.92,let[i],weight='bold',fontsize=fs+2,transform=ax[a[i,0],a[i,1]].transAxes)

# Filter NaN instances from the quiescent OLR values. Normalize them by count.
olr_icon_situ_nonnan = olr_icon_situ[~np.isnan(olr_icon_situ)]
wgts = np.ones_like(olr_icon_situ_nonnan)/float(len(olr_icon_situ_nonnan))*100
ax[1,1].hist(olr_icon_situ_nonnan,bins=np.linspace(d,u,b),facecolor='blue',alpha=0.6,weights=wgts)
ax[1,1].set_xlabel(r'OLR (w < 0 m s$^{-1})$',fontsize=fs)
ax[1,1].set_ylim([y1,y2])
ax[1,1].text(0.05,0.92,'(d)',weight='bold',fontsize=fs+2,transform=ax[1,1].transAxes)

fig.savefig('olrDist_w.pdf',bbox_inches='tight')
plt.show()

sys.exit()



# Extract 8 August 2017 at midnight / 6 am
exttime = datetime(2017,8,8,0,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))
exttime = datetime(2017,8,8,6,0)
ii = np.argwhere(zeit_icon >= np.datetime64(exttime))[0,0]
olr_icon_sub1 = np.reshape(olr_icon[ii,0],(olr_icon.shape[2]*olr_icon.shape[3],))

print('CERES   ERA5    ICON')
print('Shape: ' + str(olr_ceres_sub.shape) + ' ' + str(olr_era5_sub.shape) + ' ' + str(olr_icon_sub.shape))
print(np.nanmean(olr_ceres_sub),np.nanmean(olr_era5_sub),np.nanmean(olr_icon_sub))
print(np.nanmedian(olr_ceres_sub),np.nanmedian(olr_era5_sub),np.nanmedian(olr_icon_sub))
print(np.nanstd(olr_ceres_sub),np.nanstd(olr_era5_sub),np.nanstd(olr_icon_sub))

def kl_divergence(p,q):
    return np.nansum(np.where(((p != 0) & (q != 0)), p*np.log2(p / q),0))

titre = ['CERES TOA OLR: 2017-0808T00:30:14', 'ERA5 TOA OLR: 2017-08-08T00:00:00', \
         'ICON TOA OLR: 2017-08-08T00:00:00', 'CERES TOA OLR: 2017-0808T06:30:14', \
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
print(np.stack(bar_edges,h[:,2]))
sys.exit()
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
