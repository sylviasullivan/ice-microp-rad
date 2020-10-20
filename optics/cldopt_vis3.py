# Script to plot Baum-Yang optical properties (ice only) and their deviation
# from standard RRTM ones (Fu scheme).

import sys
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import xarray as  xr
from matplotlib import cm,colors

# Print the min, mean, max of these properties?
print_opt = True
# Plot absorption coeff, single-scattering albedo, asymmetry parameter as a function of wavelength?
print_plot1 = True
# Plot absorption coeff, single-scattering albedo, asymmetry parameter as a function of effective radius?
print_plot2 = False

basedir = '/work/bb1131/b380873/rrtm/'
#fi = basedir + 'baum_yang_cldopt.nc'
fi = basedir + 'rrtm_cldopt.nc'
cldopts = xr.open_dataset(fi)

# Dimensions -- effective radius and wavelength
r_i = cldopts['re_crystal'].values   #  4 - 124 um
wavel = cldopts['wavelength'].values

# Actual optical properties dim = [30 wavelength x 61 reff]
e_i = cldopts['extinction_per_mass_crystal']
a_i = cldopts['co_albedo_crystal']
g_i = cldopts['asymmetry_factor_crystal']

if print_opt == True:
   print('min,mean,max mass_ext_crystal: ' + str(np.nanmin(e_i)) + ' ' + str(np.nanmean(e_i)) + ' ' + str(np.nanmax(e_i)))
   print('min,mean,max co_albedo_crystal: ' + str(np.nanmin(a_i)) + ' ' + str(np.nanmean(a_i)) + ' ' + str(np.nanmax(a_i)))
   print('min,mean,max asymmetry_crystal: ' + str(np.nanmin(g_i)) + ' ' + str(np.nanmean(g_i)) + ' ' + str(np.nanmax(g_i)))

# Stack the optical properties into one matrix.
opt = np.stack((e_i,a_i,g_i))   # (3,30,61)
lbl = [r'$\epsilon_i$',r'$\alpha_i$',r'$g_i$']

fs = 13
s = 231

if(print_plot1 == True):
   fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(12,5))
   # Sort the wavelengths from smallest to largest (necessary for the default RRTM file in paritcular)
   iorg = np.argsort(wavel)
   # Generate a colormap based on effective radii.
   mappi = cm.ScalarMappable(colors.Normalize(np.nanmin(r_i),np.nanmax(r_i)),cm.jet)
   yl = [[0,0.5],[0,1],[0,1]]
   for i in np.arange(3):
       for j in np.arange(0,61,3):
           c = mappi.to_rgba(r_i[j])
           ax[i].plot(wavel[iorg],opt[i,iorg,j],color=c)
       #ax[i].set_xscale('log')
       ax[i].set_ylabel(lbl[i],fontsize=fs)
       #ax[i].set_xlim([8,14]) # Infrared atmospheric window
       ax[i].set_xlim([0,100]) # Full atmospheric window
       ax[i].set_ylim(yl[i])
       ax[i].set_xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
       if i == 2:
          mappi.set_array([])
          plt.colorbar(mappi,label=r'$r_i$ [$\mu$m]',ax=ax[i])
   #fig.savefig('cldopt_wavelength.pdf',bbox_inches='tight')
   plt.show()

if(print_plot2 == True):
   fig2 = plt.figure(figsize=(11,7))
   iorgi = np.argsort(r_i)
   iorgl = np.argsort(r_l)
   #farbe = cm.jet(np.linspace(0,1,opt.shape[1]))  #(30, 4)
   mappa = cm.ScalarMappable(colors.Normalize(np.nanmin(np.log10(wavel)),np.nanmax(np.log10(wavel))),cm.jet)
   #mappa = cm.ScalarMappable(colors.Normalize(np.nanmin(wavel),30),cm.jet) # 30 is the second largest value
   for i in np.arange(6):
       ax = fig2.add_subplot(s+i)
       if i < 3:
          for j in np.arange(0,30):
              ax.plot(r_l[iorgl],opt[i,j,iorgl],color=mappa.to_rgba(np.log10(wavel[j]))) # or farbe[j])
          ax.set_xlabel(r'Liquid radius [$\mu$m]',fontsize=fs) 
          ax.set_xlim([2,32])
       if i >= 3:
          for j in np.arange(0,30,3):
              ax.plot(r_i[iorgi],opt[i,j,iorgi],color=mappa.to_rgba(np.log10(wavel[j]))) # farbe[j])
          ax.set_xlabel(r'Ice radius [$\mu$m]',fontsize=fs)
          ax.set_xlim([4,124])
       #ax.set_xscale('log')
       ax.set_ylabel(lbl[i],fontsize=fs)
       if i == 2 or i == 5:
          mappa.set_array([])
          cbar = plt.colorbar(mappa,label=r'log($\lambda$) [$\mu$m]') #,ticks=[-2,0,2])
          #cbar.ax.set_yticklabels(['0.001','0.1','1','10','100'])
#   fig2.savefig('cldopt_reff.pdf',bbox_inches='tight')
   plt.show()

