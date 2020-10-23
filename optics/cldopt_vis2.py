import sys
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import xarray as  xr
from matplotlib import cm,colors

# Cloud optical properties files.
print_opt = False
print_plot1 = True
print_plot2 = False

basedir = '/work/bb1018/b380873/rrtm/'
fi = basedir + 'rrtm_cldopt.nc'
cldopts = xr.open_dataset(fi)

# Dimensions -- effective radius and wavelength
r_l = cldopts['re_droplet'].values   #  2 - 32 um
r_i = cldopts['re_crystal'].values   #  4 - 124 um
wavel = cldopts['wavelength'].values

# Actual optical properties dim = [30 wavelength x 61 reff]
e_l = cldopts['extinction_per_mass_droplet']
e_i = cldopts['extinction_per_mass_crystal']
a_l = cldopts['co_albedo_droplet']
a_i = cldopts['co_albedo_crystal']
g_l = cldopts['asymmetry_factor_droplet']
g_i = cldopts['asymmetry_factor_crystal']

if print_opt == True:
   print('min,mean,max mass_ext_droplet: ' + str(np.nanmin(e_l)) + ' ' + str(np.nanmean(e_l)) + ' ' + str(np.nanmax(e_l)))
   print('min,mean,max mass_ext_crystal: ' + str(np.nanmin(e_i)) + ' ' + str(np.nanmean(e_i)) + ' ' + str(np.nanmax(e_i)))
   print('min,mean,max co_albedo_droplet: ' + str(np.nanmin(a_l)) + ' ' + str(np.nanmean(a_l)) + ' ' + str(np.nanmax(a_l)))
   print('min,mean,max co_albedo_crystal: ' + str(np.nanmin(a_i)) + ' ' + str(np.nanmean(a_i)) + ' ' + str(np.nanmax(a_i)))
   print('min,mean,max asymmetry_droplet: ' + str(np.nanmin(g_l)) + ' ' + str(np.nanmean(g_l)) + ' ' + str(np.nanmax(g_l)))
   print('min,mean,max asymmetry_crystal: ' + str(np.nanmin(g_i)) + ' ' + str(np.nanmean(g_i)) + ' ' + str(np.nanmax(g_i)))

# Stack the optical properties into one matrix.
opt = np.stack((e_l,a_l,g_l,e_i,a_i,g_i))   # (6,30,61)
lbl = [r'$\epsilon_l$',r'$\alpha_l$',r'$g_l$',r'$\epsilon_i$',r'$\alpha_i$',r'$g_i$']
# (6, 30, 61)

fs = 13
s = 231

if(print_plot1 == True):
   fig = plt.figure(figsize=(12,5))
   iorg = np.argsort(wavel)
   #farbe = cm.jet(np.linspace(0,1,opt.shape[2]))  #(61, 4)
   mappl = cm.ScalarMappable(colors.Normalize(np.nanmin(r_l),np.nanmax(r_l)),cm.jet)
   mappi = cm.ScalarMappable(colors.Normalize(np.nanmin(r_i),np.nanmax(r_i)),cm.jet)
   yl = [[0,0.5],[0,1],[0,1],[0,0.5],[0,1],[0,1]]
   for i in np.arange(6):
       ax = fig.add_subplot(s+i)
       for j in np.arange(0,61,3):
           if i < 3:
              c = mappl.to_rgba(r_l[j])
           else:
              c = mappi.to_rgba(r_i[j])
           ax.plot(wavel[iorg],opt[i,iorg,j],color=c) #,color=farbe[j])
       #ax.set_xscale('log')
       ax.set_ylabel(lbl[i],fontsize=fs)
       #ax.set_xlim([8,14]) # Infrared atmospheric window
       ax.set_xlim([0,100]) # Full atmospheric window
       ax.set_ylim(yl[i])
       if i > 2:
          ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
       if i == 2:
          mappl.set_array([])
          plt.colorbar(mappl,label=r'$r_l$ [$\mu$m]')
       if i == 5:
          mappi.set_array([])
          plt.colorbar(mappi,label=r'$r_i$ [$\mu$m]')
       print(opt[i,:,30])
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

