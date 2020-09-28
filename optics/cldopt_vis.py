import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Cloud optical properties files.
basedir = '/work/bb1131/b380873/rrtm/'
fi = basedir + 'rrtm_cldopt.nc'
cldopts = xr.open_dataset(fi)

# Dimensions -- effective radius and wavelength
r_l = cldopts['re_droplet']
r_i = cldopts['re_crystal']
wavel = cldopts['wavelength']
print(wavel)

# Actual optical properties dim = [30 wavelength x 61 reff]
e_l = cldopts['extinction_per_mass_droplet']
e_i = cldopts['extinction_per_mass_crystal']

# Make matrices for surface plotting
wavel_mat = np.matlib.repmat(wavel,r_l.shape[0],1).T
r_l_mat = np.matlib.repmat(r_l,wavel.shape[0],1)
r_i_mat = np.matlib.repmat(r_i,wavel.shape[0],1)

fs = 13
fig = plt.figure()
ax = fig.add_subplot(121,projection='3d')
ax.plot_surface(np.log10(wavel_mat),np.log10(r_l_mat),e_l,cmap=cm.viridis)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
ax.set_ylabel(r'$r_{eff,l}$ [$\mu$m]',fontsize=fs)
ax.set_zlabel('$\epsilon_l$',fontsize=fs)

ax = fig.add_subplot(122,projection='3d')
ax.plot_surface(np.log10(wavel_mat),np.log10(r_i_mat),e_i,cmap=cm.viridis)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
ax.set_ylabel(r'$r_{eff,i}$ [$\mu$m]',fontsize=fs)
ax.set_zlabel('$\epsilon_i$',fontsize=fs)

plt.show()

