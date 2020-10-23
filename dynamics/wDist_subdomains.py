import pickle,gzip
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys

# Read in the w file.
basedir = '/work/bb1018/b380873/tropic_run2_output/'
wfi = xr.open_dataset(basedir + 'W_3D_icon_tropic_0051_remapdis_global0.025.nc')
w = wfi.w[0]
w_swap = np.transpose(np.array(w),axes=(1,2,0))
lats = wfi.lat
lons = wfi.lon

# Read in the different domain masks.
basedir2 = '/work/bb1018/b380873/SilviaBucci_data/'
f = open(basedir2 + 'MaskCartopy.pkl','rb')
maskk = pickle.load(gzip.open(f))
f.close()
print(maskk['regcode'])
print(maskk['regcode_inv'])
print('~~~~~~~~~~~~~~~~~~~~~~~~~~`')

# Iterate through the subdomains and extract the w values at those points.
for j in np.arange(20):
    kindx = np.argwhere(maskk['mask'] == float(j))
    print(kindx.shape)
    print(kindx[0])
    w_swap_sub = np.array(w_swap[kindx])
    print(w)
sys.exit()

# Extract what is roughly the 500 hPa level = 64 for index 2.
# Pressure levels pulled from here
# https://www.ecmwf.int/en/forecasts/documentation-and-support/correspondence-between-60-and-91-model-level-definitions
w500 = np.array(w[0])
w500 = np.reshape(w500,(1800*4600,))
w500_pos = w500[w500 > 0]

fig = plt.figure()
plt.hist(w500_pos)
#fig.savefig('wDist.pdf',bbox_inches='tight')
plt.show()
