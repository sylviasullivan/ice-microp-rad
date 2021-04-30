# This is a utility to retain only the whole-second measurements in the StratoClim data.
from netCDF4 import num2date, Dataset
import xarray as xr
import matplotlib.pyplot as plt
import sys, time

basedir = '/work/bb1018/b380873/tropic_vis/obs/'
fi = basedir + 'stratoclim2017.geophysika.0808_1.master.ci_eval.nc'
Stratoclim = Dataset(fi, 'r+')

daten = num2date(times=Stratoclim.variables['time'][:],units='seconds since 2000-01-01 00:00:00 UTC')
# indices to retain associated with whole-second measurements
indx = [i for i, d in enumerate(daten) if d.microsecond == 0]

# recast Stratoclim as an xarray dataset now; Stratoclim2 will hold only whole-second measurements
Stratoclim = xr.open_dataset(fi)
Stratoclim2 = xr.Dataset()

# iterate over the variables in the StratoClim file
for v in Stratoclim.variables:
    Stratoclim2[v] = Stratoclim[v].isel(time=indx)
Stratoclim2.to_netcdf(basedir + 'stratoclim2017.geophysika.0808_1.filtered_per_sec.nc')
